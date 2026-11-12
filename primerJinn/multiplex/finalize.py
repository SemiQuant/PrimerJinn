"""Final checks on a chosen pool: coverage, unintended products, dimers,
annealing temperature advice and amplicon secondary structure."""
from __future__ import annotations

import itertools
import math
import multiprocessing as mp
import re
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List

import primer3

from . import candidates as cand
from . import thermo
from .targets import coverage_gaps


def _ctx():
    return mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()


def primer_names(targets, chosen) -> Dict[str, Dict[str, str]]:
    out = {}
    for t in targets:
        named = [s.name for s in t.sites if s.name]
        label = named[0] if named else ""
        if len(named) > 1:
            label += f"+{len(named) - 1}"
        label = re.sub(r"[^A-Za-z0-9_.+-]", "_", label)[:24]
        base = f"{t.id}_{label}" if label else t.id
        out[t.id] = {"F": f"{base}_F", "R": f"{base}_R"}
    return out


# ----------------------------------------------------------------------------
# annealing temperature from structure-aware binding curves
# ----------------------------------------------------------------------------
def _site_open_curve(args):
    """dG_open of one primer site in genomic context at each temperature."""
    context, end1_on_strand, length, temps, cfg = args
    out = []
    for T in temps:
        up = cand._plfold_up(context, T, length, cfg)
        out.append(cand._dg_open(up, end1_on_strand, length, T))
    return out


def _site_context(p: cand.Primer, store, pad: int = 300):
    """(context strand string, 1-based end of the site window on that strand)."""
    region, rs = store.region(p.seqid, p.start - pad, p.end + pad)
    n = len(region)
    if p.side == "F":
        bottom = cand.revcomp(region)
        a = p.start - rs
        return bottom, n - a
    return region, p.end - rs + 1


def annealing_advice(primers: Dict[str, cand.Primer], cfg, store) -> dict:
    """Vendor-rule Ta and the equilibrium (model) Ta with structure.

    model Ta = highest temperature on a 1 C grid where every primer keeps
    theta >= 0.9 (dG_dup + hairpin cost + site opening cost, primer in
    excess at cfg.primer_nM), minus 1 C."""
    tm_min = min(p.tm for p in primers.values())
    tm_neb_min = min(p.tm_neb for p in primers.values())
    vendor = None
    if cfg.polymerase == "q5":
        # NEB rule at the REAL per-primer concentration (the 500 nM scale is only the
        # design reference); NEB's calculator has no DMSO term, so subtract it here
        tm_rxn = min(thermo.tm_neb(p.seq, cfg.primer_nM) for p in primers.values())
        vendor = min(round(tm_rxn + 1.0 - cfg.dmso_pct * 0.6, 1), 72.0)
    elif cfg.polymerase == "kapa":
        # KAPA2G Robust TDS: 55 C for primer Tm < 55, 60 C for 55-65, 65 C above 65
        vendor = 55.0 if tm_min < 55 else (60.0 if tm_min <= 65 else 65.0)
    lo = int(max(35, math.floor(min(cfg.ta, vendor or cfg.ta) - 24)))
    hi = int(math.ceil(max(cfg.ta, vendor or cfg.ta) + 4))
    grid = list(range(lo, hi + 1))
    coarse = grid[::2] + ([grid[-1]] if grid[-1] not in grid[::2] else [])
    names = list(primers)
    jobs = []
    for nm in names:
        ctx, end1 = _site_context(primers[nm], store)
        jobs.append((ctx, end1, primers[nm].length, coarse, cfg))
    if cfg.threads > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(min(cfg.threads, len(jobs)), mp_context=_ctx()) as ex:
            curves = list(ex.map(_site_open_curve, jobs))
    else:
        curves = [_site_open_curve(j) for j in jobs]
    phys = cfg.phys()
    table = {}
    for nm, curve in zip(names, curves):
        p = primers[nm]
        r = thermo.analyze(p.seq, phys)
        row = {}
        for T in grid:
            # linear interpolation of the coarse site-opening curve
            k = max(i for i, c in enumerate(coarse) if c <= T) if T >= coarse[0] else 0
            if k + 1 < len(coarse) and coarse[k] != T:
                f = (T - coarse[k]) / (coarse[k + 1] - coarse[k])
                so = curve[k] * (1 - f) + curve[k + 1] * f
            else:
                so = curve[k]
            hp = primer3.calc_hairpin(p.seq, **cfg.thal_kw(T))
            hc = thermo.structure_cost(hp.dg / 1000.0, T) if hp.structure_found else 0.0
            dg = r.dG(T) + hc + so
            row[T] = dict(dg_dup=r.dG(T), hairpin=hc, site_open=so, dg_eff=dg,
                          theta=thermo.theta_from_dg(dg, T, cfg.primer_nM),
                          theta_nostruct=thermo.theta_from_dg(r.dG(T), T, cfg.primer_nM))
        table[nm] = row

    def t_at(level, key="theta"):
        ok = [T for T in grid if min(table[nm][T][key] for nm in names) >= level]
        return max(ok) if ok else None

    t90 = t_at(0.9)
    t50 = t_at(0.5)
    t90_ns = t_at(0.9, "theta_nostruct")
    # Model Ta: every primer >= 90% bound on an unfolded template (the classical
    # Tm - 4 to 5 C), capped where template and primer structure leave the weakest
    # primer under half bound. Requiring 90% binding WITH structure is too strict:
    # GC-rich templates fold more as T falls, which drives the answer far below
    # vendor practice (54 C vs KAPA's 60 C on the M. tuberculosis example).
    cands = [x for x in (t90_ns, t50) if x is not None]
    model = min(cands) if cands else None
    lim = None
    if model is not None:
        nxt = model + 1 if model + 1 in grid else model
        key = "theta" if (t50 is not None and model == t50) else "theta_nostruct"
        lim = min(names, key=lambda nm: table[nm][nxt][key])
    else:
        lim = min(names, key=lambda nm: table[nm][grid[0]]["theta"])
    if cfg.polymerase == "q5":
        recommended = vendor
    else:
        recommended = model if model is not None else (vendor if vendor is not None else cfg.ta)
    rec_grid = min(grid, key=lambda T: abs(T - recommended))
    per_primer = {nm: table[nm][rec_grid] for nm in names}
    spread = (max(v["dg_eff"] for v in per_primer.values()) - min(v["dg_eff"] for v in per_primer.values()))
    return dict(rec_grid=rec_grid, vendor_ta=vendor, model_ta=model, t50=t50, t90=t90, t90_without_structure=t90_ns,
                limiting_primer=lim, recommended=recommended, at_recommended=per_primer,
                dg_eff_spread=spread, tm_min=tm_min, tm_max=max(p.tm for p in primers.values()),
                tm_neb_min=tm_neb_min, tm_neb_max=max(p.tm_neb for p in primers.values()),
                table=table, grid=grid)


# ----------------------------------------------------------------------------
# amplicon structure reports
# ----------------------------------------------------------------------------
def _amp_report(args):
    from . import structure

    amp, lf, lr, ta, cond = args
    try:
        return structure.amplicon_report(amp, lf, lr, temps=(ta, 72.0), ta=ta, cond=cond)
    except Exception as exc:  # one failed report must not end a finished design
        return exc


def finish(targets, chosen, pools, sites, cfg, store, spec_data, dimer_cache, offtarget_by_seq, verbose):
    warnings: List[str] = []
    names = primer_names(targets, chosen)

    # coverage of every requested base by the bases between primers
    inserts = [(pr.fwd.seqid, *pr.insert) for pr in chosen.values()]
    cov = coverage_gaps(sites, inserts)
    for s, gaps in cov:
        warnings.append(f"NOT COVERED: {s.label()} gaps {gaps}")

    # primers by name
    primers: Dict[str, cand.Primer] = {}
    pool_of: Dict[str, int] = {}
    for tid, pr in chosen.items():
        primers[names[tid]["F"]] = pr.fwd
        primers[names[tid]["R"]] = pr.rev
        pool_of[names[tid]["F"]] = pool_of[names[tid]["R"]] = pools[tid]

    # dimers within each pool, every pair including self
    dimers = []
    keyed = {nm: ((p.seq, p.side), cand.tailed(p, cfg)) for nm, p in primers.items()}
    for a, b in itertools.combinations_with_replacement(sorted(primers), 2):
        if pool_of[a] != pool_of[b]:
            continue
        ka, kb = keyed[a][0], keyed[b][0]
        k = (ka, kb) if ka <= kb else (kb, ka)
        if k in dimer_cache:
            d3, dany = dimer_cache[k]
        else:
            d3, dany = cand.thal_pair(keyed[a][1], keyed[b][1], cfg)
        cost = cand.dimer_cost(d3, dany, cfg)
        dimers.append(dict(pool=pool_of[a], primer1=a, primer2=b, d3_kcal=round(d3, 2), dany_kcal=round(dany, 2),
                           call=("hard" if cost >= 1000 else "soft" if cost > 0 else "ok")))
    dimers.sort(key=lambda d: (d["d3_kcal"], d["dany_kcal"]))
    for d in dimers:
        if d["call"] == "hard":
            warnings.append(f"strong dimer in pool {d['pool']}: {d['primer1']} / {d['primer2']} "
                            f"(3' dG {d['d3_kcal']}, any {d['dany_kcal']} kcal/mol at 37 C)")

    # unintended products within pools, and off-target sites of chosen primers
    products: List[dict] = []
    offtarget: Dict[str, List[dict]] = {}
    if spec_data is not None:
        from . import specificity as spec

        seqs, index, intended, scored, params, g = spec_data
        idx_names: Dict[int, List[str]] = {}
        for nm, p in primers.items():
            idx_names.setdefault(index[p.seq], []).append(nm)
        chosen_idx = set(idx_names)
        sub = [s for s in scored if s["primer"] in chosen_idx]
        prods = spec.enumerate_products(sub, seqs, g, params, names=[str(i) for i in range(len(seqs))])
        from .design import classify_product

        pair_of, tid_target = {}, {}
        tmap = {t.id: t for t in targets}
        for tid, pr in chosen.items():
            pair_of[names[tid]["F"]] = pair_of[names[tid]["R"]] = pr
            tid_target[names[tid]["F"]] = tid_target[names[tid]["R"]] = tmap[tid]
        for d in prods:
            i, j = int(d["fwd"]), int(d["rev"])
            seen_pools = set()
            for ni in idx_names[i]:
                for nj in idx_names[j]:
                    pl = pool_of[ni]
                    if pool_of[nj] != pl or pl in seen_pools:
                        continue
                    kind = classify_product(d, i, j, pair_of[ni], pair_of[nj], seqs, d["contig"])
                    if kind == "designed":
                        continue
                    seen_pools.add(pl)
                    products.append(dict(pool=pl, kind={"repeat": "repeated locus", "outer": "outer product",
                                                        "other": "unintended"}[kind],
                                         fwd=ni, rev=nj, contig=d["contig"], start=d["start"] + 1,
                                         end=d["end"], size=d["size"], fwd_mismatches=d["fwd_mm"],
                                         rev_mismatches=d["rev_mm"], score=round(d["score"], 4), call=d["call"]))
        for nm, p in primers.items():
            rows = []
            for s in offtarget_by_seq.get(index[p.seq], []):
                rows.append(dict(primer=nm, contig=g.names[s["contig"]], start=s["start"] + 1,
                                 strand="+" if s["strand"] == 1 else "-", mismatches=s["mm"],
                                 mismatch_positions_from_3prime=s["mm_pos"], dg_kcal=s["dg_any"],
                                 score=round(s["score"], 4)))
            if rows:
                offtarget[nm] = sorted(rows, key=lambda r: -r["score"])
        for d in products:
            if d["call"] != "likely":
                continue
            if d["kind"] == "repeated locus":
                warnings.append(f"{d['fwd']} / {d['rev']} also amplify a second copy of the same locus at "
                                f"{d['contig']}:{d['start']}-{d['end']} ({d['size']} bp, {d['fwd_mismatches']}+"
                                f"{d['rev_mismatches']} mismatches): reads from the copies cannot be told apart")
            elif d["kind"] == "outer product":
                same_range = bool(tid_target[d["fwd"]].origins & tid_target[d["rev"]].origins)
                if d["size"] < cfg.outer_min and not same_range:
                    warnings.append(f"outer product of {d['size']} bp in pool {d['pool']} ({d['fwd']} + {d['rev']})")
            else:
                warnings.append(f"unintended product in pool {d['pool']}: {d['fwd']} + {d['rev']} "
                                f"{d['contig']}:{d['start']}-{d['end']} ({d['size']} bp, score {d['score']})")

    # annealing temperature
    ta = annealing_advice(primers, cfg, store)

    # amplicon structure at the recommended Ta and at 72 C
    cond = cfg.struct_cond()
    jobs, order = [], []
    for tid, pr in chosen.items():
        # gene-specific amplicon: the tails are the same on every amplicon and their
        # panhandle would mask the template's own structure
        amp = store.region(pr.fwd.seqid, pr.fwd.start, pr.rev.end)[0]
        jobs.append((amp, pr.fwd.length, pr.rev.length, float(ta["recommended"]), cond))
        order.append(tid)
    if cfg.threads > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(min(cfg.threads, len(jobs)), mp_context=_ctx()) as ex:
            reps = list(ex.map(_amp_report, jobs))
    else:
        reps = [_amp_report(j) for j in jobs]
    amplicons = {}
    for tid, rep in zip(order, reps):
        if isinstance(rep, Exception):
            warnings.append(f"structure report for {tid} failed: {rep!r}")
        else:
            amplicons[tid] = rep
    reps = list(amplicons.values()) or []

    # reaction advice from amplicon composition and structure
    advice = []
    gcs = [r.gc for r in reps] or [pr.gc for pr in chosen.values()]
    tms = [r.duplex_tm for r in reps] or [0.0]
    from .config import POLYMERASES

    gc_rich = max(gcs) >= 0.65 or max(tms) > 92.0 or max([r.max_gc_50 for r in reps] or [0]) >= 0.80
    longest = max(pr.length for pr in chosen.values()) + len(cfg.tail_fwd) + len(cfg.tail_rev)
    if cfg.polymerase == "kapa":
        initial = "5 min" if max(gcs) > 0.70 else "3 min"
        advice.append(f"Cycling (KAPA2G Robust data sheet): 95 C {initial}; then 95 C 15 s, anneal at the "
                      f"recommended Ta, 72 C 15 s (amplicons < 1 kb); final 72 C 1 min. KAPA advise a 15 s "
                      f"anneal because long anneals give non-specific products with this enzyme; multiplex pools "
                      f"often need longer, so compare 15 s with 60 s on the first pool.")
        if gc_rich:
            advice.append(f"GC-rich amplicons (max GC {100 * max(gcs):.0f}%, max product Tm {max(tms):.1f} C): KAPA "
                          f"recommend 5% DMSO for amplicons above 65% GC. DMSO lowers primer Tm by about 0.6 C "
                          f"per %; re-run with --dmso 5 so primers and Ta are scored with it.")
    elif cfg.polymerase == "q5":
        ext = max(10, int(math.ceil(longest / 1000 * 20 / 5.0)) * 5)
        advice.append(f"Cycling (NEB Q5): 98 C 30 s; then 98 C 10 s, anneal at the recommended Ta 20 s, 72 C {ext} s; "
                      f"final 72 C 2 min. Keep the extension short to suppress longer outer products.")
    else:
        advice.append("Cycling (Taq-type): 95 C 2 min; then 95 C 15 s, anneal at the recommended Ta 30 s, 72 C "
                      "20 to 30 s; final 72 C 2 min. Keep the extension short to suppress longer outer products.")
    if gc_rich and cfg.polymerase != "kapa":
        den = "98 C" if POLYMERASES[cfg.polymerase]["proofreading"] else "95 C"
        advice.append(f"GC-rich amplicons (max GC {100 * max(gcs):.0f}%, max product Tm {max(tms):.1f} C): denature "
                      f"at {den} with a longer initial step and consider 3 to 5% DMSO or 1 M betaine. DMSO lowers "
                      f"primer Tm by about 0.6 C per %; re-run with --dmso to rescore.")
    if any(r.g4 for r in reps):
        n = sum(1 for r in reps if r.g4)
        advice.append(f"{n} amplicon(s) carry G-quadruplex motifs: 7-deaza-dGTP can help if they drop out.")
    if min(gcs) < 0.30:
        advice.append("AT-rich amplicons: a 60 to 65 C extension and extra MgCl2 improve their yield.")
    if cfg.polymerase == "q5":
        advice.append(f"Annealing: {ta['vendor_ta']:.1f} C by the NEB Q5 rule (lowest primer Tm + 1 C). On the first "
                      f"pool run a gradient from {ta['vendor_ta'] - 3:.0f} to {ta['vendor_ta'] + 3:.0f} C. The "
                      f"equilibrium binding numbers ignore the polymerase's DNA-binding domain, so for Q5 they rank "
                      f"primers rather than set Ta.")
    elif ta["model_ta"] is not None:
        why = ("template structure leaves the weakest primer at half binding"
               if ta["t50"] is not None and ta["model_ta"] == ta["t50"] and ta["t50"] < (ta["t90_without_structure"] or 999)
               else "every primer is at least 90% bound on unfolded template")
        msg = (f"Annealing: {ta['model_ta']} C, where {why} (limiting primer {ta['limiting_primer']}). "
               f"Run a gradient from {ta['model_ta'] - 3} to {ta['model_ta'] + 3} C.")
        if ta["vendor_ta"] is not None:
            msg += f" The KAPA data sheet rule for this primer Tm gives {ta['vendor_ta']:g} C."
        advice.append(msg)
    # structure at the primer sites at the recommended temperature
    rec_t = min(ta["grid"], key=lambda T: abs(T - ta["recommended"]))
    blocked = sorted((nm for nm, row in ta["table"].items() if row[rec_t]["site_open"] > cfg.site_open_reject),
                     key=lambda nm: -ta["table"][nm][rec_t]["site_open"])
    if blocked:
        shown = ", ".join(f"{nm} ({ta['table'][nm][rec_t]['site_open']:.1f} kcal/mol)" for nm in blocked[:5])
        shown += f" and {len(blocked) - 5} more" if len(blocked) > 5 else ""
        better = None
        for dT in range(1, 6):
            T = rec_t + dT
            if T not in ta["grid"]:
                break
            if all(ta["table"][nm][T]["theta"] > ta["table"][nm][rec_t]["theta"] for nm in blocked) and \
                    all(row[T]["theta"] >= min(0.5, row[rec_t]["theta"]) for row in ta["table"].values()):
                better = T
                break
        msg = f"Template structure covers the binding site of {shown} at {rec_t} C."
        if better:
            msg += f" Their net binding improves at {better} C without weakening the rest of the pool."
        else:
            msg += (" Raising Ta weakens binding more than it opens the structure; if these amplicons are weak, "
                    "add 5% DMSO or 1 M betaine, or redesign those primers.")
        advice.append(msg)
    ta["advice"] = advice
    return dict(coverage=cov, products=products, dimers=dimers, ta=ta, amplicons=amplicons,
                offtarget=offtarget, warnings=warnings, names=names)
