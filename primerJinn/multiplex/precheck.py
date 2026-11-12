"""Input precheck: is everything there, readable, and consistent with the genome?

Run before any design (and alone with --check). Each finding is one line with a level:
OK, INFO, WARN (design runs but may not do what you expect) or ERROR (design stops).

Checks
  environment  numba (genome-wide specificity), ViennaRNA, primer3, pysam
  genome       resolves; size, sequences, GC; resolver notes and warnings
  positions    file format; every unreadable row (not just the first); unknown sequence
               names; positions beyond a sequence end; renamed chromosomes; VCF REF
               alleles against the genome bases and VCF ##contig lengths against the
               genome (a different assembly or version shows up here); BED-style
               coordinates in a non-BED table; positions in N gaps or ambiguity codes;
               too little sequence left at a sequence end for an amplicon; grouping
               into amplicon targets
  variants     readable; sequence names; VCF REF alleles; variants at requested positions
  settings     the design settings are valid (Config.resolve)
  output       the output folder can be created and written
"""
from __future__ import annotations

import copy
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

LEVELS = ("OK", "INFO", "WARN", "ERROR")
SHOW_MAX = 15            # findings of one kind listed before "... and N more"


@dataclass
class Item:
    level: str
    topic: str
    msg: str


@dataclass
class Precheck:
    items: List[Item] = field(default_factory=list)
    meta: Optional[dict] = None
    rows: list = field(default_factory=list)       # parsed positions (dicts)
    sites: list = field(default_factory=list)      # targets.Site
    targets: list = field(default_factory=list)    # targets.Target (multiplex grouping)
    variants: object = None                        # variants.VariantMask
    variant_notes: List[str] = field(default_factory=list)

    def add(self, level: str, topic: str, msg: str) -> None:
        self.items.append(Item(level, topic, msg))

    def add_many(self, level: str, topic: str, msgs: List[str], what: str) -> None:
        for m in msgs[:SHOW_MAX]:
            self.add(level, topic, m)
        if len(msgs) > SHOW_MAX:
            self.add(level, topic, f"... and {len(msgs) - SHOW_MAX} more {what}")

    @property
    def errors(self) -> List[Item]:
        return [i for i in self.items if i.level == "ERROR"]

    @property
    def warnings(self) -> List[Item]:
        return [i for i in self.items if i.level == "WARN"]

    @property
    def ok(self) -> bool:
        return not self.errors

    def text(self, only_problems: bool = False) -> str:
        lines = ["Input check"]
        for i in self.items:
            if only_problems and i.level not in ("WARN", "ERROR"):
                continue
            lines.append(f"  {i.level:<5}  {i.topic:<9}  {i.msg}")
        ne, nw = len(self.errors), len(self.warnings)
        if ne:
            lines.append(f"Result: {ne} error(s), {nw} warning(s). Fix the errors before designing.")
        elif nw:
            lines.append(f"Result: no errors, {nw} warning(s).")
        else:
            lines.append("Result: all checks passed.")
        return "\n".join(lines)


# ----------------------------------------------------------------------------
def positions_rows(source, meta, errors=None, **kw):
    """Parse --positions (a file or an inline list such as '7570,761155') against the genome."""
    from . import genome as gmod

    if isinstance(source, str) and not os.path.exists(source) \
            and not re.search(r"\.(txt|tsv|csv|xlsx|xls|bed|vcf)(\.b?gz)?$", source, re.I):
        if re.search(r"(?<![\d,])\d{1,3}(?:,\d{3})+(?![\d,])", source) and \
                re.search(r":\d{1,3},\d{3}\b|(?<![\d,])\d{1,3}(?:,\d{3}){2,}(?![\d,])", source):
            raise ValueError("thousands separators are not supported in positions; write 7570, not 7,570 "
                             "(separate several positions with spaces or commas)")
        source = [x for x in re.split(r"[\s,;]+", source.strip()) if x]   # inline list
    return gmod.parse_positions(source, meta, errors=errors, **kw)


def _describe_source(source) -> str:
    from . import genome as gmod

    if isinstance(source, str) and os.path.isfile(source):
        low = source.lower()
        if gmod.is_vcf(source):
            kind = "VCF (one target per record, covering its REF allele)"
        elif low.endswith((".bed", ".bed.gz")):
            kind = "BED (0-based half-open, converted to 1-based)"
        elif low.endswith((".xlsx", ".xls")):
            kind = "Excel table"
        else:
            kind = "text table"
        return f"{source}: {kind}"
    return "inline list"


def _ref_check(pc: Precheck, topic: str, recs, meta) -> None:
    """Compare VCF REF alleles with the genome. recs: dicts with seqid, start, ref, line."""
    from . import genome as gmod

    checked, bad = 0, []
    for r in recs:
        ref = r.get("ref") or ""
        if not ref or not re.fullmatch(r"[ACGT]+", ref):
            continue                       # N, symbolic or missing REF: nothing to compare
        try:
            g = gmod.fetch(meta, r["seqid"], r["start"], r["start"] + len(ref) - 1)
        except ValueError:
            continue                       # out of range is reported elsewhere
        checked += 1
        if g != ref:
            bad.append(f"VCF line {r.get('line', '?')}: REF {ref} but the genome has {g} at "
                       f"{r['seqid']}:{r['start']}")
    if not checked:
        return
    if not bad:
        pc.add("OK", topic, f"REF alleles match the genome at all {checked} checked record(s)")
    elif len(bad) > 0.5 * checked:
        pc.add("ERROR", topic, f"REF alleles differ from the genome at {len(bad)} of {checked} records: the VCF "
                               f"was probably made on a different reference (assembly or version) or with "
                               f"shifted coordinates")
        pc.add_many("ERROR", topic, bad, "mismatches")
    else:
        pc.add("WARN", topic, f"REF alleles differ from the genome at {len(bad)} of {checked} records "
                              f"(the rest match): check these records")
        pc.add_many("WARN", topic, bad, "mismatches")


def _group_unknown(errs: List[str], meta) -> List[str]:
    """One line per unknown sequence name instead of one per row, with a hint."""
    unknown, rest = {}, []
    for e in errs:
        m = re.match(r"(VCF line \d+|Position row \d+) .*?unknown sequence '([^']+)'", e)
        if m:
            unknown.setdefault(m.group(2), []).append(m.group(1))
        else:
            rest.append(e)
    out = []
    for name, where in unknown.items():
        hint = ""
        if meta:
            ids = [x["seqid"] for x in meta["sequences"]]
            base = name.split(".")[0].lower()
            same = [i for i in ids if i.split(".")[0].lower() == base]
            if same:
                hint = f" The genome has {same[0]}: a different version of the same sequence, so positions may differ."
            elif len(ids) == 1:
                hint = f" The genome has one sequence, {ids[0]}; if this is the same sequence, rename it in the file."
            else:
                hint = f" Known: {', '.join(ids[:10])}{' ...' if len(ids) > 10 else ''}."
        locs = ", ".join(where[:5]) + (f" and {len(where) - 5} more" if len(where) > 5 else "")
        out.append(f"unknown sequence '{name}' ({len(where)} row(s): {locs}).{hint}")
    return out + rest


def _contig_check(pc: Precheck, topic: str, path, meta) -> None:
    from . import genome as gmod

    try:
        recs, hdr = gmod.read_vcf(path, [])
    except (ValueError, OSError):
        return
    exact = {}
    for x in meta["sequences"]:                     # exact ids and aliases, no fuzzy chromosome matching
        for k in [x["seqid"], *x.get("aliases", [])]:
            exact.setdefault(k.lower(), x)
    used = {r["chrom"] for r in recs}
    for cid, ln in hdr["contigs"].items():
        level = "ERROR" if cid in used else "INFO"
        note = "" if cid in used else " (no record uses it)"
        s = exact.get(cid.lower())
        if s is None:
            same = [x["seqid"] for x in meta["sequences"]
                    if x["seqid"].split(".")[0].lower() == cid.split(".")[0].lower()]
            if same:
                pc.add(level, topic, f"VCF header contig {cid} is a different version of the genome's {same[0]}"
                                     f"{note}: positions and alleles may not match")
            continue
        if ln and ln != s["length"]:
            pc.add(level, topic, f"VCF header says contig {cid} is {ln:,} bp but the genome's {s['seqid']} is "
                                 f"{s['length']:,} bp{note}: different assembly or version")
    if hdr["reference"]:
        pc.add("INFO", topic, f"VCF ##reference: {hdr['reference']}")


def _base_checks(pc: Precheck, rows, meta, amp_max: int) -> None:
    """N gaps, ambiguity codes, and room for an amplicon at sequence ends."""
    from . import genome as gmod

    look = gmod._seq_lookup(meta)
    in_n, ambig, ends, nflank = [], [], [], []
    half = amp_max // 2
    for r in rows:
        s = look.get(r["seqid"].lower())
        if s is None:
            continue
        L = s["length"]
        try:
            core = gmod.fetch(meta, r["seqid"], r["start"], r["end"])
            a, b = max(1, r["start"] - amp_max), min(L, r["end"] + amp_max)
            flank = gmod.fetch(meta, r["seqid"], a, b)
        except ValueError:
            continue
        lab = f"{r['name']} ({r['seqid']}:{r['start']}{'' if r['end'] == r['start'] else '-' + str(r['end'])})"
        if "N" in core:
            in_n.append(f"{lab} lies in an N gap ({core.count('N')} of {len(core)} bases are N)")
        elif re.search(r"[^ACGT]", core):
            ambig.append(f"{lab} has ambiguity codes in the genome ({core[:30]})")
        nbad = len(re.findall(r"[^ACGT]", flank))
        if nbad > 0.1 * len(flank):
            nflank.append(f"{lab}: {100 * nbad / len(flank):.0f}% N or ambiguous bases within {amp_max} bp; "
                          f"primers may not fit")
        if r["start"] - 1 < half or L - r["end"] < half:
            ends.append(f"{lab} is within {half} bp of the end of {r['seqid']} ({L:,} bp); an amplicon may not "
                        f"fit (circular sequences are not wrapped)")
    pc.add_many("WARN", "positions", in_n, "positions in N gaps")
    pc.add_many("WARN", "positions", ambig, "positions on ambiguity codes")
    pc.add_many("WARN", "positions", nflank, "positions with N-rich flanks")
    pc.add_many("WARN", "positions", ends, "positions near sequence ends")


def _env_checks(pc: Precheck, specificity: bool) -> None:
    for mod, label, need in (("RNA", "ViennaRNA", True), ("primer3", "primer3-py", True), ("pysam", "pysam", True)):
        try:
            m = __import__(mod)
            pc.add("OK", "software", f"{label} {getattr(m, '__version__', '')}".strip())
        except Exception as exc:   # noqa: BLE001 - any import failure is a broken install
            if need:
                pc.add("ERROR", "software", f"{label} could not be loaded ({exc}); install it, for example with "
                                            f"the conda environment in environment.yml")
    if specificity:
        try:
            from . import specificity as _spec  # noqa: F401
            import numba
            pc.add("OK", "software", f"numba {numba.__version__}")
        except ImportError as exc:
            pc.add("ERROR", "software", str(exc).replace("\n", " "))


def _out_check(pc: Precheck, out: str, must_exist: bool = False) -> None:
    """must_exist: getMultiPrimerSet writes <out>.xlsx and does not create folders."""
    d = Path(os.path.abspath(out)).parent
    if must_exist and not d.is_dir():
        pc.add("ERROR", "output", f"the output folder {d} does not exist; create it first")
        return
    probe = d
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    if not probe.is_dir():
        pc.add("ERROR", "output", f"cannot create the output folder {d}: {probe} is not a folder")
    elif not os.access(probe, os.W_OK):
        pc.add("ERROR", "output", f"cannot write results to {d}: no write permission in {probe}")
    else:
        pc.add("OK", "output", f"results go to {os.path.abspath(out)}.xlsx" if must_exist
               else f"results go to {os.path.abspath(out)}_* and {os.path.abspath(out)}.xlsx")


# ----------------------------------------------------------------------------
def run(genome, positions, cfg=None, variants: Optional[str] = None, variant_min_af: float = 0.0,
        genome_dir: Optional[str] = None, out: Optional[str] = None, check_env: bool = True,
        regions_mode: bool = False, product_max: Optional[int] = None) -> Precheck:
    """Check all inputs. cfg: an unresolved multiplex Config (None = defaults).
    regions_mode: checks for getMultiPrimerSet (one-sequence FASTA, one product per region,
    overlapping regions not allowed); product_max is its largest product."""
    from . import genome as gmod
    from .config import Config, primer_length_range
    from .targets import Site, group_sites, span_limits

    pc = Precheck()
    cfg = copy.deepcopy(cfg) if cfg is not None else Config()
    if check_env:
        _env_checks(pc, cfg.specificity and not regions_mode)

    # genome ------------------------------------------------------------
    try:
        meta = gmod.resolve_genome(genome, genome_dir)
    except gmod.GenomeResolutionError as exc:
        pc.add("ERROR", "genome", str(exc))
        meta = None
    except Exception as exc:   # noqa: BLE001 - unreadable FASTA, network failure
        pc.add("ERROR", "genome", f"could not read or fetch the genome: {exc}")
        meta = None
    if meta:
        pc.meta = meta
        what = meta.get("organism") or meta.get("query")
        acc = f" ({meta['accession']})" if meta.get("accession") else ""
        pc.add("OK", "genome", f"{what}{acc}: {meta['n_sequences']} sequence(s), {meta['total_length']:,} bp, "
                               f"GC {meta['gc_percent']}%")
        for n in meta.get("notes", []):
            pc.add("INFO", "genome", n)
        for w in meta.get("warnings", []):
            pc.add("INFO" if "bare integer positions" in w else "WARN", "genome", w)
        if regions_mode and meta["n_sequences"] != 1:
            pc.add("ERROR", "genome", f"getMultiPrimerSet needs a FASTA with one sequence; this one has "
                                      f"{meta['n_sequences']}. Use primerjinn-multiplex for multi-sequence genomes.")

    # positions -----------------------------------------------------------
    pc.add("INFO", "positions", f"input: {_describe_source(positions)}")
    errs: List[str] = []
    try:
        rows = positions_rows(positions, meta, errors=errs, single_sequence_fallback=regions_mode)
    except (ValueError, FileNotFoundError, OSError, EOFError) as exc:
        pc.add("ERROR", "positions", str(exc))
        rows = []
    pc.add_many("ERROR", "positions", _group_unknown(errs, meta), "unreadable rows")
    if meta is None:
        return pc
    pc.rows = rows
    is_vcf_file = isinstance(positions, str) and os.path.isfile(positions) and gmod.is_vcf(positions)
    if is_vcf_file:
        _contig_check(pc, "positions", positions, meta)
        try:
            nblock = sum(r["ref_block"] for r in gmod.read_vcf(positions, [])[0])
        except (ValueError, OSError):
            nblock = 0
        if nblock:
            pc.add("INFO", "positions", f"{nblock} gVCF reference block(s) or no-call record(s) skipped "
                                        f"(no alternative allele)")
    if not rows:
        if not any(i.topic == "positions" and i.level == "ERROR" for i in pc.items):
            pc.add("ERROR", "positions", "no positions found in the input")
    else:
        nseq = len({r["seqid"] for r in rows})
        pc.add("OK", "positions", f"{len(rows)} position(s) read on {nseq} sequence(s)"
                                  + (f"; {len(errs)} row(s) could not be used" if errs else ""))
        forced = sorted({r["seq_in"] for r in rows if r.get("forced")})
        if forced:
            pc.add("WARN", "positions", f"sequence name(s) {', '.join(forced[:5])} are not in the FASTA; taken as its "
                                        f"only sequence, {rows[0]['seqid']}")
        renamed = sorted({(r["seq_in"], r["seqid"]) for r in rows if r.get("seq_in") and not r.get("forced")})
        if renamed:
            pc.add("INFO", "positions", "sequence names mapped to the genome: "
                   + ", ".join(f"{a} -> {b}" for a, b in renamed[:8]) + (" ..." if len(renamed) > 8 else ""))
        if is_vcf_file:
            _ref_check(pc, "positions", rows, meta)
        elif (isinstance(positions, str) and os.path.isfile(positions)
              and not positions.lower().endswith((".bed", ".bed.gz"))):
            ranges = [r for r in rows if r["end"] > r["start"]]
            if len(ranges) >= 3 and len(ranges) == len(rows) and all(r["end"] == r["start"] + 1 for r in rows):
                pc.add("WARN", "positions", "every range is exactly 2 bp (end = start + 1). If these are 0-based, "
                                            "BED-style coordinates, save the file with a .bed extension so they are "
                                            "converted; otherwise each target is 2 bases")
        _base_checks(pc, rows, meta, product_max or cfg.amp_max)
        pc.sites = [Site(r["seqid"], r["start"], r["end"], r["name"] or "") for r in rows]

    # settings and grouping -------------------------------------------------
    gc = (meta.get("gc_percent") or 50) / 100
    if not regions_mode:
        rc = None
        try:
            rc = copy.deepcopy(cfg).resolve(gc)
            pc.add("OK", "settings", f"{rc.polymerase} preset, target primer Tm {rc.tm:.1f} C ({rc.tm_scale} scale), "
                                     f"design Ta {rc.ta:.1f} C, amplicons {rc.amp_min}-{rc.amp_max} bp")
        except ValueError as exc:
            pc.add("ERROR", "settings", str(exc))
        if pc.sites:
            lo = rc.primer_min if rc else primer_length_range(gc)[0]
            hi = rc.primer_max if rc else lo + 14
            comfort = span_limits(cfg.amp_opt, cfg.amp_max, primer_len=min(hi, lo + 4), margin=cfg.margin, flex=20)[0]
            _, hard, gap = span_limits(cfg.amp_opt, cfg.amp_max, primer_len=lo, margin=cfg.margin, flex=20)
            try:
                pc.targets = group_sites(pc.sites, comfort, max(hard, comfort), gap)
                pc.add("OK", "targets", f"{len(pc.sites)} position(s) group into {len(pc.targets)} amplicon "
                                          f"target(s)")
            except ValueError as exc:
                pc.add("ERROR", "targets", f"positions cannot be grouped into amplicons: {exc}")
    elif rows:
        _regions_checks(pc, rows, product_max)

    # variants ----------------------------------------------------------------
    if variants:
        from .variants import load_variants
        try:
            mask, notes = load_variants(variants, meta, variant_min_af)
            pc.variants, pc.variant_notes = mask, notes
            for n in notes:
                pc.add("WARN" if "ignored" in n else "OK", "variants", n)
            if gmod.is_vcf(variants):
                _contig_check(pc, "variants", variants, meta)
                recs, _ = gmod.read_vcf(variants, [])   # bad lines are already reported as ignored
                look = gmod._seq_lookup(meta)
                vr = [{"seqid": look[r["chrom"].lower()]["seqid"], "start": r["pos"], "ref": r["ref"],
                       "line": r["line"]} for r in recs if look.get(r["chrom"].lower()) is not None]
                _ref_check(pc, "variants", vr, meta)
            at = [s for s in pc.sites if mask.within(s.seqid, s.start, s.end)]
            if at:
                pc.add("INFO", "variants", f"{len(at)} requested position(s) carry a known variant (expected for "
                                           f"resistance targets; only primers are kept off variants)")
        except (ValueError, FileNotFoundError, OSError, EOFError) as exc:
            pc.add("ERROR", "variants", str(exc))

    if out:
        _out_check(pc, out, must_exist=regions_mode)
    return pc


def _regions_checks(pc: Precheck, rows, product_max: Optional[int]) -> None:
    """getMultiPrimerSet: one product per region; overlapping regions stop it."""
    rs = sorted(rows, key=lambda r: (r["seqid"], r["start"]))
    over = [f"{a['name']} ({a['start']}-{a['end']}) overlaps {b['name']} ({b['start']}-{b['end']})"
            for a, b in zip(rs, rs[1:]) if a["seqid"] == b["seqid"] and b["start"] < a["end"]]
    pc.add_many("ERROR", "regions", over, "overlaps")
    if over:
        pc.add("INFO", "regions", "merge overlapping regions, or use primerjinn-multiplex, which groups positions "
                                  "into amplicons itself")
    if product_max:
        big = [f"{r['name']} spans {r['end'] - r['start'] + 1} bp, more than --product_size_max {product_max}"
               for r in rows if r["end"] - r["start"] + 1 > product_max]
        pc.add_many("WARN", "regions", big, "long regions")
        close = sum(1 for a, b in zip(rs, rs[1:]) if a["seqid"] == b["seqid"] and 0 < b["start"] - a["end"] < product_max)
        if close:
            pc.add("INFO", "regions", f"{close} pair(s) of regions are closer than {product_max} bp, so their "
                                      f"products may overlap; run overlapping products in separate tubes")
    if not over:
        pc.add("OK", "regions", f"{len(rows)} region(s), none overlapping")
