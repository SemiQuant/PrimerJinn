"""primerjinn-multiplex: design a multiplex amplicon-sequencing primer pool.

Example:
    primerjinn-multiplex --genome "Mycobacterium tuberculosis H37Rv" \\
        --positions example/mtb_positions.tsv --out results/mtb
"""
from __future__ import annotations

import argparse
import sys

from .config import POLYMERASES, TAILS, Config


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="primerjinn-multiplex", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    io = ap.add_argument_group("input and output")
    io.add_argument("--genome",
                    help="FASTA path, assembly accession (GCF_/GCA_), nucleotide accession or organism name")
    io.add_argument("--positions",
                    help="file (txt/tsv/csv/xlsx, BED, or VCF) or a list such as '7570,7582,761155'. "
                         "1-based except BED; seqid:pos and start-end accepted. A VCF gives one target per "
                         "record, covering its REF allele")
    io.add_argument("--out", default="multiplex", help="output prefix (default: multiplex)")
    io.add_argument("--genome_dir", help="genome cache (default ~/.cache/primerjinn/genomes)")
    io.add_argument("--variants",
                    help="known strain variants to keep out of primers: VCF (.vcf/.vcf.gz), BED, or a "
                         "TSV/CSV/TXT of positions like --positions")
    io.add_argument("--variant_min_af", type=float, default=0.0,
                    help="ignore VCF variants with INFO AF below this (default: use all)")
    io.add_argument("--exclude", metavar="IDS",
                    help="positions to leave out: names, loci (seqid:pos, start-end) or target ids from --check "
                         "(T03, J01); comma list or a file with one per line")
    io.add_argument("--vcf_snps_only", action="store_true",
                    help="from a VCF use only single-base variants (REF and ALT one base)")
    io.add_argument("--vcf_max_indel", type=int, default=0,
                    help="leave out VCF records whose REF and ALT differ in length by more than this many bp "
                         "(structural variants count their span; 0 = keep all)")
    io.add_argument("--long_variant", type=int, default=100,
                    help="a VCF deletion longer than this (bp) gets a junction amplicon spanning the deletion "
                         "instead of tiling the deleted region (default 100; 0 = always tile)")
    io.add_argument("--check", action="store_true",
                    help="check the inputs (genome, positions, variants, settings, software, output folder), "
                         "report every problem found, and stop without designing")
    io.add_argument("--fasta_without_primers", action="store_true",
                    help="write amplicon FASTA without the primer-binding sequence")
    fb = ap.add_argument_group("after sequencing: rebalance and redesign (needs --previous)")
    fb.add_argument("--sequencing", metavar="counts.tsv",
                    help="read counts per amplicon from the sequencing run: a TSV/CSV with columns amplicon and "
                         "reads (optional sample), plus a row named unaligned with the reads that match no "
                         "amplicon (see the README)")
    fb.add_argument("--previous", metavar="PREV_primers.tsv",
                    help="the previous design's _primers.tsv (or .xlsx) that was sequenced")
    fb.add_argument("--conc", metavar="conc.tsv", help="primer<TAB>conc_nM actually used, per primer (or give a conc_nM "
                                                       "column in the counts table; default: the previous design's "
                                                       "setting for every primer)")
    fb.add_argument("--replace", metavar="IDS", help="target ids to redesign whatever the data say (comma list "
                                                     "or a file with one per line)")
    fb.add_argument("--keep", metavar="IDS", help="target ids to keep whatever the data say")
    fb.add_argument("--no_redesign", action="store_true",
                    help="only rebalance: write the suggested concentrations and stop")
    fb.add_argument("--redesign_below", type=float, default=0.2,
                    help="an amplicon below this fold of the median is weak (default 0.2): it is redesigned when "
                         "its primers would need more than --conc_max to reach the median, raised otherwise")
    fb.add_argument("--conc_max", type=float, default=500.0,
                    help="highest primer concentration (nM) to suggest; above it the pair is redesigned")
    fb.add_argument("--alpha", type=float, default=0.5, help="damping of the concentration change (default 0.5)")
    fb.add_argument("--cap", type=float, default=4.0, help="largest fold change of one concentration (default 4)")
    fb.add_argument("--reaction_ul", type=float, default=25.0, help="reaction volume for the stock volumes")
    fb.add_argument("--stock_uM", type=float, default=10.0, help="primer stock concentration for the volumes")
    am = ap.add_argument_group("amplicons and primers")
    am.add_argument("--amp_min", type=int, default=225)
    am.add_argument("--amp_opt", type=int, default=300)
    am.add_argument("--amp_max", type=int, default=375)
    am.add_argument("--margin", type=int, default=3, help="clear bases between a primer and a required position")
    am.add_argument("--tm", type=float, help="target primer Tm on the design scale (default from genome GC)")
    am.add_argument("--tm_scale", choices=["physical", "neb"],
                    help="physical (default for Taq-type enzymes) or neb (NEB Q5 calculator, default for q5)")
    am.add_argument("--tm_window", type=float, default=1.5)
    am.add_argument("--primer_min", type=int)
    am.add_argument("--primer_max", type=int)
    am.add_argument("--max_poly", type=int, help="longest mononucleotide run in a primer (default 4, 5 below 25%% GC)")
    am.add_argument("--primer_gc_min", type=float,
                    help="lowest primer GC fraction (default min(0.20, genome GC - 0.08), at least 0.10)")
    am.add_argument("--tails", choices=sorted(TAILS) + ["custom"], default="none")
    am.add_argument("--tail_fwd", default="", help="5' tail for forward primers (with --tails custom)")
    am.add_argument("--tail_rev", default="", help="5' tail for reverse primers (with --tails custom)")
    rx = ap.add_argument_group("reaction")
    rx.add_argument("--polymerase", choices=sorted(POLYMERASES), default="kapa",
                    help="kapa = generic Taq master mix (default); taq = standard Taq, 1.5 mM Mg2+; q5; hifi")
    rx.add_argument("--mono_mM", type=float, default=50.0, help="monovalent cations (K+ + Na+ + Tris+/2)")
    rx.add_argument("--mg_mM", type=float, help="total Mg2+ (default from polymerase preset)")
    rx.add_argument("--dntp_mM", type=float, default=0.8, help="total dNTP (all four)")
    rx.add_argument("--primer_nM", type=float, default=100.0, help="each primer in the pool")
    rx.add_argument("--dmso", type=float, default=0.0, help="DMSO %% v/v")
    rx.add_argument("--ta", type=float, help="design annealing temperature (default derived)")
    am.add_argument("--no_balance", action="store_true",
                    help="do not extend weakly binding primers at the 5' end (lowercase bases)")
    am.add_argument("--max_ext", type=int, default=8, help="most bases added to one primer by binding balance")
    sc = ap.add_argument_group("search")
    sc.add_argument("--no_specificity", action="store_true", help="skip the genome-wide off-target check")
    sc.add_argument("--outer_min", type=int, default=1000,
                    help="targets closer than this (outer product length) go to different pools")
    sc.add_argument("--one_pool", action="store_true",
                    help="force a single pool: targets that may not share a tube are merged into one amplicon; only "
                         "those get a larger size class (same spread), the rest keep --amp_min..--amp_max")
    sc.add_argument("--one_pool_max", type=int, default=1000,
                    help="largest amplicon allowed with --one_pool (default 1000)")
    sc.add_argument("--one_pool_drop", action="store_true",
                    help="with --one_pool: targets that still need a second pool are dropped (listed in the "
                         "report) instead of going to pool 2")
    sc.add_argument("--candidates", type=int, default=40, help="candidate pairs kept per target")
    sc.add_argument("--threads", type=int, default=0)
    sc.add_argument("--seed", type=int, default=1)
    sc.add_argument("--quiet", action="store_true")
    return ap


def config_from_args(a) -> Config:
    return Config(amp_min=a.amp_min, amp_opt=a.amp_opt, amp_max=a.amp_max, margin=a.margin, tm=a.tm,
                  tm_scale=a.tm_scale, balance=not a.no_balance, max_ext=a.max_ext, max_poly=a.max_poly, primer_gc_min=a.primer_gc_min,
                  tm_window=a.tm_window, primer_min=a.primer_min, primer_max=a.primer_max, tails=a.tails,
                  tail_fwd=a.tail_fwd.upper(), tail_rev=a.tail_rev.upper(), polymerase=a.polymerase,
                  mono_mM=a.mono_mM, mg_mM=a.mg_mM, dntp_mM=a.dntp_mM, primer_nM=a.primer_nM, dmso_pct=a.dmso,
                  ta=a.ta, specificity=not a.no_specificity, outer_min=a.outer_min,
                  one_pool=a.one_pool, one_pool_max=a.one_pool_max, one_pool_drop=a.one_pool_drop,
                  vcf_snps_only=a.vcf_snps_only, long_variant=a.long_variant, vcf_max_indel=a.vcf_max_indel,
                  candidates_per_target=a.candidates, threads=a.threads, seed=a.seed)


def _id_list(spec):
    """Comma-separated ids, or a file with one id per line."""
    if not spec:
        return []
    import os
    if os.path.exists(spec):
        with open(spec) as fh:
            return [x.strip() for x in fh if x.strip() and not x.startswith("#")]
    return [x.strip() for x in spec.split(",") if x.strip()]


def sites_from_positions(source, genome_meta):
    from .precheck import positions_rows
    from .targets import Site

    rows = positions_rows(source, genome_meta)
    return [Site(r["seqid"], r["start"], r["end"], r["name"] or "", gap=r.get("gap")) for r in rows]


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    from . import report
    from .design import design_pool
    from .precheck import run as precheck

    cfg = config_from_args(a)

    # ---- sequencing feedback: rebalance, and decide what to redesign
    replace, keep = _id_list(a.replace), _id_list(a.keep)
    if (a.sequencing or replace or keep) and not a.previous:
        print("ERROR: --sequencing, --replace and --keep need --previous (the sequenced design's _primers.tsv)",
              file=sys.stderr)
        return 2
    fb = None
    if a.sequencing:
        from . import feedback

        try:
            run = feedback.load_run(a.sequencing)
            prev = feedback.load_previous(a.previous, a.conc, default_nM=a.primer_nM)
            th = feedback.Thresholds(redesign_below=a.redesign_below, conc_max=a.conc_max, alpha=a.alpha, cap=a.cap,
                                     reaction_ul=a.reaction_ul, stock_uM=a.stock_uM)
            fb = feedback.analyse(run, prev, th, replace=replace, keep=keep)
        except (ValueError, FileNotFoundError, KeyError) as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        files = feedback.write(fb, a.out, default_nM=cfg.primer_nM)
        if not a.quiet:
            for k, v in fb.summary:
                print(f"{k}: {v}")
        replace = fb.replace
        if a.no_redesign or not replace:
            print("")
            print("no target to redesign; wrote the rebalance tables:" if not replace and not a.no_redesign
                  else "wrote the rebalance tables:")
            for f in files:
                print(f"  {f}")
            if replace:
                print(f"{len(replace)} target(s) would be redesigned without --no_redesign: {', '.join(replace)}")
            return 0
    if a.previous and not (a.genome and a.positions):
        print("ERROR: a redesign needs --genome and --positions (the same as the previous design)", file=sys.stderr)
        return 2
    if a.previous:
        from . import redesign

        try:
            cfg.resolve(0.5)   # tails only; GC-dependent values are redone by the designer
            if keep and not replace:   # --keep alone: everything else is redesigned
                ids = redesign.read_primers_table(a.previous)["target"].astype(str)
                replace = [t for t in dict.fromkeys(ids) if t not in set(keep)]
            prev_design = redesign.read_previous(a.previous, replace, keep or None, cfg.tail_fwd, cfg.tail_rev)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        cfg = config_from_args(a)
        redesign.apply(cfg, prev_design)
        if not replace:
            print("nothing to redesign (--replace is empty); keeping every primer", file=sys.stderr)
    if not (a.genome and a.positions):
        print("ERROR: --genome and --positions are required", file=sys.stderr)
        return 2
    pc = precheck(a.genome, a.positions, cfg, variants=a.variants, variant_min_af=a.variant_min_af,
                  genome_dir=a.genome_dir, out=a.out, exclude=_id_list(a.exclude))
    if a.check or not a.quiet or not pc.ok or pc.warnings:
        print(pc.text(only_problems=not a.check and a.quiet), file=sys.stdout if pc.ok else sys.stderr)
    if not pc.ok:
        return 2
    if a.check:
        return 0
    meta, sites = pc.meta, pc.sites
    cfg.variants = pc.variants
    try:
        d = design_pool(sites, meta, cfg, verbose=not a.quiet)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    files = report.write(d, a.out, fasta_with_primers=not a.fasta_without_primers)
    if fb is not None:
        from . import feedback

        new = [nm for tid, sides in d.names.items() if tid not in cfg.fixed for nm in sides.values()]
        files += feedback.write(fb, a.out, new_primers=new, default_nM=cfg.primer_nM)
    print("")
    print(f"{len(d.chosen)} amplicons in {max(d.pools.values())} pool(s); recommended Ta {d.ta['recommended']:g} C")
    if d.coverage_problems:
        print(f"WARNING: {len(d.coverage_problems)} position(s) NOT covered, see the report")
    for w in d.warnings:
        print(f"warning: {w}")
    print("wrote:")
    for f in files:
        print(f"  {f}")
    if d.dropped_sites:
        print(f"WARNING: {len(d.dropped_sites)} position(s) are N in the reference and were not designed")
    return 1 if (d.coverage_problems or d.dropped_sites) else 0


if __name__ == "__main__":
    sys.exit(main())
