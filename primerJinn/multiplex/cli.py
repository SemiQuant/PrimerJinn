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
    io.add_argument("--genome", required=True,
                    help="FASTA path, assembly accession (GCF_/GCA_), nucleotide accession or organism name")
    io.add_argument("--positions", required=True,
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
    io.add_argument("--check", action="store_true",
                    help="check the inputs (genome, positions, variants, settings, software, output folder), "
                         "report every problem found, and stop without designing")
    io.add_argument("--fasta_without_primers", action="store_true",
                    help="write amplicon FASTA without the primer-binding sequence")
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
                    help="kapa = KAPA2G Robust HotStart ReadyMix (default); taq; q5; hifi")
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
                  one_pool=a.one_pool, one_pool_max=a.one_pool_max,
                  candidates_per_target=a.candidates, threads=a.threads, seed=a.seed)


def sites_from_positions(source, genome_meta):
    from .precheck import positions_rows
    from .targets import Site

    rows = positions_rows(source, genome_meta)
    return [Site(r["seqid"], r["start"], r["end"], r["name"] or "") for r in rows]


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    from . import report
    from .design import design_pool
    from .precheck import run as precheck

    cfg = config_from_args(a)
    pc = precheck(a.genome, a.positions, cfg, variants=a.variants, variant_min_af=a.variant_min_af,
                  genome_dir=a.genome_dir, out=a.out)
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
