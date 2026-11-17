# Changelog

## 2.1.0 (2026-11-17)

New: `primerjinn-multiplex`, a multiplex amplicon-sequencing pool designer, and `primerjinn-gui`, its local web interface.

- Design from a genome (NCBI name, accession or FASTA) and a list of positions (table, BED or VCF); every position is covered, amplicons 225 to 375 bp by default.
- Primers matched on Tm and on binding at the annealing temperature, with template secondary structure (ViennaRNA), dimers, unintended products and genome-wide specificity.
- Pools by graph colouring, or `--one_pool` with size classes; junction amplicons for long deletions; VCF filters for SNPs only and for indel length.
- Recommended annealing temperature and PCR advice for a generic Taq master mix or Q5.
- After sequencing: `--sequencing counts.tsv --previous` rebalances primer concentrations and redesigns only the pairs that cannot be fixed by concentration, keeping the working primers fixed.
- Outputs: XLSX and TSV tables, order sheet, FASTA, BED, uncovered-position table, amplicon map (PNG) and a Markdown report.
- Packaging: `pyproject.toml`, version in `primerJinn.__version__`, no BLAST+ download at install time (BLAST+ is needed only by `getMultiPrimerSet` and `PCRinSilico`; install it with conda or from NCBI).

## 2.0.3

Last release of the original tools (`getMultiPrimerSet`, `PCRinSilico`, `primerJinn`).
