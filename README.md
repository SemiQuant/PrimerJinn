# PrimerJinn

<div>
    <img src="https://user-images.githubusercontent.com/8179171/236663567-94d1f5dc-2ac6-49de-9fc1-a99c7a13945d.png" width="20%" height="20%">
    <p>primerJinn has two main functions: it designs primers for multiplex PCR of given target regions in a DNA sequence (FASTA file) and performs in silico PCR given a list of primers and a reference FASTA file.</p>
</div>

[![DOI](https://doi.org/10.1186/s12859-023-05609-1)](https://doi.org/10.1186/s12859-023-05609-1)

## Features

- **Multiplex Primer Design**: Design optimal primer sets for multiple target regions
- **In Silico PCR**: Simulate PCR reactions with your primers against a reference sequence
- **Primer Interaction Analysis**: Check for primer dimers and cross-reactions
- **Q5 Polymerase Support**: Optimized settings for NEB Q5 hotstart polymerase
- **Illumina Compatibility**: Option to add Illumina partial adapters
- **Comprehensive Output**: Excel workbooks with detailed primer information and interactions

## Installation

```bash
pip install primerJinn
```

### System Requirements

- Python 3.6+
- [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK569861/)

## Usage

### 1. Multiplex Primer Design

```bash
getMultiPrimerSet \
    --region_file "./example/primer_regions.tsv" \
    --input_file "./example/ref.fasta" \
    --target_tm 65 \
    --primer_len 20 \
    --product_size_min 400 \
    --product_size_max 800 \
    --ret 100 \
    --Q5 \
    --output "example" \
    --output_fasta \
    --exclude_primers
```

### 2. In Silico PCR

```bash
PCRinSilico \
   --primer_seq ./example/primers.txt \
   --target_tm 50 \
   --input_file ./example/ref.fasta \
   --output "in_silico_PCR" \
   --output_fasta \
   --exclude_primers
```

### 3. Multiplex designer from genome positions (new, `primerjinn-multiplex`)

Designs a multiplex amplicon-sequencing pool from a genome name (or accession, or FASTA) and a list of positions. It uses modern thermodynamics (SantaLucia 2004 + Owczarzy 2008 Mg2+ correction), template and amplicon secondary structure (ViennaRNA, DNA parameters), genome-wide specificity, primer-dimer optimisation and automatic pooling. Default polymerase: KAPA2G Robust HotStart ReadyMix.

```bash
python3 -m primerJinn.multiplex.cli \
    --genome "Mycobacterium tuberculosis H37Rv" \
    --positions example/mtb_positions.tsv \
    --out results/mtb
```

Positions are 1-based: `7570`, `NC_000962.3:7570`, `761095-761277`, or a TSV/CSV/XLSX/BED file. Key options: `--polymerase {kapa,taq,q5,hifi}`, `--primer_nM`, `--tails {none,illumina,nextera,custom}`, `--amp_min/--amp_opt/--amp_max` (default 225/300/375), `--tm` (default from genome GC). Outputs: `<out>.xlsx`, TSV tables, `<out>_amplicons.fasta`, `<out>.bed` and `<out>_report.md` with the recommended PCR conditions. Design decisions and evidence: `docs/DESIGN.md`, `docs/research/`.

## Parameters

### Common Parameters (Both Tools)

| Parameter | Required | Description | Default |
|-----------|----------|-------------|----------|
| `--input_file` | Yes | Reference FASTA file path | NA |
| `--target_tm` | No | Target melting temperature (°C) | 60 |
| `--salt_concentration` | No | Salt concentration (nM, ignored if Q5=True) | 50 |
| `--output` | No | Output file prefix | 'MultiPlexPrimerSet' or 'in_silico_PCR' |
| `--Q5` | No | Use Q5 polymerase settings | False |
| `--output_fasta` | No | Output amplicon sequences in FASTA format | False |
| `--exclude_primers` | No | Exclude primer sequences from FASTA output | False |

### Primer Design Specific Parameters

| Parameter | Required | Description | Default |
|-----------|----------|-------------|----------|
| `--region_file` | Yes | TSV/XLSX file with regions (name, start, end) | NA |
| `--primer_len` | No | Primer length | 20 |
| `--product_size_min` | No | Minimum amplicon size | 400 |
| `--product_size_max` | No | Maximum amplicon size | 800 |
| `--ret` | No | Maximum primer pairs to return | 100 |
| `--background` | No | Mispriming library FASTA | None |
| `--ill_adapt` | No | Add Illumina partial adapters | False |
| `--clamp` | No | Require GC clamp | 0 |
| `--poly` | No | Max mononucleotide repeat length | 3 |
| `--no_self_background` | No | Skip self-mispriming check | False |

### In Silico PCR Specific Parameters

| Parameter | Required | Description | Default |
|-----------|----------|-------------|----------|
| `--primer_seq` | Yes | File with primer sequences | NA |
| `--product_size_max` | No | Maximum amplicon size | 2000 |
| `--req_five` | No | Require 5' end binding | True |
| `--output_fasta` | No | Output amplicons in FASTA | False |
| `--exclude_primers` | No | Exclude primers from FASTA output | False |

## Output Files

### Multiplex Primer Design
- `MultiPlexPrimerSet.xlsx`: Contains primer pairs, Tm values, and product sizes
- `MultiPlexPrimerSet_amplicons.fasta`: (Optional) FASTA file of amplicon sequences

### In Silico PCR
- `in_silico_PCR.xlsx`: Excel workbook with multiple sheets:
  - `amplicons`: Predicted PCR products
  - `dimers`: Potential primer-dimer interactions
  - `interactions`: Cross-reactions between amplicons
- `in_silico_PCR_amplicons.fasta`: (Optional) FASTA file of predicted amplicons

The FASTA output files contain:
- Headers in format: `>target_name_start_end_length`
- Sequences either:
  - Include primers (default)
  - Exclude primers (with `--exclude_primers`)

## Testing Recommendations

When validating primers for diagnostic assays:
1. Perform qPCR with individual primer pairs
2. Use EVA Green plus (20X in water) and ROX (50X)
3. Analyze melt curves
4. Confirm single bands via agarose gel electrophoresis

## Online Version

A web version of primerJinn's getMultiPrimerSet is available at [DrDx.Me](https://drdx.ucsf.edu/)

## Citation

Limberis, J.D., Metcalfe, J.Z. primerJinn: a tool for rationally designing multiplex PCR primer sets for amplicon sequencing and performing in silico PCR. BMC Bioinformatics 24, 468 (2023). [https://doi.org/10.1186/s12859-023-05609-1](https://doi.org/10.1186/s12859-023-05609-1)

## Author

Jason Limberis (JasonLimberis@ucsf.edu)

## License

This project is licensed under the MIT License.