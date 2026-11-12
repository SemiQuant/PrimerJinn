# Resolving a genome name or accession to a FASTA (PrimerJinn multiplex rebuild)

All requests below were run live on 2026-11-10 and re-run on 2026-11-11 from macOS arm64 with curl and Python 3.13.2 urllib (pysam 0.23.0). No API key, no email, no contact detail was sent; the only identifying header is `User-Agent: primerjinn-genome-resolver/1.0` and E-utilities calls carry `tool=primerjinn`. Files:

- Module (version 2, the one to use): `/private/tmp/claude-501/-Users-semiquant--claude-spares-spare-n9Plyg/5e124d61-7974-4f46-b4a7-ff4fbdb28e2f/scratchpad/research/genome_reference.py`
- Version 1, kept for comparison: `/private/tmp/claude-501/-Users-semiquant--claude-spares-spare-n9Plyg/5e124d61-7974-4f46-b4a7-ff4fbdb28e2f/scratchpad/research/genome_reference_v1.py`
- Test script (2026-11-11 run, fresh cache): `/private/tmp/claude-501/-Users-semiquant--claude-spares-spare-n9Plyg/5e124d61-7974-4f46-b4a7-ff4fbdb28e2f/scratchpad/research/gf_v3/test_v3.py`
- Full test output: `/private/tmp/claude-501/-Users-semiquant--claude-spares-spare-n9Plyg/5e124d61-7974-4f46-b4a7-ff4fbdb28e2f/scratchpad/research/gf_v3/test_output.txt` (exit 0, 44 s wall time with all downloads from a cold cache)
- Adversarial tests of version 2: `/private/tmp/claude-501/-Users-semiquant--claude-spares-spare-n9Plyg/5e124d61-7974-4f46-b4a7-ff4fbdb28e2f/scratchpad/research/gf_v2/adv_f.py` and `adv_f.out`

Version 2 fixes found by adversarial testing of version 1 (listed in the module docstring as F1 to F13): influenza names containing "/" were taken for file paths; thousands separators and a numeric name column could silently produce wrong or megabase-wide targets (now errors, plus `max_span` 50,000 bp); `pos` + `end` header columns; empty, non-FASTA and duplicate-id files; per-query notes on cache hits; bounds check in `fetch()`; BED track lines; a `ref` column wrongly taken as seqid; a per-genome `fcntl` lock file for concurrent downloads; scaffold counts by NCBI role; a missing position file now raises `FileNotFoundError`.

## 1. Decisions in brief

1. Resolution order: existing local file, then assembly accession (`GCF_`/`GCA_`), then nucleotide accession, then organism name or taxid through NCBI Taxonomy to the RefSeq reference genome.
2. Names go through NCBI Datasets v2. E-utilities is used only for single nucleotide records (esummary + efetch) and as a last-resort name search (esearch db=assembly).
3. One reference genome: use it. Several: stop and list them. None: use the only RefSeq assembly if there is exactly one, otherwise stop and list.
4. Cache under `~/.cache/primerjinn/genomes/<accession>/` (override with `--genome_dir` or `PRIMERJINN_CACHE`), plain upper-case FASTA plus `.fai` plus `metadata.json`.
5. Positions are 1-based inclusive. Bare integers are accepted only when the genome has one sequence or exactly one chromosome (plus plasmids). Otherwise `seqid:pos` is required.
6. Genome GC is computed from the downloaded FASTA, not taken from NCBI, because NCBI rounds `gc_percent` to 0.5.

## 2. NCBI Datasets v2 REST API

Base: `https://api.ncbi.nlm.nih.gov/datasets/v2`. Docs: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/

### 2.1 Taxonomy name lookup

```
GET /taxonomy/taxon/{name-or-taxid}        (URL-encode the name)
```

Fields used: `taxonomy_nodes[0].taxonomy.{tax_id, organism_name, rank, blast_name}`. On an unknown name the HTTP status is still **200** and the node holds `errors[0].reason` instead of `taxonomy`.

| Query | Result |
|---|---|
| `Mycobacterium tuberculosis H37Rv` | tax_id 83332, rank STRAIN |
| `Mycobacterium tuberculosis` | tax_id 1773 |
| `SARS-CoV-2` | tax_id 2697049 (matched as acronym) |
| `Homo sapiens`, `human`, `9606` | tax_id 9606 |
| `Plasmodium falciparum 3D7` | tax_id 36329, rank ISOLATE |
| `Staphylococcus aureus NCTC 8325` | tax_id 93061, organism_name "Staphylococcus aureus subsp. aureus NCTC 8325" |
| `E. coli` | tax_id 562 (common name) |
| `Escherichia coli K-12 MG1655` | **error**: "not a recognized NCBI Taxonomy name" |
| `Escherichia coli str. K-12 substr. MG1655` | tax_id 511145 |
| `tuberculosis`, `notanorganism xyz` | error, HTTP 200 |

So the exact lookup accepts scientific names, common names, acronyms and taxids, but informal strain strings fail. Fallback:

```
GET /taxonomy/taxon_suggest/{text}
-> sci_name_and_ids[]: {sci_name, tax_id, matched_term, rank, group_name}
```

`taxon_suggest/Escherichia coli K-12 MG1655` returns "Escherichia coli str. K-12 substr. MG1655" (511145) first, followed by unrelated phages. The suggestions are fuzzy, so the module accepts a suggestion only when every alphanumeric token of the query occurs in `sci_name` (or `matched_term`) and exactly one suggestion passes. Otherwise it lists the suggestions and stops.

### 2.2 Genome report by taxon

```
GET /genome/taxon/{taxid}/dataset_report?filters.reference_only=true&page_size=20
```

Fields used per `reports[]` item: `accession`, `paired_accession`, `organism.{organism_name, tax_id, infraspecific_names}`, `assembly_info.{assembly_name, assembly_level, refseq_category, assembly_status}`, `assembly_stats.{total_sequence_length, gc_percent, gc_count, atgc_count, total_number_of_chromosomes, number_of_component_sequences}`. An empty result is `{}` with HTTP 200 (no `reports` key).

| taxid | `reference_only=true` result |
|---|---|
| 83332 (M. tuberculosis H37Rv) | GCF_000195955.2, ASM19595v2, Complete Genome, "reference genome", 4,411,532 bp, gc_percent 65.5 |
| 1773 (M. tuberculosis) | the same single record, GCF_000195955.2 |
| 511145 (E. coli MG1655) | GCF_000005845.2, 4,641,652 bp, gc 51 |
| 562 (E. coli) | **2 records**: GCF_000005845.2 (MG1655) and GCF_000008865.2 (O157:H7 Sakai, 3 sequences) |
| 93061 (S. aureus NCTC 8325) | GCF_000013425.1, 2,821,361 bp, gc 33 |
| 36329 (P. falciparum 3D7) | GCF_000002765.6, 23,292,622 bp, 14 chromosomes, gc 19.5 |
| 9606 (human) | GCF_000001405.40, GRCh38.p14, 3,099,441,038 bp, gc 41 |
| 2697049 (SARS-CoV-2) | **empty** |
| 83331 (M. tuberculosis CDC1551) | **empty** |

Answer to the H37Rv question: "Mycobacterium tuberculosis H37Rv" and "Mycobacterium tuberculosis" both resolve to GCF_000195955.2 (sequence NC_000962.3). The species name works because the taxon query includes child taxa and only H37Rv carries the reference flag.

When `reference_only` is empty the module retries with:

```
filters.assembly_source=refseq&filters.assembly_level=complete_genome&filters.exclude_atypical=true&filters.assembly_version=current
```

- SARS-CoV-2: 1 record, GCF_009858895.2 (contains NC_045512.2, 29,903 bp). Accepted with a note.
- M. tuberculosis CDC1551: 1 record, GCF_000008585.1. Accepted with a note.
- Influenza A virus (11320): 7 records. The module stops and lists them.

The separate virus endpoint `GET /virus/taxon/{taxid}/dataset_report?filter.refseq_only=true` also works (returns NC_045512.2 for SARS-CoV-2, 56 segment records for influenza A) but returns nucleotide records, not assemblies, so segmented viruses are harder to group. The genome endpoint is enough.

### 2.3 Genome report and sequences by accession

```
GET /genome/accession/{GCF_or_GCA}/dataset_report
GET /genome/accession/{acc}/sequence_reports?page_size=1000   (next_page_token for more)
```

- `GCF_000195955` (no version) resolves to the current GCF_000195955.2.
- `GCA_000195955.2` works and reports `paired_accession` GCF_000195955.2.
- Unknown (`GCF_999999999.9`) and superseded (`GCF_000195955.1`) accessions return `{}` with HTTP 200.

`sequence_reports` fields used: `refseq_accession`, `genbank_accession`, `chr_name`, `assigned_molecule_location_type` (Chromosome, Plasmid, Mitochondrion ...), `role` (assembled-molecule, unlocalized-scaffold, ...), `length`, `gc_percent`. Examples:

```
GCF_000008865.2 (E. coli O157:H7 Sakai)
NC_002695.2 BA000007.3 chromosome Chromosome assembled-molecule 5498578 50.5
NC_002127.1 AB011548.2 pOSAK1     Plasmid    assembled-molecule    3306 43.5
NC_002128.1 AB011549.2 pO157      Plasmid    assembled-molecule   92721 47.5
```

P. falciparum 3D7 has 14 chromosome records (NC_004325.2 chr 1, 640,851 bp ... NC_037283.1 chr 14, 3,291,936 bp). Human GRCh38.p14 returns 705 sequence records: 24 chromosomes, 1 mitochondrion, 680 other records by `role` (260 alt-scaffold, 164 fix-patch, 90 novel-patch, 40 unlocalized-scaffold, 126 unplaced-scaffold).

### 2.4 Download

```
GET /genome/accession/{acc}/download_summary?include_annotation_type=GENOME_FASTA
GET /genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA
```

`download` returns `content-type: application/zip`, `content-disposition: attachment; filename=ncbi_dataset.zip`. For GCF_000195955.2: 1,281,230 bytes in 0.69 s, containing:

```
README.md
ncbi_dataset/data/assembly_data_report.jsonl
ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna   (4,466,740 bytes)
ncbi_dataset/data/dataset_catalog.json
md5sum.txt
```

The module extracts the single `*_genomic.fna`, checks it against `md5sum.txt`, and upper-cases sequence lines.

`download_summary` for human GCF_000001405.40 reports `hydrated.estimated_file_size_mb` 927 and `available_files.all_genomic_fasta.size_mb` 927.8. The module checks this before downloading (`max_download_mb`, default 1500). The human genome was not downloaded.

### 2.5 Rate limit and failure modes (Datasets)

- Response headers: `x-ratelimit-limit: 5`, `x-ratelimit-remaining: n`. The docs state 5 requests per second without a key and 10 with one (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys/).
- 12 simultaneous keyless requests gave 7 x HTTP 429 and 5 x HTTP 200 (same result on both days).
- Doc quote checked 2026-11-11: "By default, this rate limit is set at 5 requests per second (rps)"; a key raises it to 10 rps.
- Unknown taxon name: HTTP 200 with `errors`.
- Unknown or superseded assembly: HTTP 200 with `{}`.
- Download of an unknown accession: HTTP 200 with a 971 byte body that has a zip magic number but is not a valid archive (`zipfile.is_zipfile` is False). Always validate.

## 3. NCBI E-utilities

Base: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils`. Docs: https://www.ncbi.nlm.nih.gov/books/NBK25497/ and https://www.ncbi.nlm.nih.gov/books/NBK25499/ (the NBK25497 page returned a reCAPTCHA check to automated fetches, which I did not try to bypass). The limit is documented in the NCBI announcement https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ : any IP that "submits more than three E-utility requests per second will receive an error message", and a key raises this to 10 per second. It is also confirmed by the live header and 429 responses below.

```
GET /esummary.fcgi?db=nuccore&id=NC_000962.3&retmode=json
GET /efetch.fcgi?db=nuccore&id=NC_000962.3&rettype=fasta&retmode=text
GET /esearch.fcgi?db=assembly&term=<name> AND "reference genome"[filter] AND latest[filter]&retmode=json
GET /esummary.fcgi?db=assembly&id=<uid>&retmode=json      -> assemblyaccession
```

esummary nuccore fields used: `accessionversion`, `title`, `slen`, `taxid`, `organism`, `strain`, `moltype`, `topology`, `genome`, `status`, `replacedby`.

| Request | Result |
|---|---|
| esummary NC_000962.3 | slen 4,411,532, taxid 83332, sourcedb refseq, moltype dna |
| esummary MN908947.3 | slen 29,903, taxid 2697049, sourcedb insd, moltype rna |
| efetch NC_000962.3 | 4,474,618 bytes, 1.07 s, header `>NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome` |
| efetch MN908947.3 | 30,428 bytes, 0.34 s |
| efetch NC_000962 (no version) | returns the current NC_000962.3 |
| efetch NC_999999.9 | **HTTP 400**, URL-encoded text "Failed to retrieve sequence" |
| esummary NC_999999.9 | HTTP 200, `"error":"Invalid uid NC_999999.9 at position= 0"` |
| esearch assembly "Escherichia coli K-12 MG1655 AND reference genome[filter]" | 1 hit, uid 79781, assemblyaccession GCF_000005845.2 |

Rate limit: headers `x-ratelimit-limit: 3`. 12 simultaneous keyless requests gave 3 x HTTP 200 and 9 x HTTP 429 with body `{"error":"API rate limit exceeded","count":"4","limit":"3"}` (the body also echoes the caller IP in an `api-key` field). So 3 requests per second without a key is confirmed. Re-run 2026-11-11 with 10 simultaneous `einfo` calls: 2 x 200, 8 x 429.

efetch gives no checksum, so the module compares the downloaded length with esummary `slen` and rejects truncated downloads.

Identity checks run on the downloaded files: MN908947.3 equals NC_045512.2 (29,903 bp, identical strings); NC_000962.3 equals AL123456.3 (identical).

## 4. Resolution logic implemented

```
resolve_genome(name_or_path, cache_dir=None, download=True, max_download_mb=1500, refresh=False)
```

1. **Local file** if `Path(query).is_file()`. If the string looks like a path (contains `/` or ends in .fa/.fna/.fasta/.fas[.gz]) but does not exist: error "File not found", no network lookup.
2. **Assembly accession**: `^GC[AF]_\d{9}(\.\d+)?$`. dataset_report, sequence_reports, zip download.
3. **Nucleotide accession**: `^(?:[A-Z]{2}_[A-Z]{0,6}\d{5,}|[A-Z]\d{5}|[A-Z]{2}\d{6}|[A-Z]{2}\d{8})(\.\d+)?$` (RefSeq NC_/NZ_/NG_/NT_/NW_/AC_, and INSDC 1+5, 2+6, 2+8 formats). esummary then efetch. Only that one sequence is fetched, because the user named one sequence.
4. **Name or taxid**: taxonomy exact lookup, else taxon_suggest with the token rule, else esearch db=assembly. Then reference_only, then RefSeq complete, then any RefSeq assembly.

Errors: `GenomeResolutionError` (not found, too large, truncated) and its subclass `AmbiguousGenomeError` with `.alternatives`. Example message:

```
'Escherichia coli' (taxid 562) has 2 RefSeq reference genomes. Pass one accession:
  GCF_000005845.2  Escherichia coli str. K-12 substr. MG1655  [ASM584v2; Complete Genome; reference genome; 4641652 bp]
  GCF_000008865.2  Escherichia coli O157:H7 str. Sakai  [ASM886v2; Complete Genome; reference genome; 5594605 bp]
```

HTTP layer: stdlib urllib, User-Agent `primerjinn-genome-resolver/1.0`, per-host limiter (4 req/s Datasets, 2.5 req/s E-utilities), retry with exponential backoff (1, 2, 4, 8, 16 s) on 429, 5xx and network errors, honouring `Retry-After` when present.

Returned metadata keys: `query, source, accession, organism, tax_id, assembly_name, assembly_level, refseq_category, fasta, fai, sequences[{seqid, length, gc_percent, type, role, name, aliases}], n_sequences, total_length, gc_percent, gc_percent_ncbi, default_seqid, default_unambiguous, notes, warnings`.

## 5. Test results (code actually run, 2026-11-11, version 2, cold cache)

| Query | Source | Accession | Sequences | Length (bp) | GC % (computed) | Time (s, cold) |
|---|---|---|---|---|---|---|
| Mycobacterium tuberculosis H37Rv | taxonomy | GCF_000195955.2 | 1 (NC_000962.3) | 4,411,532 | 65.61 | 3.5 |
| Mycobacterium tuberculosis | taxonomy | GCF_000195955.2 | 1 | 4,411,532 | 65.61 | 1.4 |
| NC_000962.3 | nuccore | NC_000962.3 | 1 | 4,411,532 | 65.61 | 2.0 |
| GCF_000195955.2 | assembly | GCF_000195955.2 | 1 | 4,411,532 | 65.61 | 0.6 |
| MN908947.3 | nuccore | MN908947.3 | 1 | 29,903 | 37.97 | 0.9 |
| SARS-CoV-2 | taxonomy | GCF_009858895.2 | 1 (NC_045512.2) | 29,903 | 37.97 | 2.6 |
| Escherichia coli K-12 MG1655 | taxonomy (suggest) | GCF_000005845.2 | 1 (NC_000913.3) | 4,641,652 | 50.79 | 3.2 |
| Staphylococcus aureus NCTC 8325 | taxonomy | GCF_000013425.1 | 1 (NC_007795.1) | 2,821,361 | 32.87 | 2.3 |
| Plasmodium falciparum 3D7 | taxonomy | GCF_000002765.6 | 14 | 23,292,622 | 19.34 | 3.9 |
| GCF_000008865.2 | assembly | GCF_000008865.2 | 3 (1 chromosome, 2 plasmids) | 5,594,605 | 50.48 | 2.5 |
| NC_000962 | nuccore | NC_000962.3 | 1 | 4,411,532 | 65.61 | 0.2 |
| GCA_000195955 | assembly | GCA_000195955.2 | 1 (AL123456.3) | 4,411,532 | 65.61 | 4.1 |
| Homo sapiens / human / 9606 (download=False) | taxonomy | GCF_000001405.40 (GRCh38.p14) | 705 (24 chromosomes, 1 mitochondrion, 260 alt-scaffold, 164 fix-patch, 90 novel-patch, 40 unlocalized, 126 unplaced) | 3,099,441,038 | 41 (NCBI) | n/a |

Human FASTA zip is about 927 MB; `max_download_mb=100` correctly refused it. A cache hit for a repeated query takes under 1 ms. The whole test cache (9 genome folders) is 66 MB.

Failure cases, all with a clear message: `Escherichia coli` (2 references listed), `Influenza A virus` (7 RefSeq assemblies listed), `tuberculosis` (suggestions listed), `notanorganism xyz`, `GCF_999999999.9`, `NC_999999.9`, `GCF_000195955.1` (superseded), `/no/such/file.fa`.

Notable: the species-level name returns a note, for example `'Mycobacterium tuberculosis' (taxid 1773) resolved to the reference genome of 'Mycobacterium tuberculosis H37Rv' (taxid 83332)`. The tool should print notes and warnings to the user.

## 6. Multi-sequence genomes

- `seqid:pos` and `seqid:start-end` are the explicit forms. `seqid` matches case-insensitively against the FASTA id, the versionless accession (`NC_000962`), the GenBank accession (`AL123456.3`), and the NCBI molecule name (`pO157`, `7`, and `chr7` for numeric or X/Y/Z/W chromosome names). Tested: `chr7:403625` and `13:1725259` on P. falciparum, `pO157:500` on E. coli Sakai.
- Bare integer rule (implemented as `default_seqid` + `default_unambiguous`):
  - 1 sequence: use it, no warning.
  - Exactly 1 chromosome plus plasmids or organelles: use the chromosome, with a warning naming it and its length.
  - 2 or more chromosomes (P. falciparum 14, human 24): **error** that asks for `seqid:pos` and names the longest sequence as an example. Silently choosing chromosome 14 of P. falciparum would nearly always be wrong. `allow_ambiguous_default=True` restores "longest sequence" behaviour if the implementer prefers the looser rule from the brief.
- For local FASTA files without NCBI metadata, a sequence is treated as a plasmid or organelle when its header description contains "plasmid", "mitochondri", "chloroplast", "apicoplast" or "plastid". This is a heuristic.
- A nucleotide accession and the assembly that contains it give different FASTA ids only when one is GenBank and the other RefSeq (MN908947.3 versus NC_045512.2). Positions are the same because the sequences are identical.

## 7. GC computation

`GC % = 100 x (G + C) / (A + C + G + T)`, per sequence and for the whole genome, N and IUPAC codes excluded, after upper-casing. Computed with `str.count` in 50 Mb chunks through pysam.

| Genome | Computed | NCBI `gc_percent` |
|---|---|---|
| M. tuberculosis H37Rv | 65.61 | 65.5 (gc_count 2,894,612 / atgc_count 4,411,532 = 65.61) |
| E. coli MG1655 | 50.79 | 51 |
| S. aureus NCTC 8325 | 32.87 | 33 |
| P. falciparum 3D7 | 19.34 | 19.5 |
| SARS-CoV-2 | 37.97 | not checked |

NCBI rounds `gc_percent` to the nearest 0.5, so the tool should use its own value. For the Tm target the whole-genome GC is the right default, but per-sequence GC is also stored (E. coli Sakai plasmid pOSAK1 is 43.41 % against 50.54 % for the chromosome), and local GC around each target matters more for individual primers.

## 8. Caching, gzip and indexing

- Location: `$PRIMERJINN_CACHE`, else `$XDG_CACHE_HOME/primerjinn/genomes`, else `~/.cache/primerjinn/genomes`. Expose as `--genome_dir`.
- Layout: `<cache>/<accession>/<accession>.fna`, `.fna.fai`, `metadata.json`; `<cache>/queries.json` maps the normalised query to an accession.
- Versioned accessions (`NC_000962.3`, `GCF_000195955.2`) are cached permanently (`pinned`). Names and versionless accessions expire after 30 days, since NCBI can change the reference genome or release a new version. `refresh=True` forces a new lookup.
- Writes are atomic (`.part<pid>` then `os.replace`).
- Stored FASTA is plain and upper-case. NCBI eukaryote FASTA is soft-masked in lower case; upper-casing avoids case bugs in Tm and structure code. The masking information is lost; if repeat masking is wanted later it needs a separate file.
- Local gzip input: `pysam.FastaFile` cannot index plain gzip (only BGZF). Both plain gzip and BGZF inputs are decompressed once into `<cache>/local/<stem>-<hash>.fna`, where the hash covers path, size and mtime. Tested with both.
- Local plain FASTA: indexed in place with `pysam.faidx`. If the folder is read-only the index goes to `<cache>/local_index/` using `pysam.faidx(path, "--fai-idx", alt)` and is opened with `FastaFile(path, filepath_index=alt)`. Tested on a read-only folder. (`FastaFile(..., filepath_index=...)` alone does not create a missing index; it raises OSError.)
- Retrieval: `fetch(genome, seqid, start, end)` is 1-based inclusive and maps to `fa.fetch(seqid, start-1, end)`. Check: `NC_000962.3:1-20` = `TTGACCGATGACCCCGGTTC`.
- pysam 0.23.0 is already a dependency and is faster than Biopython `SeqIO.index` for random access; Biopython is not needed here.

## 9. Position input formats

All coordinates 1-based inclusive, matching `example/primer_regions.tsv` (header `start end`) and `example/mtb_positions.tsv` (header `name pos`). Implemented in `parse_positions(source, genome)` returning `[{name, seqid, start, end}]`.

| Form | Example | Why |
|---|---|---|
| Bare integer, one per line | `761155` | Simplest case for a single-chromosome bacterium or virus |
| Range | `761007-761277` (also `a..b`) | A target wider than one base, such as a resistance region |
| With sequence | `NC_000962.3:761155`, `chr7:403625-403700` | Needed for plasmids and multi-chromosome genomes; same syntax as samtools regions |
| Name plus locus on a line | `katG_S315<TAB>2155168`, `rpoB,761155` | Names flow into primer and amplicon names |
| Header table (TSV, CSV, semicolon, XLSX) | columns `name`, `seqid`/`chrom`/`chr`/`contig`, `pos`/`position`, `start`, `end` | Matches the existing example files and spreadsheets from the lab |
| BED (`.bed`, `.bed.gz`) | `NC_000962.3 9 20 a` | 0-based half-open, converted to 1-based (start + 1); applied only by file extension or explicit flag |

Rules: lines starting with `#` and blank lines are skipped; a UTF-8 BOM is tolerated; the first row is a header when it contains a known column name and no locus token; a missing `end` means a single position; duplicates (same seqid, start, end) are dropped keeping the first name; a missing name becomes `seqid_start` or `seqid_start_end`. Every row is validated (`1 <= start <= end <= sequence length`, known seqid) and errors name the row:

```
Position row 1 ['1000']: the genome has 14 sequences, so the position needs a sequence id, for example NC_037283.1:1000
Position row 1 ['5000000']: 5000000 is beyond the end of NC_000962.3 (4411532 bp)
Position row 1 ['chrX:5']: unknown sequence 'chrX'. Known: NC_000962.3
Position row 1 ['200-100']: start 200 and end 100 must satisfy 1 <= start <= end
```

Tested: `example/mtb_positions.tsv` gives 24 targets, `example/primer_regions.tsv` gives 8 regions, plus list input, CSV with empty cells, BED and XLSX.

Not supported on purpose: thousands separators (`761,155`), because a comma is also the CSV delimiter.

## 10. Limits and things not verified

- The E-utilities documentation text for the 3 requests per second limit could not be retrieved by curl here; the limit is confirmed from the live header and the 429 body only.
- Names with no exact taxonomy match depend on `taxon_suggest`, whose ranking is not documented as stable. The token rule keeps it safe (it refuses rather than guesses), and the esearch db=assembly fallback was verified only for the E. coli MG1655 string by hand, not through the module path (the suggest path caught it first).
- Segmented viruses resolve only by assembly accession or a strain-level name; a species-level name such as "Influenza A virus" lists the 7 RefSeq assemblies.
- Large genomes: the human path was tested to resolution and size guard only. The 50 Mb chunked GC scan and the 705-sequence index were not timed on a real 3.1 Gb file.
- Concurrency: version 2 takes an `fcntl.flock` lock per genome folder, so a second process waits and then reuses the file. On Windows `fcntl` is missing and the lock is skipped (atomic rename still prevents a corrupt file). Not tested with two real processes.
- `taxon_suggest` output changed between the two days for "tuberculosis" (the listed alternatives differ), which confirms it should only ever be used to make suggestions or accept a unique token-complete match.
- NCBI's choice of reference genome can change over time; the 30 day name cache and the printed accession make this visible, and users who need reproducibility should pass a versioned accession.
