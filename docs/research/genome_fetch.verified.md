# Fact check of genome_fetch.md

Checked 2026-11-11 with live curl against api.ncbi.nlm.nih.gov and eutils.ncbi.nlm.nih.gov (generic User-Agent, no key, no contact detail), WebFetch of the two NCBI docs pages, Biopython 1.85 and pysam 0.23.0. Work files: `gf_check/` in this folder.

## Verdicts on the key claims

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| 1 | Datasets keyless limit 5 req/s, 10 with key | Confirmed | Docs page quote "By default, this rate limit is set at 5 requests per second (rps)", key gives 10 rps. Live header `x-ratelimit-limit: 5`. 12 simultaneous requests: 5 x 200, 7 x 429. |
| 2 | E-utilities keyless 3 req/s per IP, 10 with key | Confirmed | NCBI Insights 2017 post contains "more than three E-utility requests per second" and "10 requests/second" with a key. Live header `x-ratelimit-limit: 3`. 10 simultaneous einfo: 3 x 200, 7 x 429 (report saw 2 x 200; same conclusion). |
| 3 | H37Rv name -> 83332 STRAIN; both names give GCF_000195955.2 | Confirmed | taxonomy/taxon: 83332 STRAIN, 1773 SPECIES. reference_only report for both: single record GCF_000195955.2. |
| 4 | NC_000962.3 4,411,532 bp, GC 65.61 %, NCBI 65.5 | Confirmed | Biopython on the freshly downloaded zip: 4,411,532 bp, GC 65.6147 % (C 1,449,998, G 1,444,614, A 758,552, T 758,368). NCBI gc_count 2,894,612 / atgc_count 4,411,532. |
| 5 | "E. coli K-12 MG1655" not exact; suggest gives 511145 -> GCF_000005845.2 | Confirmed | HTTP 200 with `errors` "not a recognized NCBI Taxonomy name". taxon_suggest first hit 511145; reference_only report GCF_000005845.2. |
| 6 | Taxid 562 has 2 reference genomes | Confirmed | GCF_000005845.2 and GCF_000008865.2 (3 sequences). |
| 7 | SARS-CoV-2 reference_only empty `{}` HTTP 200; RefSeq complete is GCF_009858895.2; NC_045512.2 identical to MN908947.3 | Confirmed | `{}` HTTP 200 seen. Fallback filter returns 1 record, GCF_009858895.2, 29,903 bp. Both efetch sequences 29,903 bp and string-identical. |
| 8 | Download returns application/zip with *_genomic.fna and md5sum.txt; unknown accession gives HTTP 200 and invalid 971-byte zip | Confirmed | Headers `application/zip`, `filename=ncbi_dataset.zip`, 1,281,230 bytes, 5 entries as listed, fna 4,466,740 bytes, md5 matches md5sum.txt. GCF_999999999.9: HTTP 200, 971 bytes, starts with `PK`, `zipfile.is_zipfile` False. |
| 9 | efetch unknown id HTTP 400; esummary HTTP 200 "Invalid uid"; versionless NC_000962 resolves to .3 | Confirmed | All three reproduced. |
| 10 | GRCh38.p14 GCF_000001405.40, 3,099,441,038 bp, 705 sequences (24 chr, 1 MT, 680 other), 927 MB | Corrected (length) | 705 records confirmed: 25 assembled-molecule (24 + MT), 260 alt, 164 fix-patch, 90 novel-patch, 40 unlocalized, 126 unplaced. 927 MB (927.83) confirmed. But 3,099,441,038 bp is NCBI `total_sequence_length` of the primary assembly; the 705 records in the FASTA sum to 3,298,430,636 bp. Note `assigned_molecule_location_type` is "Chromosome" for 704 records, so counting chromosomes must use `role`, as the module does. Also 927 MB is the compressed zip; the extracted FASTA is about 3.3 GB on disk. |
| 11 | P. falciparum 3D7 GCF_000002765.6, 14 chr, 23,292,622 bp, GC 19.34 %; S. aureus NCTC 8325 GCF_000013425.1, 2,821,361 bp, GC 32.87 % | Confirmed | Datasets report and Biopython recount: 19.34 % and 32.87 %. NCBI gc_count/atgc_count give the same (4,505,858/23,292,622; 927,332/2,821,360, one non-ACGT base in S. aureus). |
| 12 | pysam cannot index plain gzip; filepath_index alone does not build missing index; faidx `--fai-idx` writes elsewhere | Confirmed | Reproduced in pysam 0.23.0: plain gzip raises OSError ("Cannot index files compressed with gzip, please use bgzip"), BGZF works, missing filepath_index raises OSError, `pysam.faidx(path,"--fai-idx",alt)` works and is usable. |
| 13 | Repo examples are 1-based, headers `start end` (8 regions) and `name pos` (24 positions) | Confirmed | Row counts 8 and 24. 1-based checked biologically on H37Rv: 761154-761156 = TCG (rpoB S450), 7569-7571 = GCG (gyrA A90), 2155167-2155169 = GCT, reverse complement AGC (katG S315, minus strand). The repo `--region_file` help also says 1-based. Note `primer_regions.tsv` has CRLF line endings; the parser handled them. |

Other numbers in the report re-checked live and correct: Sakai per-sequence GC (chromosome 50.54, pOSAK1 43.41, pO157 47.60), E. coli MG1655 50.79 %, SARS-CoV-2 37.97 %, influenza A (11320) 7 RefSeq complete assemblies, CDC1551 (83331) reference_only empty, unknown and superseded GCF accessions return `{}` with HTTP 200, `NC_000962.3:1-20` = TTGACCGATGACCCCGGTTC, AL123456.3 identical to NC_000962.3. NCBI `gc_percent` rounding to 0.5 is consistent with all 7 genomes checked (it is an inference from data, not documented).

## Module run on fresh inputs (genome_reference.py v2)

All cross-checked against Biopython GC and sequence slicing:

- `Mycoplasma genitalium G37` -> GCF_000027325.1, 580,076 bp, 31.69 % (Biopython 31.69 %), note about no designated reference.
- `Mycobacterium bovis BCG str. Pasteur 1173P2` -> GCF_000009445.1, 4,374,522 bp, 65.64 % (match).
- `NC_001802.1` (HIV-1) -> 9,181 bp, 42.12 % (match). `NZ_CP009100.1` and lower-case `nc_000962.3` resolve.
- `GCF_000006745.1` (Vibrio cholerae, 2 chromosomes): GC 47.49 % (NCBI 47.5). Bare `1000` correctly refused. `II:5-9`, `NC_002506:100`, `II:10..20`, last base 1,072,315 accepted; `II:0`, 1,072,316, `5-`, `-5`, `1e3` rejected. `fetch` slices equal Biopython slices.
- The strain name `Vibrio cholerae O1 biovar El Tor str. N16961` stops with 2 RefSeq candidates, as designed.

## Bugs found and fixed (genome_reference_fixed.py)

The previous `genome_reference_fixed.py` (a fix of v1 from 2026-11-10) was kept as `genome_reference_fixed_prev_v1.py`.

1. INSDC WGS/TSA contig accessions (4 letters + 8 to 10 digits, 6 letters + 9 to 11 digits, for example `AAAA02000001.1`) did not match `NUCCORE_RE`, fell through to the taxonomy path and failed with a misleading "not a recognised NCBI Taxonomy name". Fixed regex; `AAAA02000001.1` now resolves (1,710 bp).
2. `_open` slept after the last failed attempt (16 s wasted) before raising. Fixed.

## Missed or understated points

- `max_download_mb` guards the compressed zip only. Human needs about 0.93 GB for the zip plus about 3.3 GB for the extracted upper-cased FASTA, so a disk-space check is also needed.
- Human GENOME_FASTA includes 680 alt, patch and unplaced records. For primer design these create false multi-mapping; offer a primary-assembly filter (role `assembled-molecule` plus unlocalized/unplaced, excluding alt and patch).
- The rate limiter is per process. Parallel workers (16 cores) or several CLI runs share one IP and can exceed 3 or 5 req/s; network calls should stay in the main process.
- Chromosome aliases cover numeric and X/Y/Z/W names only; `chrII` (Roman numerals, yeast, Vibrio) is not an alias, though `II` works.
- The legacy repo code mixes conventions (`region_start - flank_size` used as a 0-based index, `target_length = end - start`), so carrying 1-based input into old code needs care.
- Datasets `sequence_reports` returned `total_count` 191 for 705 records; rely on `next_page_token`, not `total_count`.
- The esearch db=assembly fallback was not exercised by the module (the report says so); it remains unverified.
- The `taxon_suggest` result for "E. coli K-12 MG1655" lists "Escherichia coli B7A" second; the module's all-tokens rule correctly rejects it, but that ranking can change.

## Recommendations

All 12 recommendations are consistent with what was verified. Two additions: use `genome_reference_fixed.py`, and add a free-disk check of about 4 x the zip size before large downloads.
