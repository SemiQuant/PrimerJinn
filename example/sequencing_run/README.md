# Example read-count table

`counts.tsv` is what `primerjinn-multiplex --sequencing` (and the GUI's "From sequencing" section) reads: one row per amplicon with the number of reads that are its designed product, and a row named `unaligned` (and here also `off_target`) with the reads that match no designed amplicon. The numbers are from one nanopore run of a 23-amplicon M. tuberculosis panel; that panel was not designed by primerJinn, so the amplicon names are gene names rather than primerJinn target ids. To try the feedback step, name the rows after your own design's target ids (or primer names without `_F`/`_R`).

See the README section "After sequencing" for the format.
