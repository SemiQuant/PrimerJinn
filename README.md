<div>
    <img src="https://user-images.githubusercontent.com/8179171/236663567-94d1f5dc-2ac6-49de-9fc1-a99c7a13945d.png" width="20%" height="20%">
    <p>primerJinn has two main functions, it designs primers for a multiplex PCR of given target regions of a DNA sequence in a FASTA file and it performs in silico PCR given a list of primers and a reference FASTA file.</p>
</div>

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7903629.svg)](https://doi.org/10.5281/zenodo.7903629) -->


The scripts takes several arguments:

| **Program**       | **Required** | **Parameter**      | **Description**                                                                                                                                                           | **Default**                             |
|-------------------|--------------|--------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------|
| **Both**          | Y            | input_file         | The path to the input FASTA file.                                                                                                                                         | NA                                      |
| **primer design** | Y            | region_file        | The path to the primer design region file. This file should have three columns: name, start position, and end position (1-based). It can be in either tsv or xlxs format. | NA                                      |
| **Both**          | N            | salt_concentration | Salt concentration (in nM, Ignored if Q5 True).                                                                                                                           | 50                                      |
| **Both**          | N            | target_tm          | The desired melting temperature (Tm) for the primers                                                                                                                      | 60C                                     |
| **Both**          | N            | output             | The name of the output file.                                                                                                                                              | 'MultiPlexPrimerSet' or 'in_silico_PCR' |
| **primer design** | N            | primer_len         | The desired length of the primers.                                                                                                                                        | 20                                      |
| **primer design** | N            | product_size_min   | The desired min size for the PCR product.                                                                                                                                 | 400                                     |
| **primer design** | N            | product_size_max   | The desired max size for the PCR product.                                                                                                                                 | 800                                     |
| **primer design** | N            | ret                | The maximum number of primer pairs to return                                                                                                                              | 100                                     |
| **primer design** | N            | Q5                 | A boolean indicating whether to use NEB Q5 hotstart polymerase settings for primer3                                                                                       | TRUE                                    |
| **primer design** | N            | background         | The path to the mispriming library FASTA file                                                                                                                             |                                         |
| **primer design** | N            | ill_adapt         | Add Illumina partial adapters                                                                                                                              | FALSE                                         |
| **primer design** | N            | clamp         | Require GC clamp                                                                                                                               | 0                                         |
| **primer design** | N            | poly         | Maximum allowable length of a mononucleotide repeat (poly-X) in the primer sequence                                                                                                                               | 3                                         |
| **primer design** | N            | no_self_background         | If specified, primers are NOT checked for mispriming against the input genome as for large genomes this is very slow                                                                                                                               | False                                         |
| **in silico PCR** | N            | product_size_max   | Maximum length of PCR products in nucleotides.                                                                                                                            | 2000                                    |
| **in silico PCR** | N            | req_five           | Require the 5' end of the primer to bind?                                                                                                                                 | TRUE                                    |


### Dependencies

-   Python 3
-   [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK569861/)

### Installation

```
pip install primerJinn
```

### Multiplex primer generation example usage

```
getMultiPrimerSet \
    --region_file "./example/primer_regions.tsv" \
    --input_file "./example/ref.fasta" \
    --target_tm 65 \
    --primer_len 20 \
    --product_size_min 400 \
    --product_size_max 800 \
    --ret 100 \
    --Q5 \
    --background "" \
    --output "example"
```

###MultiPlexPrimerSet.xlxs

| Forward Primer       | Reverse Primer       | Forward tm | Reverse tm | Product Size | Name                   | Forward Primer TM NEB | Reverse Primer TM NEB |
|----------------------|----------------------|------------|------------|--------------|------------------------|-----------------------|-----------------------|
| GAGGAACACCACTAGTACCG | CTCGATGACTTTACGGCCAT | 65.08334   | 64.69482   | 675          | Target 1461045-1461291 | 64                    | 64                    |
| CATGGGATATGGAGCGATCG | GGGGTCGTAGGAGATCTTGA | 65.95586   | 65.59987   | 791          | Target 490900-491416   | 65                    | 65                    |
| CCGGTTGTCCATTCCGTTTA | CTGTACGTATTTGGGTTGCG | 64.17967   | 65.57867   | 454          | Target 1303831-1303911 | 65                    | 64                    |
| GGATGCGAGCTATATCTCCG | AATACGCCGAGATGTGGATG | 65.15669   | 64.9819    | 458          | Target 1674182-1674222 | 65                    | 65                    |
| CAACAGTTCATCCCGGTTCG | GACGGATTTGTCGCTCACTA | 64.88889   | 66.57304   | 759          | Target 2288681-2289242 | 66                    | 64                    |
| GCCACCATCGAATATCTGGT | GCTCCAGGAAGGGAATCATC | 65.63736   | 65.01728   | 778          | Target 761007-761277   | 64                    | 65                    |
| GCATACCGAACGTCACAGAT | ACGGTCACCTACAAAAACGG | 65.68941   | 65.23466   | 665          | Target 778990-779488   | 65                    | 65                    |
| GCTCTTAAGGCTGGCAATCT | CGGTCACACTTTCGGTAAGA | 64.82047   | 65.53067   | 577          | Target 2154831-2154873 | 65                    | 64                    |


There is an online version of primerJinn getMultiPrimerSet available at [DrDx.Me](https://drdx.ucsf.edu/)


#### Notes
When testing the primers for a targeted diagnostic assay, we would recommend doing a qPCR using individual primers and reagents intended for the PCR, and adding in EVA Green plus (20X in water) and ROX (50X). This will allow you see the relative efficiency of each primer, and investigate the melt curves, the products can also be run on an agarose gel to confirm single bands.


### PCR in silico example usage

```
PCRinSilico \
   --primer_seq ./example/primers.txt \
   --target_tm 50 \
   --input_file ./example/ref.fasta
```

#### in_silico_PCR.tsv
| qseq1 | qseq1_input          | qstart1 | qend1 | direction1 | mismatch1 | qseq2 | qseq2_input          | qstart2 | qend2 | direction2 | mismatch2 | binding_pos_diff | reference   | ref_region |         |
|-------|----------------------|---------|-------|------------|-----------|-------|----------------------|---------|-------|------------|-----------|------------------|-------------|------------|---------|
| p1    | GAGGAACACCACTAGTACCG | 1       | 20    | +          | 0         | p9    | CTCGATGACTTTACGGCCAT | 1       | 20    | -          | 0         | 655              | NC_000962.3 | 1461637    | 1460982 |
| p2    | CATGGGATATGGAGCGATCG | 1       | 20    | +          | 0         | p10   | GGGGTCGTAGGAGATCTTGA | 1       | 20    | -          | 0         | 771              | NC_000962.3 | 491474     | 490703  |
| p3    | CCGGTTGTCCATTCCGTTTA | 1       | 20    | +          | 0         | p11   | CTGTACGTATTTGGGTTGCG | 1       | 20    | -          | 0         | 434              | NC_000962.3 | 1304111    | 1303677 |
| p4    | GGATGCGAGCTATATCTCCG | 1       | 20    | +          | 0         | p12   | AATACGCCGAGATGTGGATG | 1       | 20    | -          | 0         | 438              | NC_000962.3 | 1674558    | 1674120 |
| p5    | CAACAGTTCATCCCGGTTCG | 1       | 20    | +          | 0         | p13   | GACGGATTTGTCGCTCACTA | 1       | 20    | -          | 0         | 739              | NC_000962.3 | 2289391    | 2288652 |
| p6    | GCCACCATCGAATATCTGGT | 1       | 20    | +          | 0         | p14   | GCTCCAGGAAGGGAATCATC | 1       | 20    | -          | 0         | 758              | NC_000962.3 | 761564     | 760806  |
| p7    | GCATACCGAACGTCACAGAT | 1       | 20    | +          | 0         | p15   | ACGGTCACCTACAAAAACGG | 1       | 20    | -          | 0         | 645              | NC_000962.3 | 779595     | 778950  |
| p8    | GCTCTTAAGGCTGGCAATCT | 1       | 20    | +          | 0         | p16   | CGGTCACACTTTCGGTAAGA | 1       | 20    | -          | 0         | 557              | NC_000962.3 | 2155289    | 2154732 |


#### in_silico_PCR_amplicon_interactions.tsv
| amplicon1_PF | amplicon1_PR | amplicon2_PF | amplicon2_PR | tm          |
|--------------|--------------|--------------|--------------|-------------|
| kkd_F_2      | sge_R        | kkd_R        | sge_R        | 90.22726832 |


#### in_silico_PCR_primer_dimears.tsv
| Sequence1 | Sequence2 | MeltingTemp |
|-----------|-----------|-------------|
| p1        | p2        | 73          |
| p1        | p3        | 74          |
| p1        | p4        | 73          |
| p1        | p5        | 73          |
| p1        | p6        | 73          |
| p1        | p7        | 72          |
| p1        | p8        | 73          |
| p1        | p9        | 73          |
| p1        | p10       | 74          |


# Citation
(Limberis, J.D., Metcalfe, J.Z. primerJinn: a tool for rationally designing multiplex PCR primer sets for amplicon sequencing and performing in silico PCR. BMC Bioinformatics 24, 468 (2023).)[https://doi.org/10.1186/s12859-023-05609-1]