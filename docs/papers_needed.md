# Papers the research could not read in full

Update 2026-11-11: the user supplied all of these, plus the SADDLE, Pythia and Owczarzy 2008 supplements (`/Users/semiquant/Downloads/papers`, not copied into the repo). Read and applied: Owczarzy 2008 + SI (out-of-sample test), SantaLucia and Hicks 2004 (dangling ends), SantaLucia 2007, SADDLE + SI S9, Pythia, Su 1996, Musso 2006, PRISM. Still to read: Huang 1992, Lefever 2013, Stadhouders 2010, Rejali 2018 (Taq 3'-mismatch factors). See DESIGN.md for what changed.

The research agents worked from abstracts, code reimplementations, patents or secondary sources for these. Full texts (and supplements) would let us check or calibrate the numbers named in the last column. Put PDFs in `docs/papers/` (not committed) and note them in `PROGRESS.md`.

| Priority | Paper | DOI | What we would check |
|---|---|---|---|
| 1 | Owczarzy et al. 2008, Biochemistry 47:5336, Mg2+ salt correction, **with supplementary data** | 10.1021/bi702363u | Coefficients were read from the authors' patent (US 8,055,451). The supplement's KCl + Mg2+ melting data would give an out-of-sample accuracy test in PCR-like buffer (the in-sample MAE is 1.25 C). Also the Mg2+:dNTP binding constant |
| 1 | SantaLucia 2007, Methods Mol Biol 402:3-34, "Physical principles and Visual-OMP software for optimal PCR design" | 10.1007/978-1-59745-528-2_1 | The percent-bound targets and multi-state model behind designing on fraction bound at Ta rather than Tm. Our theta and dG_eff thresholds have no primary anchor yet |
| 1 | Xie et al. 2022, Nat Commun 13:1881, SADDLE, **Supplementary section S9** | 10.1038/s41467-022-29500-4 | Temperature and salt behind their -11.5 kcal/mol binding target and their dimer data, to transfer their thresholds to our model |
| 2 | Mann et al. 2009, NAR 37:e95, Pythia, **supplementary data** | 10.1093/nar/gkp443 | Production threshold and simulation conditions for equilibrium primer design |
| 2 | Su et al. 1996, NAR 24:1574, extension temperature for A+T-rich DNA | 10.1093/nar/24.8.1574 | Basis of the 60 to 65 C extension advice for AT-rich genomes such as P. falciparum (only seen through a secondary source) |
| 2 | Huang, Arnheim and Goodman 1992, NAR 20:4567, Taq extension of 3' mismatches | 10.1093/nar/20.17.4567 | The Taq mismatch factors in the specificity model (read from the abstract only) |
| 2 | Lefever et al. 2013, Clin Chem 59:1470 | 10.1373/clinchem.2013.203653 | How many 3' mismatches block priming (abstract only); sets the off-target score cut-offs |
| 2 | Rejali, Moric and Wittwer 2018, Clin Chem 64:801 | 10.1373/clinchem.2017.282285 | Rate loss from mismatches up to 6 nt from the 3' end (abstract only) |
| 3 | Stadhouders et al. 2010, J Mol Diagn 12:109 | 10.2353/jmoldx.2010.090035 | Mismatch effects in qPCR (abstract only) |
| 3 | Musso et al. 2006, J Mol Diagn 8:544, betaine + DMSO + 7-deaza-dGTP | 10.2353/jmoldx.2006.060058 | Exact additive concentrations behind the GC-rich advice |
| 3 | Chen et al. 2026, Bioinformatics 42(7), PRISM multiplex design | 10.1093/bioinformatics/btag478 | Algorithm details beyond the abstract, for comparison |
| 3 | SantaLucia and Hicks 2004, Annu Rev Biophys Biomol Struct 33:415 | 10.1146/annurev.biophys.32.110601.141800 | NN table was checked against three code bases and SantaLucia 1998; the original would close the loop |
