# From Genotypes to Haplotypes: A Reproducible Haplotype Workflow for investigating Genetic Structures Beyond Variants in Human Population Genetics

This repository consists of the companion code for the posted project "From Genotypes to Haplotypes: A Reproducible Framework for Investigating Genetic Structures Beyond Variants in Human Population Genetics"

The workflow initiates from the phased 1000 Genotypes for a set of selected variants, builds haplotypes per locus and consequently investigates:
1. Haplotype frequencies across populations and superpopulations.
2. Relationships among haplotype sequences, through utilization of Hamming Distance and PCoA.
3. Population structure with PCA & DAPC.
4. Co-occurrence events of neighboring variant alleles (here: rs7782412 & rs7799652)

The pipeline examines three distinct loci, one mapped to FOXP2 and two to CNTNAP2 coordinates. The analysis is descriptive, resulting in frequency tables, distances and it does not include statistical tests.

### Requirements

- R (developed with R 4.5)
- R packages: vcfR, ggplot2, reshape2, tidyr, tidyverse, dplyr, randomcoloR, adegenet, ade4, hierfstat, stringdist, ape, plotly, GWLD, pheatmap


### Input Data

| File | Content |
|---|---|
| `contr_phased.vcf.gz` | Phased genotypes (`GT`) for the selected variants, for all 2,504 individuals. Variant IDs must be unique (the script checks this). |
| `contrib_ids.txt` | One variant ID per line: the variants that were previously identified as contributing to population differentiation. |
| `pop.csv` | Comma-separated, no header. Column 2 is the sample ID and column 3 is the population code (for example `YRI`, `CEU`). Column 1 is not used. |

Note: In the uploaded workflow files are downloaded by releases corresponding to repository "FOXP2_CNTNAP2_Population_Analysis" (https://github.com/evalysiova/FOXP2_CNTNAP2_Population_Analysis/).

 ### Loci tested

Coordinates tested are mapped on Chromosome 7 (GRCh38). The analysis uses each locus extended by 20kb on both sides.


| Region | Locus | Coordinates | With 20 kb margin |
|---|---|---|---|
| 1 | *FOXP2* | 7:114362480-114650436 | 7:114342480-114670436 |
| 2 | *CNTNAP2* | 7:147017367-147323381 | 7:146997367-147343381 |
| 3 | *CNTNAP2* | 7:148057724-148333434 | 7:148037724-148353434 |

### Outputs and Figures

Figures are written to `output_haplo/figures/`:
 
- `<region><superpopulation>_haplotypes.png` (for example `region1EAS_haplotypes.png`): stacked bar plots of haplotype frequencies per population
- `haplo_analysis.pdf`: the base R and heatmap plots created after the PDF device is opened
| Poster figure | Produced by |
|---|---|
| Figure 1, haplotype frequencies per superpopulation | Step 4 (PNG files) |
| Figure 2A, haplotype frequencies clustered by sequence distance | Step 8 (`pheatmap`) |
| Figure 2B, DAPC by population | Step 5 (`dapc3`) |
| Figure 2C, PCoA of haplotype sequences | Step 7 (2D and 3D plots) |
| Figure 3B, two-variant haplotypes | Step 9 |
| Figure 3A, linkage (RMI) heatmaps | Not produced by this script |
 

### Authors and Funding

Evangelia Sotiria Lysiova, Marianthi Georgitsi and Aristotelis C. Papageorgiou
Department of Molecular Biology and Genetics, Democritus University of Thrace, Alexandroupolis, Greece.
 
This work is partially funded by the MSc Programme "Applied Bioinformatics and Data Analysis" of the Department of Molecular Biology and Genetics (D.U.Th.).

