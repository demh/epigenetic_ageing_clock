# Biological insights into the epigenetic ageing clock

Created by Daniel E. Martin-Herranz.

Copyright (C) 2019 D.E. Martin-Herranz

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2559588.svg)](https://doi.org/10.5281/zenodo.2559588)


## Introduction

[Epigenetic clocks](https://www.nature.com/articles/s41576-018-0004-3) are mathematical models that predict the biological age of an individual using DNA methylation data, and which have emerged in the last few years as the most accurate biomarkers of the ageing process. However, little is known about the molecular mechanisms that control the rate of such clocks. In this work, we have studied the behaviour of the human epigenetic clock in patients with a variety of developmental disorders, harbouring mutations in proteins of the epigenetic machinery. We found a role for the H3K36 methylation machinery in the control of the ticking rate of the epigenetic clock.  


## Analyses related to the epigenetic ageing clock

This repository contains all the code that you need to run the following analyses:

* Pre-processing raw DNA methylation data (IDAT files, 450K array) to obtain beta-values.

* Calculating the biological (epigenetic) age of your samples using [Horvath's epigenetic clock](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115).

* Estimating blood cell type composition using the [Houseman method](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-86) and the [IDOL reference](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0943-7). 

* Correcting for batch effects in the context of the epigenetic clock (using the control probes from the 450K array).

* Calculating the epigenetic age acceleration (EAA) of your samples (i.e. whether their biological age differs from their chronological age). 

* Calculating the mitotic age of your samples according to the [epigenetic mitotic clock](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1064-3).

* Calculating methylation Shannon entropy.  

* Identifying differentially methylated probes (DMPs) during ageing and calculating enrichment for certain (epi)genomic features.

* Reproducing all the analyses and the figures in our manuscript.


## Quick workflow to calculate DNAmAge (biological age from Horvath's clock) starting with raw IDAT files.

1. Obtain the matrix of beta-values from raw IDAT files (IDAT files --> idat_to_DNAmAge_methylation_matrix.R --> beta-values matrix):

`Rscript ~/epigenetic_ageing_clock/utils/idat_to_DNAmAge_methylation_matrix.R ~/path/to/idats ~/path/to/output ~/epigenetic_ageing_clock/utils/  'noob'`

2. Calculate DNAmAge (beta-values matrix --> calculate_DNAm_age.R --> DNAmAge calculations):

`Rscript ~/epigenetic_ageing_clock/utils/calculate_DNAm_age.R -i ~/path/to/output/DNAmAge_methylation_matrix_from_idat_xxxxx.csv -o ~/path/to/output/ -a ~/epigenetic_ageing_clock/utils/`


## Citing us 

If you found this code useful for you research, please cite our original publication:

> Martin-Herranz, D.E. et al. Screening for genes that accelerate the epigenetic aging clock in humans reveals a role for the H3K36 methyltransferase NSD1. Genome Biol 20, 146 (2019). https://doi.org/10.1186/s13059-019-1753-9


## Contacting us

If you experience any issues with the code or have any suggestions, please contact us at daniel@chronomics.com.


