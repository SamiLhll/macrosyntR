# macrosyntR <a><img src='https://github.com/SamiLhll/macrosyntR/blob/5db3512442b0ca271b69dc34dc5989190f0fbcc0/inst/img/macrosyntR.png' align="right" height="190" /></a>


An R package for evaluation of pair-wise synteny conservation at the genome-wide scale.
It takes a table of orthologs and genome annotation files formatted as BED to automatically
infer significantly conserved linkage groups, and order them on an oxford grid.   

<!-- badges: start -->
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

-----------------------------------------------------------------------   


# Installation

```{r}

# Before this package is submitted to CRAN,
# get the development version from GitHub using devtools :
# install.packages("devtools")
devtools::install_github("SamiLhll/macrosyntR")

```


# Usage

Check out the vignette for a comprehensive step-by-step tutorial illustrating how the package works using publicly available data, and how to customize the analysis. 

### Preparing input data :

To start drawing an oxford grid using this package, you'll need three files :   
* 1 - A two columns table of orthologous pairs of genes between species to compare (as generated by [mbhxpress](https://github.com/SamiLhll/mbhXpress))   
* 2 - A bed file listing the genomic coordinates and sequence names of all the orthologs of the 1st species (with names matching the first column of file 1)   
* 3 - A bed file, same as the other, for the 2nd species. (with names matching those of the first column of file 2)

### Get an automatically ordered and colored oxford grid :

To illustrate the results of the package we compare the publicly available data from the lancelet Branchiostoma floridae ([Simakov et al. 2020](https://doi.org/10.1038/s41559-020-1156-z)) with the Siboglinidae Paraescarpia echinospica ([Sun et al. 2021](https://doi.org/10.1093/molbev/msab203))   
Once you have your pairs of orthologs, getting an ordered oxford grid using this package is achieved as following :   

```{r}

library(macrosyntR)

# load table of mutual best hits and integrate genomic coordinates :
my_orthologs_table <- load_mbh_df(mbh_table = "DATA/Bflo_vs_Pech.tab",
                                  sp1_bed = "DATA/Bfloridae.protein_products.bed",
                                  sp2_bed = "DATA/Pechinospica.protein_products.bed")

# draw an oxford grid with auto-clustering and coloring of the linkage groups :
plot_oxford_grid(mbh_df = my_orthologs_table,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica",
                 auto_order_clusters = TRUE,
                 color_clusters = TRUE)

```

<a><img src='https://github.com/SamiLhll/macrosyntR/blob/1c6dd09793e41bcbe2850669614cef3905925b3c/inst/img/Bflo_Pech_oxf_grid.png' align="center" height="400" /></a>   


# Getting help

Need help, Identified a bug, or want to see other features implemented ?   
Feel free to open an issue here or send an email at :   
elhilali.sami@gmail.com

