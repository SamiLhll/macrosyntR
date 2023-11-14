# macrosyntR <a><img src='https://github.com/SamiLhll/macrosyntR/blob/f7c23587e2ae4c08b3709d8ac046128497a8fe60/inst/img/macrosyntR_logo.150pp.png' align="right" height="190" /></a>


An R package for evaluation of synteny conservation at the genome-wide scale.
It takes a table of orthologs and genome annotation files formatted as BED to automatically
infer significantly conserved linkage groups, and order them on an Oxford grid or a chord diagram using a network based greedy algorithm.   

<!-- badges: start -->
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
  [![CRAN](http://www.r-pkg.org/badges/version/macrosyntR)](https://cran.r-project.org/package=macrosyntR) 
  [![downloads](https://cranlogs.r-pkg.org/badges/grand-total/macrosyntR)](https://cran.r-project.org/package=macrosyntR)
  [![bioRxiv:10.1101/2023.01.26.525673](https://img.shields.io/badge/bioRxiv-10.1101/2023.01.26.525673-abcfed.svg)](https://doi.org/10.1101/2023.01.26.525673)
<!-- badges: end -->

-----------------------------------------------------------------------   


# Installation

```{r}

# A stable version is available on CRAN and can be downloaded using :
install.packages("macrosyntR")
# get the development version from GitHub using devtools :
# install.packages("devtools")
devtools::install_github("SamiLhll/macrosyntR",build_vignettes = TRUE)
# building the vignette makes the installation a bit longer but its mandatory so ou can access it by doing :   
vignette("macrosyntR")

```


# Usage

Check out the vignette for a comprehensive step-by-step tutorial illustrating how the package works using publicly available data, and how to customize the analysis. 

### Preparing input data :

To start drawing an oxford grid using this package, you'll need three files :   
* 1 - A two columns table of orthologous pairs of genes between species to compare (as generated by [rbhxpress](https://github.com/SamiLhll/rbhXpress))   
* 2 - A bed file listing the genomic coordinates and sequence names of all the orthologs of the 1st species (with names matching the first column of file 1)   
* 3 - A bed file, same as the other, for the 2nd species. (with names matching those of the first column of file 2)

### Get an automatically ordered and colored Oxford grid :

To illustrate the results of the package we compare the publicly available data from the lancelet Branchiostoma floridae ([Simakov et al. 2020](https://doi.org/10.1038/s41559-020-1156-z)) with the Siboglinidae Paraescarpia echinospica ([Sun et al. 2021](https://doi.org/10.1093/molbev/msab203))   
Once you have your pairs of orthologs, getting an ordered Oxford grid using this package is achieved as following :   

```{r}

library(macrosyntR)

# Load table of orthologs and integrate with genomic coordinates :
my_orthologs <- load_orthologs(orthologs_table = "DATA/Bflo_vs_Pech.tab",
                               sp1_bed = "DATA/Bfloridae.protein_products.bed",
                               sp2_bed = "DATA/Pechinospica.protein_products.bed")

# Draw an oxford grid :
p1 <- plot_oxford_grid(my_orthologs,
                       sp1_label = "B.floridae",
                       sp2_label = "P.echinospica")
p1

# Automatically reorder the Oxford grid and color the detected clusters (communities):
p2 <- plot_oxford_grid(my_orthologs,
                       sp1_label = "B.floridae",
                       sp2_label = "P.echinospica",
                       reorder = TRUE,
                       color_by = "clust")
p2

# Plot the significant linkage groups :
my_macrosynteny <- compute_macrosynteny(my_orthologs)
p3 <- plot_macrosynteny(my_macrosynteny)
p3


# Call the reordering function, test significance and plot it :
my_orthologs_reordered <- reorder_macrosynteny(my_orthologs)
my_macrosynteny <- compute_macrosynteny(my_orthologs_reordered)
p4 <- plot_macrosynteny(my_macrosynteny)
p4

```

<a><img src='https://github.com/SamiLhll/macrosyntR/blob/a5f008b3b596f6ec1c4a73952fde7bf3fbdad57c/inst/img/example.png' align="center" height="650" /></a>   
# Additional ressources

### Computing orthologs as reciprocal best hits :
If you don't know how to do it, I masde an util called rbhXpress that uses Diamond blast to generate an output compatible with macrosyntR.
Please find more details in the following repository : [rbhXpress](https://github.com/SamiLhll/rbhXpress)

### Using OrthoFinder output : 
I implemented a shell script that extracts all single copy orthologs from a run of OrthoFinder and turns it into pairwise tables for all the species. Please find more details in the following repository : [OF_to_macrosyntR](https://github.com/SamiLhll/OF_to_macrosyntR)

# Getting help

Need help, Identified a bug, or want to see other features implemented ?   
Feel free to open an issue here or send an email at :   
elhilali.sami@gmail.com

# Citation

If used in your research, please cite :   

* El Hilali, S., Copley R., "macrosyntR : Drawing automatically ordered Oxford Grids from standard genomic files in R", bioRxiv (2023). [doi:10.1101/2023.01.26.525673](https://doi.org/10.1101/2023.01.26.525673)

