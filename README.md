# macrosyntR

<!-- badges: start -->
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->


An R package to identify, automatically order and visualize the macro-synteny blocks between two species.

-----------------------------------------------------------------------   
<img src="https://github.com/SamiLhll/macrosyntR/blob/3eaae713fe29fd7e9b7894d77ded4007022d223c/inst/img/oxford_grid.png" alt="oxford_grid" width="800"/>


# Installation

```{r}

# Get the development version from GitHub:
# install.packages("devtools")
devtools::install_github("SamiLhll/macrosyntR")

# remove the package :
remove.packages("macrosyntR")

```


# Getting started :

### Foreword : 

This package won't calculate the orthologs from the sequences. As you will have to compute it yourself, I need to clarify that each dot on the oxford grid correspond to one ortholog. It wouldn't work to use orthogroups as defined by OrthoFinder() as it would pack together more than one gene per species.   
Instead, for linux users I would recommend that you run [mbhXpress](https://github.com/SamiLhll/mbhXpress).   
It uses diamond blast so it doesn't take more than few minutes to calculate mutual best hits with a couple of CPUs. 

### Overview :

|     Function          |         description                                                                                          | 
|-----------------------|--------------------------------------------------------------------------------------------------------------|
| load_mbh_df()         | integrates genomic coordinates (bedfiles) of the orthologs of the two species to compare                     |
| calculate_macrosynt() | compares all the chromosomes to each other and identifies the significantly conserved macro synteny blocks   |
| reorder_synteny()     | takes an mbh_table (from load_mbh_df()) and outputs an mbh_table ordered for plot_oxford_grid()              |
| plot_macrosynt()      | draws a dotplot displaying the significant macro-synteny blocks with their relative amount of orthologs      |
| plot_oxford_grid()    | draws an oxford grid from an mbh_table (output of either load_mbh_df() or reorder_synteny()                  |    


# Step-by-step tutorial using publicly available data :

### 0 - download and pre-process the data :

I'm going to show how this package can help visualizing the macro-synteny blocks by comparing the data of the lancelet Branchiostoma floridae ([Simakov et al. 2020](https://doi.org/10.1038/s41559-020-1156-z)) with the data of the vestimentifera (deep-sea tubeworm) Paraescarpia echinospica ([Sun et al. 2021](https://doi.org/10.1093/molbev/msab203)).   

Download the sequences of proteins (fasta format) and their genomic coordinates :    

 - B.floridae : 
 The data are available on ncbi at https://www.ncbi.nlm.nih.gov/genome/?term=txid7739   
 get the protein sequences at  by clicking the "Download sequences in FASTA format for protein".   
 get the genomic coordinates by clicking "Download genome annotation in tabular format" and further click download as csv.   
 
 - P.echinospica :   
 The data are available on figshare.   
 get the protein sequences here : https://figshare.com/ndownloader/files/28945467   
 get the genome annotation here : https://figshare.com/ndownloader/files/28945458   
 
 
 Compute the mutual best hits of the fasta sequences. Using [mbhXpress](https://github.com/SamiLhll/mbhXpress) you can achieve it by typing the following in your terminal :
 
 ```{bash,eval = FALSE}
 
 # call mbhXpress with using 6 threads :
 bash mbhXpress -a GCF_000003815.2_Bfl_VNyyK_protein.faa -b Pec_ragoo_v1.0.pep.fasta -o Bflo_vs_Pech.tab -t 6
 
 ```
 
To convert the genome annotation to the [bed file format](https://www.ensembl.org/info/website/upload/bed.html), I'm using the following command lines (if unfamiliar with this you can use a spreadsheet software). The concept is to keep the chrom, chromStart, chromEnd mandatory fields plus the name optional field that links the genomic region with the protein product :   
 
  ```{bash, eval = FALSE}
 
 # B.floridae CSV file to bed
tail -n +2 proteins_75_971579.csv | cut -d "," -f1,3,4,9  | sed -e 's/\"//g' -e 's/,/\t/g' -e 's/chromosome /BFL/g' > Bfloridae.protein_products.bed

 # P.echinospica gff file to bed
fgrep "gene" Pec_genes.ragoo_v1.0.gff | cut -f1,4,5,9 | cut -d ";" -f 1 | fgrep "Superscaffold" | sed -e 's/ID=//g' -e 's/Superscaffold/PEC/g' > Pechinospica.protein_products.bed
 
 ```
 
### 1 - Draw an oxford grid to visualize the data :

### 2 - Calculate and plot the significant macrosynteny blocks :

### 3 - Put some order into it :

### 4 - Summary :


