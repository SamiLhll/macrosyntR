# macrosyntR
An R package to plot oxford grided plots, and results of one-tailed fischer test
in order to investigate the chromosomal evolution of whole assembled genomes.

-----------------------------------------------------------------------   

# Install

You can install this package by entering the following within R:

```{r}

devtools::install_github("SamiLhll/macroSyntR")

```

# Dependencies

This package was implemented under R version = 4.1.0
It depends on the following R packages :   
- dplyr  
- ggplot2   
- ggthemes   

# Usage

You need to calculate the orthologs between two species before considering generating any plot with this package.
I implemented a bashscript that calculates the mutual best hits of two fasta files (peptides) using blastp in my GenomicUtils repo :
[https://github.com/SamiLhll/GenomicUtils.git](https://github.com/SamiLhll/GenomicUtils.git)
When it's done you have a table with two columns, each corresponding to a species. On each row you'll have the protein IDs of two mutual best hits. 

### load MBH table into R :

The first step is to load the mutual best hit table (created with the [generate_MBH_table](https://github.com/SamiLhll/GenomicUtils/blob/a8803782f64c7ff31f0723d9e11f8f7d1a57e907/MacroSynteny/Generate_blastp_MBH) bashscript) along with the genomic coordinates of the genes coding for the proteins on their respective species in BED format.

```{r}

# most basic usage :
orthologs <- load_MBH_table("inst/extdata/sp1_vs_sp2.bed","sp1.bed","sp2.bed")
head(orthologs)

```

This command above gives the following output :


|   |sp2_pep | sp1_pep | sp1_chr | sp1_start | sp1_stop | sp2_chr | sp2_start | sp2_stop | sp1_index | sp2_index |
|---|--------|---------|---------|-----------|----------|---------|-----------|----------|-----------|-----------|
| 1 |Y1.13.p1|  T15550 |     X9  |37370904   |37389905  |    Y1   |  94074    |98562     | 1780      |   1       |
| 2 |Y1.15.p1|  T15167 |     X9  |21098862   |21127221  |    Y1   | 100166    |108471    |  1679     |    2      |
| 3 |Y1.17.p1|  T14952 |     X9  |13460779   |13475763  |    Y1   | 113555    |116613    |  1623     |    3      |
| 4 |Y1.18.p1|  T15515 |     X9  |36170511   |36209562  |    Y1   | 116920    |146404    |  1766     |    4      |
| 5 |Y1.21.p1|  T15485 |     X9  |34971016   |35007104  |    Y1   | 154209    |168520    |  1757     |    5      |
| 6 |Y1.23.p1|  T15212 |     X9  |23300128   |23313252  |    Y1   | 170863    |175425    |  1696     |    6      |

### Plot the oxford grided plot :

plot_synteny_oxfod_grid(orthologs,sp1_name="sp1",sp2_name="sp2")

