# macrosyntR

<!-- badges: start -->
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->


An R package to calculate, automatically order and visualize the macro-synteny blocks between two species.

-----------------------------------------------------------------------   
<img src="https://github.com/SamiLhll/macrosyntR/blob/55a9b126baee3e88335dd8c834eb5e70b7863aec/inst/img/oxford_grid.png" alt="oxford_grid" width="800"/>


# Installation

```{r}

# Get the development version from GitHub:
# install.packages("devtools")
devtools::install_github("SamiLhll/macrosyntR")

```


# Getting started :

### Foreword : 

This package won't calculate the orthologs from the sequences. As you will have to compute it yourself, I need to clarify that each dot on the oxford grid plot corrspond to one ortholog. It wouldn't work to use orthogroups as defined by OrthoFinder() as it would pack together more than one gene per species. 
Instead, for linux users I would recommend that you use my shellscript [mbhXpress](https://github.com/SamiLhll/mbhXpress)
It uses diamond blast so it doesn't take more than few minutes to calculate mutual best hits with a couple of CPUs. 

### Overview :

This package provides the user with 5 functions :

|     Function          |         description                                                                                          | 
|-----------------------|--------------------------------------------------------------------------------------------------------------|
| load_mbh_df()         | integrates genomic coordinates (bedfiles) of the orthologs of the two species to compare                     |
| calculate_macrosynt() | calculates all the chromosomes to each other and calculate the significantly conserved macro synteny blocks  |
| reorder_synteny()     | takes an mbh_table (from load_mbh_df()) and outputs an mbh_table ordered for plot_oxford_grid()              |
| plot_macrosynt()      | draws a dotplot displaying the significant macro-synteny blocks with their relative amount of orthologs      |
| plot_oxford_grid()    | draws an oxford grid from an mbh_table (output of either load_mbh_df() or reorder_synteny()                  |    


# Step-by-step tutorial with publicly available data :

### 1 - Generate the input files :

The first step is to load the mutual best hit table (that will have been created with the [generate_MBH_table](https://github.com/SamiLhll/GenomicUtils/blob/a8803782f64c7ff31f0723d9e11f8f7d1a57e907/MacroSynteny/Generate_blastp_MBH) bashscript) along with the genomic coordinates of the genes coding for the proteins on their respective species in BED format.

### 2 - Draw an oxford grid to visualize my data :

### 3 - Calculate and plot the significant macrosynteny blocks :

### 4 - Put some order :

### 5 - Summary :

```{r}

# most basic usage :
orthologs <- load_MBH_table("inst/extdata/sp1_vs_sp2.tab","sp1.bed","sp2.bed")
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

```{r}
# most basic usage :
plot_synteny_oxfod_grid(orthologs,sp1_name="sp1",
                                  sp2_name="sp2")

```

![](inst/img/Rplot1.png)


```{r}
# filter to keep only the meaningfull associations using the sp2_keep_chr_names argument :

plot_synteny_oxford_grid(orthologs,sp1_name = "sp1",
                                   sp2_name = "sp2",
                                   sp2_keep_chr_names = c("Y1","Y3","Y7","Y9","Y14","Y19"))

```

![](inst/img/Rplot2.png)


```{r}

# Play with other rendering arguments to 
# reorder the chromosomes using the sp2_chr_order, 
# add some colors by setting the colors argument to TRUE
# or change the size of the dots : 

plot_synteny_oxford_grid(orthologs,sp1_name = "sp1",
                                   sp2_name = "sp2",
                                   sp2_keep_chr_names = c("Y1","Y3","Y7","Y9","Y14","Y19"),
                                   colors =TRUE,
                                   sp2_chr_order = c("Y3","Y7","Y14","Y19","Y9","Y1"),
                                   dot_size = 0.6)

```

![](inst/img/Rplot3.png)

### Calculate contingency table and plot result of fischer test :

```{r}

# calculate contingency table :

contingency_table <- calculate_contingency_table(orthologs)

# most basic usage of plotting function :
plot_fischer_test(contingency_table,sp1_name = "sp1",
                                    sp2_name = "sp2)

```
![](inst/img/Rplot4.png)

```{r}
# use arguments to filter and reorder the chromosomes on the plot :
plot_fischer_test(contingency_table,sp1_name = "sp1",
                                    sp2_name = "sp2,
                                    sp2_keep_chr_names = c("Y3","Y7","Y14","Y19","Y9","Y1"),
                                    sp2_chr_order = c("Y3","Y7","Y14","Y19","Y9","Y1"))

```

![](inst/img/Rplot5.png)

