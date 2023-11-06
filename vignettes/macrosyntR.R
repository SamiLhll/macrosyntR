## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(macrosyntR)

## -----------------------------------------------------------------------------

my_orthologs_table <- load_orthologs(orthologs_table = system.file("extdata","Bflo_vs_Pech.tab",package="macrosyntR"),
                                     bedfiles = c(system.file("extdata","Bflo.bed",package="macrosyntR"),
                                     system.file("extdata","Pech.bed",package="macrosyntR")))

head(my_orthologs_table)


## -----------------------------------------------------------------------------

my_orthologs_with_3_sp <- load_orthologs(orthologs_table = system.file("extdata","Single_copy_orthologs.tsv",package="macrosyntR"),
                                     bedfiles = c(system.file("extdata","Bflo.bed",package="macrosyntR"),
                                                  system.file("extdata","Pech.bed",package="macrosyntR"),
                                                  system.file("extdata","Pyes.bed",package="macrosyntR")))

head(my_orthologs_with_3_sp)


## -----------------------------------------------------------------------------

# compute significance :
macrosynteny_df <- compute_macrosynteny(my_orthologs_table)
head(macrosynteny_df)


## ----eval = FALSE, out.width = '450px'----------------------------------------
#  
#  # visualize the loaded data on a oxford grid :
#  plot_oxford_grid(my_orthologs_table,
#                   sp1_label = "B.floridae",
#                   sp2_label = "P.echinospica")
#  
#  # Visualize the results of the test of significance :
#  plot_macrosynteny(macrosynteny_df,
#                    sp1_label = "B.floridae",
#                    sp2_label = "P.echinospica")
#  

## ----echo = FALSE,out.width = c('300px','300px')------------------------------

# visualize the loaded data on a oxford grid :
plot_oxford_grid(my_orthologs_table,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica")

# Visualize the results of the test of significance :
plot_macrosynteny(macrosynteny_df,
                  sp1_label = "B.floridae",
                  sp2_label = "P.echinospica")


## ----eval = FALSE-------------------------------------------------------------
#  my_linkage_groups <- compute_linkage_groups(my_orthologs_with_3_sp)
#  head(my_linkage_groups)
#  

## ----echo = FALSE-------------------------------------------------------------

my_linkage_groups <- compute_linkage_groups(my_orthologs_with_3_sp)
head(my_linkage_groups)


## ----eval = FALSE-------------------------------------------------------------
#  
#  # visualize the loaded data on a oxford grid :
#  my_orthologs_table_reordered <- reorder_macrosynteny(my_orthologs_table)
#  plot_oxford_grid(my_orthologs_table_reordered,
#                   sp1_label = "B.floridae",
#                   sp2_label = "P.echinospica")
#  
#  # compute significance and visualize on a dotplot :
#  macrosynteny_df_reordered <- compute_macrosynteny(my_orthologs_table_reordered)
#  plot_macrosynteny(macrosynteny_df_reordered,
#                    sp1_label = "B.floridae",
#                    sp2_label = "P.echinospica")
#  

## ----echo = FALSE-------------------------------------------------------------
# visualize the loaded data on a oxford grid :
my_orthologs_table_reordered <- reorder_macrosynteny(my_orthologs_table)
plot_oxford_grid(my_orthologs_table_reordered,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica")
# compute significance and visualize on a dotplot :
macrosynteny_df_reordered <- compute_macrosynteny(my_orthologs_table_reordered)
plot_macrosynteny(macrosynteny_df_reordered,
                  sp1_label = "B.floridae",
                  sp2_label = "P.echinospica")


## ----eval = FALSE-------------------------------------------------------------
#  # select only the orthologs falling in the chromosomes of interest and plot:
#  subset_of_orthologs <- subset(my_orthologs_table, sp1.Chr %in% c("BFL13","BFL15","BFL2","BFL3") & sp2.Chr %in% c("PEC2","PEC5","PEC11"))
#  
#  plot_oxford_grid(subset_of_orthologs,
#                   sp1_label = "B.floridae",
#                   sp2_label = "P.echinospica")
#  
#  # reorder :
#  subset_of_orthologs$sp2.Chr <- factor(subset_of_orthologs$sp2.Chr,levels = c("PEC5","PEC11","PEC2"))
#  plot_oxford_grid(subset_of_orthologs,
#                   sp1_label = "B.floridae",
#                   sp2_label = "P.echinospica")
#  
#  # Compute and plot macrosynteny :
#  macrosynteny_of_subset <- compute_macrosynteny(subset_of_orthologs)
#  plot_macrosynteny(macrosynteny_of_subset,
#                   sp1_label = "B.floridae",
#                   sp2_label = "P.echinospica")
#  

## ----echo = FALSE,out.width = c('300px','300px')------------------------------
# select only the orthologs falling in the chromosomes of interest and plot: 
subset_of_orthologs <- subset(my_orthologs_table, sp1.Chr %in% c("BFL13","BFL15","BFL2","BFL3") & sp2.Chr %in% c("PEC2","PEC5","PEC11"))

plot_oxford_grid(subset_of_orthologs,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica")

# reorder :
subset_of_orthologs$sp2.Chr <- factor(subset_of_orthologs$sp2.Chr,levels = c("PEC5","PEC11","PEC2"))
plot_oxford_grid(subset_of_orthologs,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica")



## ----eval = FALSE-------------------------------------------------------------
#  
#  # visualize the loaded data on a oxford grid  with reordering and coloring by cluster :
#  plot_oxford_grid(my_orthologs_table,
#                   sp1_label = "B.floridae",
#                   sp2_label = "P.echinospica",
#                   reorder = TRUE,
#                   color_by = "clust")
#  
#  # redo and color by sp1.Chr instead :
#  plot_oxford_grid(my_orthologs_table,
#                   sp1_label = "B.floridae",
#                   sp2_label = "P.echinospica",
#                   reorder = TRUE,
#                   color_by = "sp1.Chr")
#  

## ----echo = FALSE,out.width = c('300px','300px')------------------------------

# visualize the loaded data on a oxford grid  with reordering and coloring by cluster :
plot_oxford_grid(my_orthologs_table,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica",
                 reorder = TRUE,
                 color_by = "clust")

# redo and color by sp2.Chr instead :
plot_oxford_grid(my_orthologs_table,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica",
                 reorder = TRUE,
                 color_by = "sp2.Chr")


## -----------------------------------------------------------------------------

library(ggplot2)

# legend on right (works also with "top" and "left") :
plot_macrosynteny(macrosynteny_df_reordered) +
  theme(legend.position = "right")


## -----------------------------------------------------------------------------

# Check how many colors are necessary :
print(length(unique(my_orthologs_table_reordered$sp2.Chr)))

# change color_palette using plot_oxford_grid option color_palette :
color_palette_Pechinospica_chromosomes <- c("#A52A2A", "#FFD39B", "#66CDAA", "#8EE5EE", "#7FFF00", "#FFD700", "#FF7F00", "#474747", "#6495ED", "#FF3030", "#0000EE", "#FF1493", "#8A2BE2", "#080808")


plot_oxford_grid(my_orthologs_table_reordered,
                 color_by = "sp2.Chr",
                 color_palette = color_palette_Pechinospica_chromosomes)


## -----------------------------------------------------------------------------

# change the colors in plot_macrosynteny using ggplot2 functions :

plot_macrosynteny(macrosynteny_df_reordered) +
scale_color_manual(values = c("gray45","darkgreen")) +
  theme(legend.position = "right")


## -----------------------------------------------------------------------------

library(dplyr)

# Let's color only the orthologs that were previously selected in the part 3.2 :
my_orthologs_table_modified <- my_orthologs_table_reordered %>%
  mutate(selected = "no") %>%
  mutate(selected = replace(selected,sp1.ID %in% subset_of_orthologs$sp1.ID,"yes"))

plot_oxford_grid(my_orthologs_table_modified,
                 color_by = "selected",
                 color_palette = c("black","firebrick"))

# set the argument shade_non_significant to FALSE to have colors on all the genes of interest :

plot_oxford_grid(my_orthologs_table_modified,
                 color_by = "selected",
                 shade_non_significant = FALSE,
                 color_palette = c("black","firebrick"))



## -----------------------------------------------------------------------------

plot_chord_diagram(my_orthologs_with_3_sp,
                   species_labels = c("B. flo","P. ech", "P. yes"),
                   color_by = "LGs") +
  theme(legend.position = "none")


## ----eval=FALSE---------------------------------------------------------------
#  
#  # Change the chromosome names to keep only numbers
#  levels(my_orthologs_with_3_sp$sp1.Chr) <- stringr::str_replace(levels(my_orthologs_with_3_sp$sp1.Chr),"BFL","")
#  levels(my_orthologs_with_3_sp$sp2.Chr) <- stringr::str_replace(levels(my_orthologs_with_3_sp$sp2.Chr),"PEC","")
#  levels(my_orthologs_with_3_sp$sp3.Chr) <- stringr::str_replace(levels(my_orthologs_with_3_sp$sp3.Chr),"chr","")
#  
#  
#  plot_chord_diagram(my_orthologs_with_3_sp,
#                     species_labels = c("B. flo","P. ech", "P. yes"),
#                     color_by = "LGs") +
#    theme(legend.position = "none")
#  

