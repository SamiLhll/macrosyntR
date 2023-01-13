## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(macrosyntR)

## -----------------------------------------------------------------------------
my_orthologs_table <- load_orthologs(orthologs_table = system.file("extdata","Bflo_vs_Pech.tab",package="macrosyntR"),
                                     sp1_bed = system.file("extdata","Bflo.protein_products.bed",package="macrosyntR"),
                                     sp2_bed = system.file("extdata","Pech.protein_products.bed",package="macrosyntR"))

head(my_orthologs_table)


## -----------------------------------------------------------------------------

# compute significance :
macrosynteny_df <- compute_macrosynteny(my_orthologs_table)
head(macrosynteny_df)


## ----eval = FALSE-------------------------------------------------------------
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


