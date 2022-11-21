## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(macrosyntR)

## -----------------------------------------------------------------------------
my_orthologs_table <- load_orthologs(orthologs_table = "../inst/extdata/Bflo_vs_Pech.tab",
                                     sp1_bed = "../inst/extdata/Bflo.protein_products.bed",
                                     sp2_bed = "../inst/extdata/Pech.protein_products.bed")

head(my_orthologs_table)


## -----------------------------------------------------------------------------
# visualize the loaded data on a oxford grid :
plot_oxford_grid(my_orthologs_table,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica")
# compute significance and visualize on a dotplot :
macrosynteny_df <- compute_macrosynteny(my_orthologs_table)
head(macrosynteny_df)
plot_macrosynteny(macrosynteny_df,
                  sp1_label = "B.floridae",
                  sp2_label = "P.echinospica")


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
# visualize the loaded data on a oxford grid :
plot_oxford_grid(my_orthologs_table,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica",
                 reorder = TRUE,
                 color_clusters = TRUE)


