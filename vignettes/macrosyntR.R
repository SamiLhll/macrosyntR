## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(macrosyntR)

## -----------------------------------------------------------------------------
my_orthologs_table <- load_mbh_df(mbh_table = "../inst/extdata/Bflo_vs_Pech.tab",
                                  sp1_bed = "../inst/extdata/Bflo.protein_products.bed",
                                  sp2_bed = "../inst/extdata/Pech.protein_products.bed")

head(my_orthologs_table)


## -----------------------------------------------------------------------------
# visualize the loaded data on a oxford grid :
plot_oxford_grid(my_orthologs_table,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica")
# compute significance and visualize on a dotplot :
macrosynteny_df <- calculate_macrosynt(my_orthologs_table)
head(macrosynteny_df)
plot_macrosynt(macrosynteny_df,
               sp1_label = "B.floridae",
               sp2_label = "P.echinospica")


## -----------------------------------------------------------------------------
# visualize the loaded data on a oxford grid :
my_orthologs_table_reordered <- reorder_synteny(my_orthologs_table)
plot_oxford_grid(my_orthologs_table_reordered,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica")
# compute significance and visualize on a dotplot :
macrosynteny_df_reordered <- calculate_macrosynt(my_orthologs_table_reordered)
plot_macrosynt(macrosynteny_df_reordered,
               sp1_label = "B.floridae",
               sp2_label = "P.echinospica")


## -----------------------------------------------------------------------------
# visualize the loaded data on a oxford grid :
plot_oxford_grid(my_orthologs_table,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica",
                 auto_order_clusters = TRUE,
                 color_clusters = TRUE)

# Change the vector of colors, 1 color per cluster, we see on the plot above that we have 10 clusters :
plot_oxford_grid(my_orthologs_table,
                 sp1_label = "B.floridae",
                 sp2_label = "P.echinospica",
                 auto_order_clusters = TRUE,
                 color_clusters = TRUE,
                 clusters_color_palette = c("#FFA500","#FFA500", "#1E90FF","#FFA500", "#1E90FF", "#1E90FF", "#1E90FF", "#228B22", "#228B22", "#228B22"))



