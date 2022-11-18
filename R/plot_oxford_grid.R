# plot_oxford_grid
#
# This is a function to plot the oxford grided plot to compare the macro synteny of two species
#' @title plot the Macro-synteny as an oxford grid.
#' @description This is a function to plot the oxford grided plot to compare the macro synteny of two species. Its input will have been loaded using load_orthologs()
#'
#' @param orthologs_df dataframe. orthologs with genomic coordinates loaded by the load_orthologs()
#' @param sp1_label string. name of 1st species to display on the plot
#' @param sp2_label string. name of 2nd species to display on the plot
#' @param dot_size float. (default = 0.5)
#' @param dot_alpha flot. (default = 0.4)
#' @param auto_order_clusters logical. (default = FALSE) tells whether to reorder the chromosomes in clusters as implemented in reorder_macrosynteny()
#' @param color_clusters logical. (default = FALSE) color the orthologous pairs in a significant linkage group
#' @param pvalue_threshold float. (default = 0.001) 
#' @param clusters_color_palette vector. (default = NULL) list of colors (as string under double quote) for the clusters. The amount of colors must match the amount of clusters.
#' 
#' @seealso [load_orthologs()]
#' @seealso [reorder_macrosynteny()]
#' 
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @import ggthemes
#' @export


plot_oxford_grid <- function(orthologs_df,
                                 sp1_label = "",
                                 sp2_label = "",
                                 dot_size = 0.5,
                                 dot_alpha = 0.4,
                                 auto_order_clusters = FALSE,
                                 color_clusters = FALSE,
                                 pvalue_threshold= 0.001,
                                 clusters_color_palette = NULL) {
  
  sp1.Index <- sp2.Index <- sp2.Chr <-significant <- clust <- NULL
  
  orthologs_df_to_plot <- orthologs_df
  ### reorder df first :
  if (auto_order_clusters) {
    # calculate clusters and reordered synteny
    orthologs_reordered <- reorder_macrosynteny(orthologs_df)
    orthologs_df_to_plot <- orthologs_reordered
  }
  # separate dots not in clusters, and dots in clusters
  if (color_clusters) {
    ### [Exception here] Check that orthologs_df_to_plot has the clust column
    if (!("clust" %in% colnames(orthologs_df_to_plot))) { 
      stop("Asked to color the clusters but the clust column couldn't be found in the data. Make sure to set auto_order_cluster = TRUE or use reorder_macrosynteny()")
    }
    ###
    # convert to character for discrete values coloring :
    temp_macrosynt <- compute_macrosynteny(orthologs_df)
    temp_orthologs_and_macrosynt <- merge(orthologs_df_to_plot,temp_macrosynt)
    final_df_with_groups <- subset(temp_orthologs_and_macrosynt,significant == "yes")
    non_linkage_df <- subset(temp_orthologs_and_macrosynt, significant == "no")
    
    # initialize the plot with colors :
    p <- ggplot2::ggplot(final_df_with_groups,ggplot2::aes(x=sp1.Index,y=sp2.Index,color = clust)) + 
      ggplot2::geom_point(size=dot_size,alpha=dot_alpha) +
      ggplot2::geom_point(data = non_linkage_df,ggplot2::aes(x=sp1.Index,y=sp2.Index),
                          color="gray",
                          size = dot_size,
                          alpha = dot_alpha)
    if (! is.null(clusters_color_palette)) {
      #### [Exception here] : check that length(clusters_color_palette) == length(levels(final_df_with_groups$clust))
      p <- p +
        ggplot2::scale_color_manual(values = clusters_color_palette)
    }
  }
  else {
    p <- ggplot2::ggplot(orthologs_df_to_plot,ggplot2::aes(x=sp1.Index,y=sp2.Index)) + ggplot2::geom_jitter(size=dot_size,alpha=dot_alpha)
  }
  ### build the plot :
  p <- p + ggthemes::theme_tufte() +
    ggplot2::facet_grid(sp2.Chr ~ sp1.Chr,scales = "free",space = "free") +
    ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   strip.text.x=ggplot2::element_text(family="sans",angle=90,size = 7),
                   strip.text.y=ggplot2::element_text(family="sans",angle=0,size = 7),
                   axis.title.x =ggplot2::element_text(family="sans",size=10),
                   axis.title.y =ggplot2::element_text(family="sans",size=12),
                   plot.title =ggplot2::element_text(family="sans",size=12,hjust = 0.5),
                   strip.background=ggplot2::element_rect(colour="white",fill="white"),
                   plot.background = ggplot2::element_rect(fill = "white",colour = "white"),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.spacing = ggplot2::unit(0.01, "lines"),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA, color = "gray",size=0.1),
                   panel.background = ggplot2::element_rect(fill = "white", colour = "white")) +
    ggplot2::labs(y=sp2_label, x= paste0("(",length(orthologs_df_to_plot$sp1.ID)," orthologs)"),title = sp1_label)
  
  return(p)
  
}