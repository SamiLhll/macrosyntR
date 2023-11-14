# plot_oxford_grid
#
# This is a function to plot the oxford grided plot to compare the macro synteny of two species
#' @title plot the Macro-synteny as an oxford grid.
#' @description This is a function to plot the oxford grided plot to compare the macro synteny of two species. It requires as input an orthologs_df loaded by load_orthologs()
#'
#' @param orthologs_df dataframe. orthologs with genomic coordinates loaded by the load_orthologs()
#' @param sp1_label character. name of 1st species to display on the plot
#' @param sp2_label character. name of 2nd species to display on the plot
#' @param dot_size numeric. (default = 0.5)
#' @param dot_alpha numeric. (default = 0.4)
#' @param reorder logical. (default = FALSE) tells whether to reorder the chromosomes in clusters as implemented in reorder_macrosynteny()
#' @param keep_only_significant logical. (default = FALSE)
#' @param color_by string/variable name. (default = NULL) column of the orthologs_df to use to color the dots.
#' @param pvalue_threshold numeric. (default = 0.001) 
#' @param color_palette vector. (default = NULL) list of colors (as string under double quote) for the clusters. The amount of colors must match the amount of clusters.
#' @param shade_non_significant logical. (default = TRUE) When TRUE the orthologs located on non-significant linkage groups are displayed in "grey"
#' @param reverse_species logical. (default = FALSE) When TRUE the x and y axis of the plot are reversed. sp1 is displayed on the y axis and sp2 is displayed on the x axis.
#' @param keep_sp1_raw_order logical.(default equals FALSE) tells if the reordering should be constrained on the species1 and change just the order of the species2
#' 
#' @seealso [load_orthologs()]
#' @seealso [reorder_macrosynteny()]
#' 
#' @return A ggplot2 object
#'
#' @import ggplot2
#' 
#' @examples 
#' # basic usage of plot_oxford_grid : 
#' 
#' orthologs_table <- system.file("extdata","my_orthologs.tab",package="macrosyntR")
#' 
#' my_orthologs <- read.table(orthologs_table,header=TRUE)
#'
#' plot_oxford_grid(my_orthologs,
#'                  sp1_label = "B. floridae",
#'                  sp2_label = "P. echinospica")
#' 
#' # plot a reordered Oxford Grid and color by cluster :
#' \donttest{
#' plot_oxford_grid(my_orthologs,
#'                  sp1_label = "B. floridae",
#'                  sp2_label = "P. echinospica",
#'                  reorder = TRUE,
#'                  color_by = "clust")
#'  }
#' @export


plot_oxford_grid <- function(orthologs_df,
                             sp1_label = "",
                             sp2_label = "",
                             dot_size = 0.5,
                             dot_alpha = 0.4,
                             reorder = FALSE,
                             keep_only_significant = FALSE,
                             color_by = NULL,
                             pvalue_threshold= 0.001,
                             color_palette = NULL,
                             shade_non_significant = TRUE,
                             reverse_species = FALSE,
                             keep_sp1_raw_order = FALSE) {
  
  sp1.Index <- sp2.Index <- sp2.Chr <-significant <- clust <- NULL
  
  orthologs_df_to_plot <- orthologs_df
  # Error check : proper format for arguments :
  if (!(is.character(sp1_label) & length(sp1_label) == 1)) { stop("Wrong format for argument 'sp1_label'. Must be a single value of type character")}
  if (!(is.character(sp2_label) & length(sp2_label) == 1)) { stop("Wrong format for argument 'sp2_label'. Must be a single value of type character")}
  if (!(is.numeric(dot_size) & length(dot_size) == 1)) { stop("Wrong format for argument 'dot_size'. Must be a single value of type numeric")}
  if (!(is.numeric(dot_alpha) & length(dot_alpha) == 1)) { stop("Wrong format for argument 'dot_alpha'. Must be a single value of type numeric")}
  if (!(is.logical(reorder) & length(reorder) == 1)) { stop("Wrong format for argument 'reorder'. Must be a single value of type logical")}
  if (!(is.logical(keep_only_significant) & length(keep_only_significant) == 1)) { stop("Wrong format for argument 'keep_only_significant'. Must be of type logical")}
  if (!(is.logical(shade_non_significant) & length(shade_non_significant) == 1)) { stop("Wrong format for argument 'shade_non_significant'. Must be of type logical")}
  if (!(is.logical(reverse_species) & length(reverse_species) == 1)) { stop("Wrong format for argument 'reverse_species'. Must be of type logical")}
  if (!(is.logical(keep_sp1_raw_order) & length(keep_sp1_raw_order) == 1)) { stop("Wrong format for argument 'keep_sp1_raw_order'. Must be of type logical")}
  
  if(!(is.null(color_by))) {
    if (reorder) {
      if (!((color_by %in% colnames(orthologs_df) | color_by == "clust") & length(color_by) == 1 & is.character(color_by))) {
        stop("Wrong format for 'color_by' argument. Must be a single value of type character and correspond to the name of one column of the orthologs_df")
      }
    }
    else if (!(color_by %in% colnames(orthologs_df) & length(color_by) == 1 & is.character(color_by))) { 
      stop("Wrong format for 'color_by' argument. Must be a single value of type character and correspond to the name of one column of the orthologs_df")
      }
  }
  if (!(is.numeric(pvalue_threshold) & length(pvalue_threshold) == 1)) {stop("Wrong format for 'pvalue_threshold' argument. Must be a single value of type numeric")}
  if (!(is.null(color_by))) {
    if (color_by != "clust") {
      if (!(is.null(color_palette)) & (length(color_palette) < length(unique(orthologs_df[[color_by]])))) {stop(paste("Wrong format in argument 'color_palette'. Must be a list of colors with as much values as the amount of unique elements in the column specified in the argument 'color_by'. color palette has ",
                                                                                                               length(color_palette)," elements, while the color_by column is ",length(unique(orthologs_df[[color_by]]))))}
    }
  }
  
  # Error check : format of orthologs_df 
  required_fields <- c("sp1.ID","sp1.Index","sp1.Chr","sp2.ID","sp2.Index","sp2.Chr")
  for (i in required_fields) {
    if (isFALSE(i %in% colnames(orthologs_df))) {
      stop("Missing fields in the provided 'orthologs_df'. All the following columns are required : sp1.ID,sp1.Index,sp1.Chr,sp2.ID,sp2.Index,sp2.Chr")
    }
  }
  # Error check : orthologs_df is empty
  if (length(orthologs_df$sp1.Chr) == 0) {stop("Table provided through the 'orthologs_df' argument is empty")}
  # Warning check : when number of chromosomes is too high
  if (length(unique(orthologs_df$sp1.Chr)) >= 300) { 
    warning(paste0("The first species in the 'orthologs_df' has ",length(unique(orthologs_df$sp1.Chr))," chromosomes. Computational time can be very long on fragmented genomes"))
  }
  if (length(unique(orthologs_df$sp2.Chr)) >= 300) { 
    warning(paste0("The second species in the 'orthologs_df' has ",length(unique(orthologs_df$sp2.Chr))," chromosomes. Computational time can be very long on fragmented genomes"))
  }
  ### reverse species if needed :
  if (reverse_species) {
    orthologs_df <- reverse_species_order(orthologs_df)
  }
  
  ### reorder df first :
  if (reorder) {
    # calculate clusters and reordered synteny
    orthologs_reordered <- reorder_macrosynteny(orthologs_df,pvalue_threshold = pvalue_threshold,
                                                keep_only_significant = keep_only_significant,
                                                keep_sp1_raw_order = keep_sp1_raw_order)
    orthologs_df_to_plot <- orthologs_reordered
  }
  # separate dots not in clusters, and dots in clusters
  if (! is.null(color_by)) {
    # ### [Exception here] Check that orthologs_df_to_plot has the clust column
    # if (!(color_ %in% colnames(orthologs_df_to_plot))) { 
    #   stop("Asked to color the clusters but the clust column couldn't be found in the data. Make sure to set reorder = TRUE or use reorder_macrosynteny()")
    # }
    ###
    # convert to character for discrete values coloring :
    temp_macrosynt <- compute_macrosynteny(orthologs_df)
    temp_orthologs_and_macrosynt <- merge(orthologs_df_to_plot,temp_macrosynt)
    final_df_with_groups <- subset(temp_orthologs_and_macrosynt,significant == "yes")
    non_linkage_df <- subset(temp_orthologs_and_macrosynt, significant == "no")
    
    # initialize the plot with colors :
    color_by_sym <- ggplot2::ensym(color_by)
    p <- ggplot2::ggplot(final_df_with_groups,ggplot2::aes(x=sp1.Index,y=sp2.Index,color = !!color_by_sym)) + 
      ggplot2::geom_point(size=dot_size,alpha=dot_alpha)
    # Check shade_non_significant:
    if (shade_non_significant) {
      p <- p + ggplot2::geom_point(data = non_linkage_df,ggplot2::aes(x=sp1.Index,y=sp2.Index),
                                   color="gray",
                                   size = dot_size,
                                   alpha = dot_alpha)
    }
    else { p <- p + ggplot2::geom_point(data = non_linkage_df,ggplot2::aes(x=sp1.Index,y=sp2.Index,color = !!color_by_sym),
                                        size = dot_size,
                                        alpha = dot_alpha)
    }
      
    if (! is.null(color_palette)) {
      #### [Exception here] : check that length(clusters_color_palette) == length(levels(final_df_with_groups$clust))
      p <- p +
        ggplot2::scale_color_manual(values = color_palette)
    }
  }
  else {
    p <- ggplot2::ggplot(orthologs_df_to_plot,ggplot2::aes(x=sp1.Index,y=sp2.Index)) + ggplot2::geom_jitter(size=dot_size,alpha=dot_alpha)
  }
  ### build the plot :
  p <- p + ggplot2::theme_void() +
    ggplot2::facet_grid(sp2.Chr ~ sp1.Chr,scales = "free",space = "free") +
    ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   strip.text.x=ggplot2::element_text(family="sans",angle=90,size = 7),
                   strip.text.y=ggplot2::element_text(family="sans",angle=0,size = 7),
                   axis.title.x =ggplot2::element_text(family="sans",size=10),
                   axis.title.y =ggplot2::element_text(family="sans",size=12),
                   plot.subtitle =ggplot2::element_text(family="sans",size=12,hjust = 0.5),
                   strip.background=ggplot2::element_rect(colour="white",fill="white"),
                   plot.background = ggplot2::element_rect(fill = "white",colour = "white"),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.spacing = ggplot2::unit(0.01, "lines"),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA, color = "gray",size=0.1),
                   panel.background = ggplot2::element_rect(fill = "white", colour = "white")) +
    ggplot2::labs(y=sp2_label, x= paste0("(",length(orthologs_df_to_plot$sp1.ID)," orthologs)"),subtitle = sp1_label)
  
  return(p)
  
}