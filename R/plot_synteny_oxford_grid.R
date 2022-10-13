# plot_synteny_oxford_grid
#
# This is a function to plot the oxford grided plot to compare the macro synteny of two species
#' @title plot the synteny oxford grid.
#' @description This is a function to plot the oxford grided plot to compare the macro synteny of two species. Its input will have been loaded using the load_MBH_table function, and the output is a ggplot2 object.
#'
#' @param MBH_table dataframe of mutual best hits with genomic coordinates loaded by the load_MBH_table function
#' @param sp1_name string. name of 1st species to display on the plot
#' @param sp2_name string. name of 2nd species to display on the plot
#' @param sp1_keep_chr_names list of chr names to keep on species1.(default is to keep all)
#' @param sp2_keep_chr_names list of chr names to keep on species 2. (default is to keep all)
#' @param sp1_chr_order ordered list of chromosome names to use for plotting
#' @param sp2_chr_order ordered list of chromosome names to use for plotting
#'
#' @return A ggplot2 object
#'
#' @import dplyr
#' @import ggplot2
#' @import ggthemes
#' @export


plot_synteny_oxford_grid <- function(MBH_table,colors = FALSE,sp1_name = "",sp2_name = "",sp1_keep_chr_names = NULL,sp2_keep_chr_names = NULL,sp1_chr_order = NULL,sp2_chr_order = NULL,dot_size = 0.2) {

  MBH_table_to_plot <- MBH_table
  ### Filter the chromosomes if necessary :
  if (! is.null(sp1_keep_chr_names)) { MBH_table_to_plot <- subset(MBH_table_to_plot,sp1_chr %in% sp1_keep_chr_names)}
  if (! is.null(sp2_keep_chr_names)) { MBH_table_to_plot <- subset(MBH_table_to_plot,sp2_chr %in% sp2_keep_chr_names)}
  ### Reorder chromosomes if necessary :
  if (! is.null(sp1_chr_order)) { MBH_table_to_plot$sp1_chr <- factor(as.character(MBH_table_to_plot$sp1_chr),levels = sp1_chr_order)}
  if (! is.null(sp2_chr_order)) { MBH_table_to_plot$sp2_chr <- factor(as.character(MBH_table_to_plot$sp2_chr),levels = sp2_chr_order)}
  ### deal with colors first :
  if (colors) {p <- ggplot(MBH_table_to_plot,aes(x=sp1_index,y=sp2_index,color = sp2_chr)) + geom_jitter(size=dot_size,alpha=0.3) }
  else {p <- ggplot(MBH_table_to_plot,aes(x=sp1_index,y=sp2_index)) + geom_jitter(size=dot_size,alpha=0.3)}
  ### build the plot :
  p <- p + theme_tufte() +
    facet_grid(sp2_chr ~ sp1_chr,scales = "free",space = "free") +
    theme(axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          strip.text.x=element_text(family="sans",angle=0,size = 5),
          strip.text.y=element_text(family="sans",angle=0,size = 6),
          axis.title.x =element_text(family="sans",size=10),
          axis.title.y =element_text(family="sans",size=12),
          plot.title = element_text(family="sans",size=12,hjust = 0.5),
          strip.background=element_rect(colour="white",fill="white"),
          plot.background = element_rect(fill = "white",colour = "white"),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.01, "lines"),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, color = "gray",size=0.1),
          panel.background = element_rect(fill = "white", colour = "white")) +
    labs(y=sp2_name, x= paste0("(",length(MBH_table_to_plot$sp1_pep)," orthologs)"),title = sp1_name)

  return(p)

}
