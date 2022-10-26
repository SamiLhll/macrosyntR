# plot_oxford_grid
#
# This is a function to plot the oxford grided plot to compare the macro synteny of two species
#' @title plot the Macro-synteny as an oxford grid.
#' @description This is a function to plot the oxford grided plot to compare the macro synteny of two species. Its input will have been loaded using load_mbh_df()
#'
#' @param mbh_df dataframe. mutual best hits with genomic coordinates loaded by the load_mbh_df()
#' @param sp1_label string. name of 1st species to display on the plot
#' @param sp2_label string. name of 2nd species to display on the plot
#' @param dot_size float. default = 0.2
#' @param colors logical. default = FALSE 
#' 
#' @seealso [load_mbh_df()]
#' 
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @import ggthemes
#' @export


plot_oxford_grid <- function(mbh_df,
                             sp1_label = "",
                             sp2_label = "",
                             dot_size = 0.2,
                             colors = FALSE) {
  
  sp1_index <- sp2_index <- sp2_chr <- NULL

  mbh_df_to_plot <- mbh_df
  ### deal with colors first :
  if (colors) {
    p <- ggplot2::ggplot(mbh_df_to_plot,aes(x=sp1_index,y=sp2_index,color = sp2_chr)) + 
      ggplot2::geom_jitter(size=dot_size,alpha=0.3) 
    }
  else {
    p <- ggplot2::ggplot(mbh_df_to_plot,aes(x=sp1_index,y=sp2_index)) + ggplot2::geom_jitter(size=dot_size,alpha=0.4)
    }
  ### build the plot :
  p <- p + ggthemes::theme_tufte() +
    ggplot2::facet_grid(sp2_chr ~ sp1_chr,scales = "free",space = "free") +
    ggplot2::theme(axis.text.y=element_blank(),
          axis.text.x=ggplot2::element_blank(),
          strip.text.x=ggplot2::element_text(family="sans",angle=90,size = 5),
          strip.text.y=ggplot2::element_text(family="sans",angle=0,size = 6),
          axis.title.x =ggplot2::element_text(family="sans",size=10),
          axis.title.y =ggplot2::element_text(family="sans",size=12),
          plot.title =ggplot2::element_text(family="sans",size=12,hjust = 0.5),
          strip.background=ggplot2::element_rect(colour="white",fill="white"),
          plot.background = ggplot2::element_rect(fill = "white",colour = "white"),
          axis.ticks = ggplot2::element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.01, "lines"),
          panel.grid = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(fill = NA, color = "gray",size=0.1),
          panel.background = ggplot2::element_rect(fill = "white", colour = "white")) +
    ggplot2::labs(y=sp2_label, x= paste0("(",length(mbh_df_to_plot$sp1_pep)," orthologs)"),title = sp1_label)

  return(p)

}
