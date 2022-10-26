# plot_macrosynt
#
# This is a function to plot the result of calculate_macrosynt()
#' @title Plot Macro-synteny
#' @description This is a function to generate the contingency table of an MBH dataframe and apply fischer test to calculate the significant associations.
#'
#' @param macrosynt_df dataframe of contingency table with p-values calculated by the calculate_contingency_table function
#' @param sp1_label string. The name of the species1 to display on the plot
#' @param sp2_label string. The name of the species2 to put on the plot
#'
#' @return ggplot2 object
#'
#' @import dplyr
#' @import ggplot2
#' @export


plot_macrosynt <- function(macrosynt_df,
                           sp1_label="",
                           sp2_label="") {
  
  sp1_chr <- sp2_chr <- orthologs <- significant <- NULL
  
  macrosynt_df_to_plot <- macrosynt_df
  
  ### Plot :
  p <- ggplot2::ggplot(macrosynt_df_to_plot,aes(x=sp1_chr,y=sp2_chr)) +
    ggplot2::geom_point(aes(size=orthologs,color=significant)) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values=c("#cf6b04","#4baaf2")) +
    ggplot2::scale_y_discrete(limits=rev) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x=sp1_label,y=sp2_label)
  
  
  p <- p + ggplot2::theme(axis.ticks = element_blank(),legend.position = "none",
                          axis.title.x = element_text(family = "sans",size = 12),
                          axis.text.x = element_text(family = "sans",size = 7,angle=90),
                          axis.text.y = element_text(family = "sans",size = 10),
                          axis.title = element_text(family = "sans",size = 12),
                          legend.text = element_text(size = 10),
                          legend.title = element_text(size = 12),
                          plot.title = element_text(family = "sans",size=12,hjust = 0.5),
                          plot.background = element_rect(fill = "white",colour = "white"))
  
  
  return(p)
  
}