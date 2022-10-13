# plot_fischer_test
#
# This is a function to plot the result of the calculate_contingency_table function
#' @title Plot fischer test
#' @description This is a function to generate the contingency table of an MBH dataframe and apply fischer test to calculate the significant associations.
#'
#' @param contingency_table dataframe of contingency table with p-values calculated by the calculate_contingency_table function
#' @param sp1_name string. The name of the species1 to display on the plot
#' @param sp2_name string. The name of the species2 to put on the plot
#' @param plot_legend boolean set to TRUE to display the legend (default = FALSE)
#' @param sp1_chr_order ordered list of chromosome names to use for plotting
#' @param sp2_chr_order ordered list of chromosome names to use for plotting
#'
#' @return A ggplot2 object
#'
#' @import dplyr
#' @import ggplot2
#' @export


plot_fischer_test <- function(contingency_table,sp1_name="",sp2_name="",plot_legend = FALSE,sp1_chr_order = NULL,sp2_chr_order = NULL) {
  contingency_table_to_plot <- contingency_table
  ### Reorder chromosomes if necessary :
  if (! is.null(sp1_chr_order)) { contingency_table_to_plot$sp1_chr <- factor(as.character(contingency_table_to_plot$sp1_chr),levels = sp1_chr_order)}
  if (! is.null(sp2_chr_order)) { contingency_table_to_plot$sp2_chr <- factor(as.character(contingency_table_to_plot$sp2_chr),levels = sp2_chr_order)}
  ### Plot :
  p <- ggplot(contingency_table_to_plot,aes(x=sp1_chr,y=sp2_chr)) +
    geom_point(aes(size=a,color=significant)) +
    theme_bw() +
    scale_color_manual(values=c("#cf6b04","#4baaf2")) +
    scale_y_discrete(limits=rev) +
    scale_x_discrete(position = "top") +
    labs(x = "",y = sp2_name,color="significant < 10-3 ",size="# of orthologs :",title = sp1_name)
  
  if (plot_legend) {
    p <- p + theme(axis.ticks = element_blank(),legend.position = "top",
                   axis.title.x = element_text(family = "sans",size = 12),
                   axis.text.x = element_text(family = "sans",size = 7),
                   axis.text.y = element_text(family = "sans",size = 10),
                   axis.title = element_text(family = "sans",size = 12),
                   legend.text = element_text(size = 10),
                   legend.title = element_text(size = 12),
                   plot.title = element_text(family = "sans",size=12,hjust = 0.5),
                   plot.background = element_rect(fill = "white",colour = "white"))
  }
  else {
    p <- p + theme(axis.ticks = element_blank(),legend.position = "none",
                   axis.title.x = element_text(family = "sans",size = 12),
                   axis.text.x = element_text(family = "sans",size = 7),
                   axis.text.y = element_text(family = "sans",size = 10),
                   axis.title = element_text(family = "sans",size = 12),
                   legend.text = element_text(size = 10),
                   legend.title = element_text(size = 12),
                   plot.title = element_text(family = "sans",size=12,hjust = 0.5),
                   plot.background = element_rect(fill = "white",colour = "white"))
  }
    
  
  return(p)
  
}