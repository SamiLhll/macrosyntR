# plot_macrosynteny
#
# This is a function to plot the result of compute_macrosynteny()
#' @title Plot Macro-synteny
#' @description This is a function to generate the contingency table of an MBH dataframe and apply fischer test to calculate the significant associations.
#'
#' @param macrosynt_df dataframe of contingency table with p-values calculated by the compute_macrosynteny() function
#' @param sp1_label character. The name of the species1 to display on the plot
#' @param sp2_label character. The name of the species2 to put on the plot
#'
#' @seealso [compute_macrosynteny()]
#'
#' @return ggplot2 object
#'
#' @import ggplot2
#' 
#' @examples 
#' # basic usage of plot_macrosynteny : 
#' 
#' orthologs_table <- system.file("extdata","my_orthologs.tab",package="macrosyntR")
#' 
#' my_orthologs <- read.table(orthologs_table,header=TRUE)
#'                                
#' my_macrosynteny <- compute_macrosynteny(my_orthologs)
#' 
#' plot_macrosynteny(my_macrosynteny,
#'                   sp1_label = "B.floridae",
#'                   sp2_label = "P.echinospica")
#' 
#' @export


plot_macrosynteny <- function(macrosynt_df,
                              sp1_label="",
                              sp2_label="") {
  
  sp1.Chr <- sp2.Chr <- orthologs <- significant <- NULL
  
  # Error check : proper format for arguments :
  if (!(is.character(sp1_label) & length(sp1_label) == 1)) { stop("Wrong format for argument 'sp1_label'. Must be a single value of type character")}
  if (!(is.character(sp2_label) & length(sp2_label) == 1)) { stop("Wrong format for argument 'sp2_label'. Must be a single value of type character")}
  # Error check : proper formatting of macrosynt_df
  required_fields <- c("sp1.Chr","sp2.Chr","orthologs","pval","significant","pval_adj")
  for (i in required_fields) {
    if (isFALSE(i %in% colnames(macrosynt_df))) {
      stop("Missing fields in the provided 'macrosynt_df'. All the following columns are required : sp1.Chr, sp2.Chr, orthologs, pval, significant, pval_adj")
    }
  }
  # Error check : macrosynt_df is empty
  if (length(macrosynt_df$sp1.Chr) == 0) {stop("Table provided through the 'macrosynt_df' argument is empty")}
  # Warning check : when number of chromosomes is too high
  if (length(unique(macrosynt_df$sp1.Chr)) >= 300) { 
    warning(paste0("The first species in the 'macrosynt_df' has ",length(unique(macrosynt_df$sp1.Chr))," chromosomes. Computational time can be very long on fragmented genomes"))
  }
  if (length(unique(macrosynt_df$sp2.Chr)) >= 300) { 
    warning(paste0("The second species in the 'macrosynt_df' has ",length(unique(macrosynt_df$sp2.Chr))," chromosomes. Computational time can be very long on fragmented genomes"))
  }
  
  macrosynt_df_to_plot <- macrosynt_df
  
  ### Plot :
  p <- ggplot2::ggplot(macrosynt_df_to_plot,aes(x=sp1.Chr,y=sp2.Chr)) +
    ggplot2::geom_point(aes(size=orthologs,color=significant)) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values=c("#cf6b04","#4baaf2")) +
    ggplot2::scale_y_discrete(limits=rev) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x=sp1_label,y=sp2_label)
  
  
  p <- p + ggplot2::theme(axis.ticks = ggplot2::element_blank(),legend.position = "none",
                          axis.title.x = ggplot2::element_text(family = "sans",size = 12),
                          axis.text.x = ggplot2::element_text(family = "sans",size = 7,angle=90),
                          axis.text.y = ggplot2::element_text(family = "sans",size = 7),
                          axis.title = ggplot2::element_text(family = "sans",size = 12),
                          legend.text = ggplot2::element_text(size = 10),
                          legend.title = ggplot2::element_text(size = 12),
                          plot.title = ggplot2::element_text(family = "sans",size=12,hjust = 0.5),
                          plot.background = ggplot2::element_rect(fill = "white",colour = "white"))
  
  
  return(p)
  
}