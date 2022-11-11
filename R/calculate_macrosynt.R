# calculate_macrosynt
#
# This is a function to generate the contingency table of an MBH dataframe and apply fischer test to calculate the significant associations.
#' @title Calculate contingency table
#' @description This is a function to generate the contingency table of an MBH dataframe and 
#'   apply fischer test to calculate the significant associations. It outputs a datafram shaped as following :
#'   sp1.Chr,sp2.Chr,orthologs,pval
#'
#' @param mbh_df dataframe. mutual best hits with genomic coordinates loaded with load_mbh_df()
#' @param pvalue_threshold numeric. threshold for significancy. (default equals 0.001)
#' 
#' @return A dataframe object
#'
#' @importFrom stats fisher.test
#' @importFrom stats p.adjust
#' @import tidyr
#' @importFrom dplyr rename mutate select
#' 
#' @export


calculate_macrosynt <- function(mbh_df,pvalue_threshold = 0.001) {
  
  ### construct the contingency table :
  final_i <- final_j <- final_a <- final_b <- final_c <- final_d <- NULL
  sp1.Chr <- sp2.Chr <- value <- variable <- pvalues_value <- pvalues_adj <- odds <- a <- significant <-NULL
  
  for (i in unique(mbh_df$sp1.Chr)) {
    subsetted_mbh_df_in <- subset(mbh_df,sp1.Chr == i)
    subsetted_mbh_df_out <- subset(mbh_df,sp1.Chr != i)
    compared_specie_scaffs <- unique(subsetted_mbh_df_in$sp2.Chr)
    for (j in compared_specie_scaffs) {
      cell_a <- length(subset(subsetted_mbh_df_in, sp2.Chr == j)$sp1.Chr)
      final_a <- c(final_a,cell_a)
      cell_b <- length(subset(subsetted_mbh_df_out, sp2.Chr == j)$sp1.Chr)
      final_b <- c(final_b,cell_b)
      cell_c <- length(subset(subsetted_mbh_df_in, sp2.Chr != j)$sp1.Chr)
      final_c <- c(final_c,cell_c)
      cell_d <- length(subset(subsetted_mbh_df_out, sp2.Chr != j)$sp1.Chr)
      final_d <- c(final_d,cell_d)
      final_j <- c(final_j,j)
    }
    final_i <- c(final_i,rep(i,times=length(compared_specie_scaffs)))
  }
  
  contingency_table <- data.frame(sp1 = final_i,
                                  sp2 = final_j,
                                  a = final_a,b = final_b,
                                  c = final_c,d = final_d)
  ### DONE
  
  rownames(contingency_table) <- paste(contingency_table$sp1,
                                       contingency_table$sp2,sep="-")
  contingency_table <- subset(contingency_table,select=c("a","b","c","d"))
  
  ### Calculate pvalues and odds :
  p_values <- apply(contingency_table, 1,
                    function(x) {
                      tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                      stats::fisher.test(tbl, alternative="greater")$p.value
                    })
  odds_values <- apply(contingency_table, 1,
                       function(x) {
                         tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                         stats::fisher.test(tbl, alternative="greater")$estimate[[1]]
                       })
  
  ### merge pvalues and odds with the contingency table :
  contingency_table$variable <- rownames(contingency_table)
  p_values.melted <- reshape2::melt(p_values) %>%
    dplyr::rename(pvalues_value = value)
  p_values.melted$variable <- rownames(p_values.melted)
  
  odds_values.melted <- reshape2::melt(odds_values) %>%
    dplyr::rename(odds_values = value)
  odds_values.melted$variable <- rownames(odds_values.melted)
  
  contingency_table <- merge(contingency_table,p_values.melted) %>%
    dplyr::mutate(pvalues_adj = stats::p.adjust(pvalues_value)) %>%
    dplyr::mutate(significant = "yes",
                  significant = replace(significant,pvalues_adj > pvalue_threshold,"no"))
  
  contingency_table <- merge(contingency_table,odds_values.melted) %>%
    dplyr::mutate(odds = "> 1",
                  odds = replace(odds,odds_values < 1,"< 1")) %>%
    tidyr::separate(variable,c("sp1.Chr","sp2.Chr"),sep="-")
  
  ### reshape before returning :
  macrosynt_df <- contingency_table %>%
    dplyr::select(sp1.Chr,sp2.Chr,a,pvalues_value,significant,pvalues_adj) %>%
    dplyr::rename(orthologs = a,pval = pvalues_value)
  # copy the levels of mbh_df to keep the same ordering when plotting :
  macrosynt_df$sp1.Chr <- factor(macrosynt_df$sp1.Chr,levels = levels(mbh_df$sp1.Chr))
  macrosynt_df$sp2.Chr <- factor(macrosynt_df$sp2.Chr,levels = levels(mbh_df$sp2.Chr))
  
  return(macrosynt_df)
  
}