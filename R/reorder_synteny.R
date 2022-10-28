# reorder_syneny
#
# This is a function to reorder an mbh_df, that was generated with load_mbh_df(). It allows for automatic reordering, or manual filtering and reordering.
#' @title Reorder the mbh_df before plotting
#' @description his is a function to reorder an mbh_df, that was generated with load_mbh_df(). It allows for automatic reordering, or manual filtering and reordering.
#'
#' @param mbh_df dataframe. mutual best hits with genomic coordinates loaded with load_mbh_df()
#' @param pvalue_threshold numeric. threshold for significancy. (default equals 0.001)
#' @return A dataframe object
#'
#' @import dplyr
#' 
#' @export


reorder_synteny <- function(mbh_df,
                            pvalue_threshold = 0.001) {
  
  sp1.Chr <- sp2.Chr <- significant <- orthologs <- NULL
  
  
  # 1 - sort the sp1 chr by decreasing order of macrosynteny group sizes :
  contingency_table <- calculate_macrosynt(mbh_df,pvalue_threshold)
  significant_entries <- subset(contingency_table, significant == "yes") %>%
    dplyr::arrange(desc(orthologs))
  
  canva_levels_sp1 <- unique(significant_entries$sp1.Chr)
  mbh_df_levels_reordered <- subset(mbh_df,sp1.Chr %in% canva_levels_sp1)
  
  # 2 - sort the sp2 chr :
  levels_sp1 <- NULL
  levels_sp2 <- NULL
  for (i in canva_levels_sp1) {
    # get all sp2 chromosomes significant with sp1 :
    temp_significant_sp2_1st_round <- subset(significant_entries,sp1.Chr == i)
    # get all sp1 chromosomes significant with these sp2 :
    temp_significant_sp1_1st_round <- subset(significant_entries,sp2.Chr %in% as.character(temp_significant_sp2_1st_round$sp2.Chr))
    # additional filtering of sp2 chromosomes that are already placed in the table :
    temp_significant_sp2_2nd_round <- temp_significant_sp2_1st_round[ !(temp_significant_sp2_1st_round$sp2.Chr %in% levels_sp2),]
    temp_significant_sp1_2nd_round <- temp_significant_sp1_1st_round[ !(temp_significant_sp1_1st_round$sp1.Chr %in% levels_sp1),]
    levels_sp2 <- c(levels_sp2,unique(as.character(temp_significant_sp2_2nd_round$sp2.Chr)))
    levels_sp1 <- c(levels_sp1,unique(as.character(temp_significant_sp1_2nd_round$sp1.Chr)))
  }
  
  mbh_df_levels_reordered <- subset(mbh_df_levels_reordered, sp2.Chr %in% levels_sp2)
  mbh_df_levels_reordered$sp1.Chr <- factor(mbh_df_levels_reordered$sp1.Chr,levels = levels_sp1)
  mbh_df_levels_reordered$sp2.Chr <- factor(mbh_df_levels_reordered$sp2.Chr,levels = levels_sp2)
  
  return(mbh_df_levels_reordered)
}