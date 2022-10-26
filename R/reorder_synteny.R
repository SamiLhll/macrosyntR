# reorder_syneny
#
# This is a function to reorder an mbh_df, that was generated with load_mbh_df(). It allows for automatic reordering, or manual filtering and reordering.
#' @title Reorder the mbh_df before plotting
#' @description his is a function to reorder an mbh_df, that was generated with load_mbh_df(). It allows for automatic reordering, or manual filtering and reordering.
#'
#' @param mbh_df dataframe. mutual best hits with genomic coordinates loaded with load_mbh_df()
#' @param auto logical. Tells if the reodering must be done automatically. When plotting the resulting reordered df, you'll have the macrosynteny groups ordered by decreasing size from left to right, top to bottom (default = TRUE)
#' @param sp1_order vector of strings. (default = NULL) For manual reordering, change the auto parameter to FALSE and list the ordered chromsomes of species 1 that you want to keep such as c("chrX","chr3",...,"chrY")
#' @param sp2_order vector of strings. (default = NULL) For manual reordering, change the auto parameter to FALSE and list the ordered chromsomes of species 2 that you want to keep such as c("chrX","chr3",...,"chrY")
#' @return A dataframe object
#'
#' @import dplyr
#' 
#' @export


reorder_synteny <- function(mbh_df,
                            auto= TRUE,
                            sp1_order= NULL,
                            sp2_order = NULL) {
  
  sp1_chr <- sp2_chr <- significant <- orthologs <- NULL
  
  # 1 - case of manual filtering :
  if (! auto) {
    mbh_df_levels_reordered <- mbh_df
    if (! is.null(sp1_order)) {
      mbh_df_levels_reordered <- subset(mbh_df_levels_reordered,sp1_chr %in% sp1_order)
      mbh_df_levels_reordered$sp1_chr <- factor(mbh_df_levels_reordered$sp1_chr,levels = sp1_order)
    }
    if (! is.null(sp2_order)){
      mbh_df_levels_reordered <- subset(mbh_df_levels_reordered,sp2_chr %in% sp2_order)
      mbh_df_levels_reordered$sp2_chr <- factor(mbh_df_levels_reordered$sp2_chr,levels = sp2_order)
    }
  }
  else {
    # 2 - sort the sp1 chr by decreasing order of macrosynteny group sizes :
    contingency_table <- calculate_macrosynt(mbh_df)
    significant_entries <- subset(contingency_table, significant == "yes") %>%
      dplyr::arrange(desc(orthologs))
    
    canva_levels_sp1 <- unique(significant_entries$sp1_chr)
    mbh_df_levels_reordered <- subset(mbh_df,sp1_chr %in% canva_levels_sp1)
    
    # 3 - sort the sp2 chr :
    levels_sp1 <- NULL
    levels_sp2 <- NULL
    for (i in canva_levels_sp1) {
      # get all sp2 chromosomes significant with sp1 :
      temp_significant_sp2_1st_round <- subset(significant_entries,sp1_chr == i)
      # get all sp1 chromosomes significant with these sp2 :
      temp_significant_sp1_1st_round <- subset(significant_entries,sp2_chr %in% as.character(temp_significant_sp2_1st_round$sp2_chr))
      # additional filtering of sp2 chromosomes that are already placed in the table :
      temp_significant_sp2_2nd_round <- temp_significant_sp2_1st_round[ !(temp_significant_sp2_1st_round$sp2_chr %in% levels_sp2),]
      temp_significant_sp1_2nd_round <- temp_significant_sp1_1st_round[ !(temp_significant_sp1_1st_round$sp1_chr %in% levels_sp1),]
      levels_sp2 <- c(levels_sp2,unique(as.character(temp_significant_sp2_2nd_round$sp2_chr)))
      levels_sp1 <- c(levels_sp1,unique(as.character(temp_significant_sp1_2nd_round$sp1_chr)))
    }
    mbh_df_levels_reordered <- subset(mbh_df_levels_reordered, sp2_chr %in% levels_sp2)
    mbh_df_levels_reordered$sp1_chr <- factor(mbh_df_levels_reordered$sp1_chr,levels = levels_sp1)
    mbh_df_levels_reordered$sp2_chr <- factor(mbh_df_levels_reordered$sp2_chr,levels = levels_sp2)
  }
  return(mbh_df_levels_reordered)
}