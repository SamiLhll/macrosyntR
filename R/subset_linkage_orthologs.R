# subset_linkage_orthologs
#
# This is a function to subset_the orthologs that are contained in the linkage groups of an orthologs_df
#' @title Subset Orthologs contained in conserved linkage groups
#' @description This is a function to subset an orthologs_df and keep only the orthologs that are within significant
#' linkage groups computed by the function compute_linkage_groups().
#'
#' @param orthologs_df dataframe. orthologs with genomic coordinates loaded with load_orthologs()
#' @param linkages dataframe. table listing the linkage groups as returned by the function compute_linkage_groups()
#' 
#' @seealso [load_orthologs()]
#' @seealso [compute_linkage_groups()]
#' 
#' @return A dataframe object
#'
#' @import tidyr
#' @importFrom dplyr group_by arrange mutate row_number ungroup
#' 
#' @examples 
#' # basic usage of compute_linkage_groups: 
#' 
#' orthologs_table <- system.file("extdata","my_orthologs.tab",package="macrosyntR")
#' 
#' my_orthologs <- read.table(orthologs_table,header=TRUE)
#'                                
#' my_macrosynteny <- compute_linkage_groups(my_orthologs)
#' 
#' @export



subset_linkage_orthologs <- function(orthologs_df,
                                     linkages=NULL) {
  
  sp1.Chr <- sp1.Start <- sp2.Chr <- sp2.Start <- NULL
  
  if (is.null(linkages)) {
    linkages <- compute_linkage_groups(orthologs_df)
  }
  
  orthologs_df_to_return <- merge(orthologs_df,linkages)
  
  # recompute indexes :
  orthologs_df_to_return <- orthologs_df_to_return %>%
    group_by(sp1.Chr) %>% arrange(sp1.Start) %>% mutate(sp1.Index = dplyr::row_number()) %>% ungroup() %>%
    group_by(sp2.Chr) %>% arrange(sp2.Start) %>% mutate(sp2.Index = dplyr::row_number()) %>% ungroup()
  
  
  return(orthologs_df_to_return)
}
