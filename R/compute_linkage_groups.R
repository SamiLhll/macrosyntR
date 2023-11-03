# compute_linkage_groups
#
# This is a function to compute the linkage groups of an table of orthologs that includes two or more species.
#' @title Compute Linkage groups
#' @description This is a function to compute the conserved linkage groups shared between two or more species.
#' It computes the significant associations between chromosomes of all species versus all (pairwise) 
#' using the fischer test implemented in compute_macrosynteny().  
#' It outputs a dataframe shaped as following :
#'   sp1.Chr,sp2.Chr,..., spN.chr,n,LGs
#'  where n is the number of shared orthologs in the group and LGs are the IDs for the linkage groups
#'
#' @param orthologs_df dataframe. orthologs with genomic coordinates loaded with load_orthologs()
#' 
#' @return A dataframe object
#'
#' @import tidyr
#' @importFrom dplyr group_by summarise ungroup arrange n
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


compute_linkage_groups <- function(orthologs_df) {
  
  
  sp1.Chr <- sp2.Chr <- significant <- n <- NULL
  
  amount_of_species <- length(which(colnames(orthologs_df) %in% paste0("sp",c(1:60),".Chr")))
  
  options(dplyr.summarise.inform = FALSE)
  all_linkages <- NULL
  # Iterate and compare every species two by two :
  for (i in c(1:amount_of_species)) {
    if (i < amount_of_species) {
      departure_species <- paste0("sp",i)
      for (j in c((i+1):amount_of_species)) {
        arrival_species <- paste0("sp",j)
        # subset orthologs_df to get only the two species under study :
        departure_colnames <- colnames(orthologs_df)[grep(departure_species,colnames(orthologs_df))]
        arrival_colnames <- colnames(orthologs_df)[grep(arrival_species,colnames(orthologs_df))]
        subset_orthologs_df <- orthologs_df[c(departure_colnames,arrival_colnames)]
        # rename to sp1 and sp2 :
        colnames(subset_orthologs_df) <- gsub(departure_species,"sp1",colnames(subset_orthologs_df))
        colnames(subset_orthologs_df) <- gsub(arrival_species,"sp2",colnames(subset_orthologs_df))
        temp_macrosynteny <- compute_macrosynteny(subset_orthologs_df) %>%
          select(sp1.Chr,sp2.Chr,significant)
        temp_macrosynteny <- subset(temp_macrosynteny,significant == "yes") %>%
          select(-significant)
        colnames(temp_macrosynteny) <- c(paste0(departure_species,".Chr"),paste0(arrival_species,".Chr"))
        if (i == 1 & j == 2) { all_linkages <- temp_macrosynteny }
        else { all_linkages <- merge(all_linkages,temp_macrosynteny) }
      }
    }
  }
  # Compute amount of orthologs
  #colnames to group_by :
  grp_cols <- colnames(all_linkages)
  dots <- lapply(grp_cols, as.symbol)
  linkages_to_return <- merge(orthologs_df,all_linkages) %>%
    group_by(.groups = !!!syms(grp_cols)) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    arrange(desc(n))
  # mutate(linkage = as.character(dplyr::row_number()))
  letters702 <- c(letters, sapply(letters, function(x) paste0(x, letters)))
  linkages_to_return$LGs <- letters702[1:length(linkages_to_return$sp1.Chr)]

  
  options(dplyr.summarise.inform = TRUE)
  
  return(linkages_to_return)
}