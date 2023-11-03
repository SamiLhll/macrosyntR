# reorder_multiple_macrosyntenies
#
# This is a function to reorder an orthologs_df, that was generated with load_orthologs(). 
# It retrieves communities using igraph::cluster_fast_greedy

#' @title Reorder the chromosomes of two or more species before plotting
#' @description This is a function to reorder an orthologs_df, same as reorder_macrosynteny, 
#' but it handles tables with more than 2 species.
#'
#' @param orthologs_df dataframe. orthologs with genomic coordinates loaded with load_orthologs()
#'  
#' @return A dataframe object
#' @seealso [load_orthologs()]
#' @seealso [compute_macrosynteny()]
#' @seealso [reorder_macrosynteny()]
#'
#' 
#'@examples 
#' # basic usage of reorder_macrosynteny : 
#' 
#' orthologs_table <- system.file("extdata","my_orthologs.tab",package="macrosyntR")
#' 
#' my_orthologs <- read.table(orthologs_table,header=TRUE)
#' 
#' my_orthologs_reordered <- reorder_multiple_macrosyntenies(my_orthologs)
#' 
#' @export


reorder_multiple_macrosyntenies <- function(orthologs_df) {
  
  amount_of_species <- length(which(colnames(orthologs_df) %in% paste0("sp",c(1:60),".Chr")))
  orthologs_df_to_return <- orthologs_df
  
  for (i in c(1:(amount_of_species-1))) {
    departure_species <- paste0("sp",i)
    arrival_species <- paste0("sp",i+1)
    # subset the dataframe to have only the two species under study :
    departure_colnames <- colnames(orthologs_df)[grep(departure_species,colnames(orthologs_df))]
    arrival_colnames <- colnames(orthologs_df)[grep(arrival_species,colnames(orthologs_df))]
    subset_orthologs_df <- orthologs_df_to_return[c(departure_colnames,arrival_colnames)]
    
    # rename departure_species into sp1 and arrival_species into sp2 :
    temp_colnames <- gsub(departure_species,"sp1",colnames(subset_orthologs_df))
    temp_colnames <- gsub(arrival_species,"sp2",temp_colnames)
    
    colnames(subset_orthologs_df) <- temp_colnames
    # use reorder_macrosynteny to reorder pairwise (keep_sp1_raw_order after the first run)
    if (i > 1) { keep_sp1_raw_order <- TRUE }
    else { keep_sp1_raw_order <- FALSE }
    subset_orthologs_reordered <- reorder_macrosynteny(subset_orthologs_df,
                                                                   keep_sp1_raw_order = keep_sp1_raw_order)
    
    # get only the resulting levels to apply on the spX.Chr of orthologs_df
    reordered_levels_departure <- levels(subset_orthologs_reordered$sp1.Chr)
    reordered_levels_arrival <- levels(subset_orthologs_reordered$sp2.Chr)
    
    orthologs_df_to_return[[paste0(departure_species,".Chr")]] <- factor(orthologs_df_to_return[[paste0(departure_species,".Chr")]], levels = reordered_levels_departure)
    orthologs_df_to_return[[paste0(arrival_species,".Chr")]] <- factor(orthologs_df_to_return[[paste0(arrival_species,".Chr")]],levels = reordered_levels_arrival)
    
    
  }
  
  return(orthologs_df_to_return)
}