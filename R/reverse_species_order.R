# reverse_species_order
#
# This is a function to reverse the species order in an orthologs_df.
#' @title Reverse order of the species in an orthologs_df.
#' @description Returns an orthologs_df (data.frame) with reversed species order compared to the inputted orthologs_df.
#'   sp1 becomes sp2 and the otherway around. It intends at facilitating the integration of more than just two datasets.
#'   It outputs a data.frame shaped as following :
#'   sp1.ID,sp1.Chr,sp1.Start,sp1.End,sp1.Index,sp2.ID,sp2.Chr,sp2.Start,sp2.End,sp2.Index 
#'
#' @param orthologs_df orthologs_df dataframe. mutual best hits with genomic coordinates loaded with load_orthologs()

#' @return dataframe composed of genomic coordinates and relative index of orthologs on both species
#' 
#' @seealso [load_orthologs()]
#' 
#' 
#' @examples 
#' # basic usage of reverse_species_order :
#' 
#' orthologs_table <- system.file("extdata","my_orthologs.tab",package="macrosyntR")
#' 
#' my_orthologs <- read.table(orthologs_table,header=TRUE)
#' 
#' my_orthologs_reversed <- reverse_species_order(my_orthologs)
#' 
#' @export


reverse_species_order <- function(orthologs_df) {

  temp_colnames <- orthologs_df_to_return <- NULL
  
  #### Error check : format of orthologs_df 
  required_fields <- c("sp1.ID","sp1.Index","sp1.Chr","sp2.ID","sp2.Index","sp2.Chr")
  for (i in required_fields) {
    if (isFALSE(i %in% colnames(orthologs_df))) {
      stop("Missing fields in the provided 'orthologs_df'. All the following columns are required : sp1.ID,sp1.Index,sp1.Chr,sp2.ID,sp2.Index,sp2.Chr")
    }
  }
  # Error check : orthologs_df is empty
  if (length(orthologs_df$sp1.Chr) == 0) {stop("Table provided through the 'orthologs_df' argument is empty")}
  
  ### RUN
  # replace sp1 by sp2 and sp2 by sp1 in colnames :
  temp_colnames <- colnames(orthologs_df)
  temp_colnames <- gsub("sp1","spx",temp_colnames)
  temp_colnames <- gsub("sp2","sp1",temp_colnames)
  temp_colnames <- gsub("spx","sp2",temp_colnames)
  
  # copy the orthologs_df and assign the modified colnames in the new table :
  orthologs_df_to_return <- orthologs_df
  colnames(orthologs_df_to_return) <- temp_colnames
  
  return(orthologs_df_to_return)
}