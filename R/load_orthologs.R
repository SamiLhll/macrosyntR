# load_orthologs
#
# This is a function to load a table with orthologous genes between two or more species.
#' @title load orthologs with their genomic coordinates.
#' @description Puts together the table of orthologous genes with their genomic coordinates
#'   in the two or more species. It outputs a data.frame shaped as following :
#'   sp1.ID,sp1.Chr,sp1.Start,sp1.End,sp1.Index,sp2.ID,sp2.Chr,sp2.Start,sp2.End,sp2.Index,... 
#'
#' @param orthologs_table character. Full path to the orthologs table (format : geneID_on_species1   geneID_on_species2   geneID_on_speciesN)
#' @param sp1_bed (deprecated) character. Full path to the genomic coordinates of the genes on species1
#' @param sp2_bed (deprecated) character. Full path to the genomic coordinates of the genes on species2
#' @param bedfiles array. List of full paths to the genomic coordinates ordered as in the appearing order of the orthologs_table (BED format)

#' @return dataframe composed of genomic coordinates and relative index of orthologs on both species
#'
#' @import utils
#' @importFrom dplyr rename mutate group_by arrange ungroup select row_number
#' 
#' @examples 
#' # basic usage of load_orthologs for two species :
#' 
#' orthologs_file <- system.file("extdata","Bflo_vs_Pyes.tab",package="macrosyntR")
#' bedfile_sp1 <- system.file("extdata","Bflo.bed",package="macrosyntR")
#' bedfile_sp2 <- system.file("extdata","Pyes.bed",package="macrosyntR")
#'
#' 
#' my_orthologs <- load_orthologs(orthologs_table = orthologs_file,
#'                                bedfiles = c(bedfile_sp1,bedfile_sp2))
#'# example with 3 species :
#' orthologs_file <- system.file("extdata","Single_copy_orthologs.tsv",package="macrosyntR")
#' bedfile_sp3 <- system.file("extdata","Pech.bed",package="macrosyntR")
#' 
#' my_orthologs <- load_orthologs(orthologs_table = orthologs_file,
#'                                bedfiles = c(bedfile_sp1,bedfile_sp2,bedfile_sp3))
#' 
#' 
#' @export


load_orthologs <- function(orthologs_table,
                           sp1_bed = NULL,
                           sp2_bed = NULL,
                           bedfiles = NULL) {
  
  V1 <- V2 <- V3 <- V4 <- NULL
  Start <- End <- Chr <- Loci <- Index <- ID <- NULL
  
  # Error check : 1 - species1_bed and species2_bed contains at least 4 fields (tab separated)
  # Error check : 2 - temp_orthologs_table contains two fields (tab separated)
  
  
  #### 1 - Open orthologs :
  orthologs_table_to_return <- utils::read.csv(orthologs_table,sep ="\t", header = FALSE)
  # Error check : the orthologs table must contain two columns
  if(length(orthologs_table_to_return) < 2) {stop("The table of orthologs must contain at least two columns separated by a \"\\t\"")}
  
  # rename columns with "spx.ID"
  number_of_species <- length(orthologs_table_to_return)
  colnames(orthologs_table_to_return) <- paste0("sp",c(1:number_of_species),".ID")
  ####
  
  #### 2 - Open bedfiles :
  
  ## When species1_bed and species2_bed are used :
  if (!is.null(sp1_bed)) { 
    if (!is.null(sp2_bed)) {
      bedfiles <- c(sp1_bed,sp2_bed)
    }
    else {
      stop("Error. Arguments sp1_bed and sp2_bed are meant to be used together")
    }
  }
  else if (!is.null(sp2_bed)) {
    stop("Error. Arguments sp1_bed and sp2_bed are meant to be used together")
  }
  
  ##
  if (!is.null(bedfiles)) {
    current_species <- 0
    ## Check amount of bedfiles matchs amount of columns in orthologs_table :
    if (length(orthologs_table_to_return) != length(bedfiles)) {
      stop("Error. Amount of species in orthologs_table differs from the amount of provided bedfiles.")
    }
    ##
    for (current_bedfile_path in bedfiles) {
      current_species <- current_species + 1
      # open bedfile :
      current_bedfile <- utils::read.csv(current_bedfile_path,sep = "\t", header = FALSE) 
      # check conformity of bedfile :
      if (length(current_bedfile) < 4) { stop(paste0("Error. Bedfiles must be tab separated containing at least 4 fields. Wrong format in : ",current_bedfile_path))}
      # compute indexes :
      current_bedfile <- current_bedfile %>%
        rename(Chr = V1, Start = V2, End = V3, ID = V4) %>%
        mutate(Loci = (Start + End)/2) %>%
        group_by(Chr) %>%
        arrange(Loci) %>%
        mutate(Index = dplyr::row_number()) %>%
        ungroup() %>%
        select(ID,Chr,Start,End,Index)
      # Convert Chr to factor
      current_bedfile$Chr <- factor(current_bedfile$Chr)
      # rename columns with species identifier :
      species_prefix <- paste0("sp",current_species,".")
      colnames(current_bedfile) <- paste0(species_prefix,colnames(current_bedfile))
      # merge with table of orthologs
      orthologs_table_to_return <- merge(orthologs_table_to_return,current_bedfile)
    }
  }
  return(orthologs_table_to_return)
}
