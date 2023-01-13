# load_orthologs
#
# This is a function to load a table with orthologous genes between two species.
#' @title load orthologs with their genomic coordinates.
#' @description Puts together the table of orthologous genes with their genomic coordinates
#'   in the two species under study. It outputs a data.frame shaped as following :
#'   sp1.ID,sp1.Chr,sp1.Start,sp1.End,sp1.Index,sp2.ID,sp2.Chr,sp2.Start,sp2.End,sp2.Index 
#'
#' @param orthologs_table character. Full path to the orthologs table (format : geneID_on_species1   geneID_on_species2)
#' @param sp1_bed character. Full path to the genomic coordinates of the genes on species1 (BED format)
#' @param sp2_bed character. Full path to the genomic coordinates of the genes on species2 (BED format)

#' @return dataframe composed of genomic coordinates and relative index of orthologs on both species
#'
#' @import utils
#' @importFrom dplyr rename mutate group_by arrange ungroup select row_number
#' 
#' @examples 
#' # basic usage of load_orthologs :
#' 
#' orthologs_file <- system.file("extdata","Bflo_vs_Pech.tab",package="macrosyntR")
#' bedfile_sp1 <- system.file("extdata","Bflo.protein_products.bed",package="macrosyntR")
#' bedfile_sp2 <- system.file("extdata","Pech.protein_products.bed",package="macrosyntR")
#' 
#' my_orthologs <- load_orthologs(orthologs_table = orthologs_file,
#'                                sp1_bed = bedfile_sp1,
#'                                sp2_bed = bedfile_sp2)
#' 
#' @export


load_orthologs <- function(orthologs_table,
                           sp1_bed,
                           sp2_bed) {
  
  V1 <- V2 <- V3 <- V4 <- NULL
  sp1.Start <- sp1.End <- sp2.Start <- sp2.End <- sp1.Chr <- sp1.Loci <- sp1.Index <- sp2.Chr <- sp2.Loci <- sp2.Index <- sp1.ID <- sp2.ID <- NULL
  
  # Error check : 1 - species1_bed and species2_bed contains at least 4 fields (tab separated)
  # Error check : 2 - temp_orthologs_table contains two fields (tab separated)

  # open and arrange species1 bedfile :
  species1_bed <- utils::read.csv(sp1_bed, sep = "\t",header = FALSE)
  if (length(species1_bed) < 4) { stop("The bed files must be tab separated and contain at least 4 fields (coordinates and name)")}
  species1_bed <- species1_bed %>%
    dplyr::rename(sp1.Chr = V1, sp1.Start = V2,sp1.End = V3, sp1.ID = V4) %>%
    dplyr::mutate(sp1.Loci = (sp1.Start + sp1.End) /2)
  # open and arrange species2 bedfile :
  species2_bed <- utils::read.csv(sp2_bed, sep = "\t",header = FALSE)
  if (length(species2_bed) < 4) { stop("The bed files must be tab separated and contain at least 4 fields (coordinates and name)")}
  species2_bed <- species2_bed %>%
    dplyr::rename(sp2.Chr = V1, sp2.Start = V2,sp2.End = V3, sp2.ID = V4) %>%
    dplyr::mutate(sp2.Loci = (sp2.Start + sp2.End) /2)
  
  # open orthologs :
  temp_orthologs_table <- utils::read.csv(orthologs_table,sep ="\t", header = FALSE)
  # Error check : the orthologs table must contain two columns
  if(length(temp_orthologs_table) != 2) {stop("The table of orthologs must contain two columns separated by a \"\\t\"")}
  temp_orthologs_table <- temp_orthologs_table %>%
    dplyr::rename(sp1.ID = V1, sp2.ID = V2)
    
  # add the genomic coordinates to the orthologs table :
  temp_orthologs_table <- merge(temp_orthologs_table,species1_bed)
  # Error check : Check that names (4th field of bed files) match the orthologs table
  if (length(temp_orthologs_table$sp1.ID) == 0) {stop("Names (4th column) in sp1_bed don't match with any field in the table of orthologs")}
  orthologs_table_to_return <- merge(temp_orthologs_table,species2_bed)
  # Error check : Check that names (4th field of bed files) match the orthologs table
  if (length(orthologs_table_to_return$sp1.ID) == 0){ stop("Names (4th column) in sp2_bed don't match with any field in the table of orthologs")}
  
  # Calculate the indexes :
  orthologs_table_to_return <- orthologs_table_to_return %>%
    dplyr::group_by(sp1.Chr) %>%
    dplyr::arrange(sp1.Loci) %>%
    dplyr::mutate(sp1.Index = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sp2.Chr) %>%
    dplyr::arrange(sp2.Loci) %>%
    dplyr::mutate(sp2.Index = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(sp1.ID,sp1.Chr,sp1.Start,sp1.End,sp1.Index,sp2.ID,sp2.Chr,sp2.Start,sp2.End,sp2.Index)
  
  orthologs_table_to_return$sp1.Chr <- factor(orthologs_table_to_return$sp1.Chr)
  orthologs_table_to_return$sp2.Chr <- factor(orthologs_table_to_return$sp2.Chr)
    
  
  return(orthologs_table_to_return)
}
