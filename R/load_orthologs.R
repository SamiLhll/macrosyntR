# load_orthologs
#
# This is a function to load a table with orthologous genes between two species.
#' @title load orthologs with their genomic coordinates.
#' @description Puts together the table of orthologous genes with their genomic coordinates
#'   in the two species under study. It outputs a data.frame shaped as following :
#'   sp1.ID,sp1.Chr,sp1.Start,sp1.End,sp1.Index,sp2.ID,sp2.Chr,sp2.Start,sp2.End,sp2.Index 
#'
#' @param orthologs_table string. Full path to the orthologs table (format : geneID_on_species1   geneID_on_species2)
#' @param sp1_bed string. Full path to the genomic coordinates of the genes on species1 (BED format)
#' @param sp2_bed string. Full path to the genomic coordinates of the genes on species2 (BED format)

#' @return dataframe composed of genomic coordinates and relative index of orthologs on both species
#'
#' @import utils
#' @importFrom dplyr rename mutate group_by arrange ungroup select row_number
#' @export


load_orthologs <- function(orthologs_table,
                           sp1_bed,
                           sp2_bed) {
  
  V1 <- V2 <- V3 <- V4 <- NULL
  sp1.Start <- sp1.End <- sp2.Start <- sp2.End <- sp1.Chr <- sp1.Loci <- sp1.Index <- sp2.Chr <- sp2.Loci <- sp2.Index <- sp1.ID <- sp2.ID <- NULL
  
  # open orthologs :
  temp_orthologs_table <- utils::read.csv(orthologs_table,sep ="\t", header = FALSE) %>%
    dplyr::rename(sp1.ID = V1, sp2.ID = V2)
  # open and arrange species1 bedfile :
  species1_bed <- utils::read.csv(sp1_bed, sep = "\t",header = FALSE) %>%
    dplyr::rename(sp1.Chr = V1, sp1.Start = V2,sp1.End = V3, sp1.ID = V4) %>%
    dplyr::mutate(sp1.Loci = (sp1.Start + sp1.End) /2)
  # open and arrange species2 bedfile :
  species2_bed <- utils::read.csv(sp2_bed, sep = "\t",header = FALSE) %>%
    dplyr::rename(sp2.Chr = V1, sp2.Start = V2,sp2.End = V3, sp2.ID = V4) %>%
    dplyr::mutate(sp2.Loci = (sp2.Start + sp2.End) /2)
  
  # add the genomic coordinates to the orthologs table :
  temp_orthologs_table <- merge(temp_orthologs_table,species1_bed)
  orthologs_table_to_return <- merge(temp_orthologs_table,species2_bed)
  
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