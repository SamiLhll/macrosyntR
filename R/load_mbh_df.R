# load_mbh_df
#
# This is a function to load a table with calculated mutual best hits between two species.
#' @title load Mutual Best Hit table.
#' @description Puts together the mutual best hits calculated beforehand with the genomic coordinates of the 
#'   orthologs on both species to be compared. It outputs a data.frame shaped as following :
#'   sp1.ID,sp1.Chr,sp1.Loci,sp1.Index,sp2.ID,sp2.Chr,sp2.Loci,sp2.Index 
#'
#' @param mbh_table string. Full path to the mutual best hits table (format : geneID_on_species1   geneID_on_species2)
#' @param sp1_bed string. Full path to the genomic coordinates of the genes on species1 (BED format)
#' @param sp2_bed string. Full path to the genomic coordinates of the genes on species2 (BED format)

#' @return dataframe composed of genomic coordinates and relative index of orthologs on both species
#'
#' @import utils
#' @importFrom dplyr rename mutate group_by arrange ungroup select row_number
#' @export


load_mbh_df <- function(mbh_table,
                        sp1_bed,
                        sp2_bed) {
  
  V1 <- V2 <- V3 <- V4 <- NULL
  sp1.Start <- sp1.End <- sp2_start <- sp2.End <- sp1.Chr <- sp1.Loci <- sp1.Index <- sp2.Chr <- sp2.Loci <- sp2.Index <- sp1.ID <- sp2.ID <- NULL
  
  # open mbh output :
  temp_MBH_table <- utils::read.csv(mbh_table,sep ="\t", header = FALSE) %>%
    dplyr::rename(sp1.ID = V1, sp2.ID = V2)
  # open and arrange species1 bedfile :
  species1_bed <- utils::read.csv(sp1_bed, sep = "\t",header = FALSE) %>%
    dplyr::rename(sp1.Chr = V1, sp1.Start = V2,sp1.End = V3, sp1.ID = V4) %>%
    dplyr::mutate(sp1.Loci = (sp1.Start + sp1.End) /2)
  # open and arrange species2 bedfile :
  species2_bed <- utils::read.csv(sp2_bed, sep = "\t",header = FALSE) %>%
    dplyr::rename(sp2.Chr = V1, sp2_start = V2,sp2.End = V3, sp2.ID = V4) %>%
    dplyr::mutate(sp2.Loci = (sp2_start + sp2.End) /2)
  
  # add the genomic coordinates to the MBH table :
  temp_MBH_table <- merge(temp_MBH_table,species1_bed)
  MBH_table_to_return <- merge(temp_MBH_table,species2_bed)
  
  # Calculate the indexes :
  MBH_table_to_return <- MBH_table_to_return %>%
    dplyr::group_by(sp1.Chr) %>%
    dplyr::arrange(sp1.Loci) %>%
    dplyr::mutate(sp1.Index = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sp2.Chr) %>%
    dplyr::arrange(sp2.Loci) %>%
    dplyr::mutate(sp2.Index = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(sp1.ID,sp1.Chr,sp1.Loci,sp1.Index,sp2.ID,sp2.Chr,sp2.Loci,sp2.Index)
  
  MBH_table_to_return$sp1.Chr <- factor(MBH_table_to_return$sp1.Chr)
  MBH_table_to_return$sp2.Chr <- factor(MBH_table_to_return$sp2.Chr)
  return(MBH_table_to_return)
}
