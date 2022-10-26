# load_mbh_df
#
# This is a function to load a table with calculated mutual best hits between two species.
#' @title load Mutual Best Hit table.
#' @description Puts together the mutual best hits calculated beforehand with the genomic coordinates of the orthologs on both species to be compared.
#'
#' @param mbh_table string. Full path to the mutual best hits table (format : geneID_on_species1   geneID_on_species2)
#' @param sp1_bed string. Full path to the genomic coordinates of the genes on species1 (BED format)
#' @param sp2_bed string. Full path to the genomic coordinates of the genes on species2 (BED format)

#' @return dataframe composed of genomic coordinates of orthologs on both species
#'
#' @import utils
#' @import dplyr
#' @export


load_mbh_df <- function(mbh_table,
                        sp1_bed,
                        sp2_bed) {
  
  V1 <- V2 <- V3 <- V4 <- NULL
  sp1_start <- sp1_stop <- sp2_start <- sp2_stop <- sp1_chr <- sp1_mid <- sp2_chr <- sp2_mid <- NULL
  
  # open mbh output :
  temp_MBH_table <- utils::read.delim2(mbh_table,sep ="\t", header = FALSE) %>%
    dplyr::rename(sp1_pep = V1, sp2_pep = V2)
  # open and arrange species1 bedfile :
  species1_bed <- utils::read.delim2(sp1_bed, sep = "\t",header = FALSE) %>%
    dplyr::rename(sp1_chr = V1, sp1_start = V2,sp1_stop = V3, sp1_pep = V4) %>%
    dplyr::mutate(sp1_mid = (sp1_start + sp1_stop) /2)
  # open and arrange species2 bedfile :
  species2_bed <- utils::read.delim2(sp2_bed, sep = "\t",header = FALSE) %>%
    dplyr::rename(sp2_chr = V1, sp2_start = V2,sp2_stop = V3, sp2_pep = V4) %>%
    dplyr::mutate(sp2_mid = (sp2_start + sp2_stop) /2)
  
  # add the genomic coordinates to the MBH table :
  temp_MBH_table <- merge(temp_MBH_table,species1_bed)
  MBH_table_to_return <- merge(temp_MBH_table,species2_bed)
  
  # Calculate the indexes :
  MBH_table_to_return <- MBH_table_to_return %>%
    dplyr::arrange(sp1_chr,sp1_mid) %>%
    dplyr::mutate(sp1_index = row_number()) %>%
    dplyr::select(-sp1_mid) %>%
    dplyr::arrange(sp2_chr,sp2_mid) %>%
    dplyr::mutate(sp2_index = row_number()) %>%
    dplyr::select(-sp2_mid)
  MBH_table_to_return$sp1_chr <- factor(MBH_table_to_return$sp1_chr)
  MBH_table_to_return$sp2_chr <- factor(MBH_table_to_return$sp2_chr)
  return(MBH_table_to_return)
}
