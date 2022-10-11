# load_MBH_table
#
# This is a function to load a table with calculated mutual best hits between two species.
#' @title load Mutual Best Hit table.
#' @description Puts together the mutual best hits calculated beforehand with the genomic coordinates of the orthologs on both species to be compared.
#'
#' @param MBH_table_path path to the mutual best hits table (format : geneID_on_species1   geneID_on_species2)
#' @param species1_bed_path path to the genomic coordinates of the genes on species1 (BED format)
#' @param species2_bed_path path to the genomic coordinates of the genes on species2 (BED format)
#' @param sp1_keep_n_chr integer to tell the n chromosomes of species1 with the higher amount of orthologs to keep
#' @param sp2_keep_n_chr integer to tell the n chromosomes of species2 with the higher amount of orthologs to keep
#' @param sp1_keep_chr_names list of chr names to keep on species1
#' @param sp2_keep_chr_names list of chr names to keep on species 2
#'
#' @return a dataframe composed of genomic coordinates of orthologs on both species
#'
#' @import dplyr
#' @export


load_MBH_table <- function(MBH_table_path,species1_bed_path,species2_bed_path,sp1_keep_n_chr = 0,sp2_keep_n_chr = 0,sp1_keep_chr_names=NULL,sp2_keep_chr_names=NULL) {
  # load Generate_MBH output :
  temp_MBH_table <- read.csv(MBH_table_path,sep ="\t", header = FALSE) %>%
    rename(sp1_pep = V1, sp2_pep = V2)
  # load and arrange species1 bedfile :
  species1_bed <- read.csv(species1_bed_path, sep = "\t",header = FALSE) %>%
    rename(sp1_chr = V1, sp1_start = V2,sp1_stop = V3, sp1_pep = V4) %>%
    mutate(sp1_mid = (sp1_start + sp1_stop) /2)
  # load and arrange species2 bedfile :
  species2_bed <- read.csv(species2_bed_path, sep = "\t",header = FALSE) %>%
    rename(sp2_chr = V1, sp2_start = V2,sp2_stop = V3, sp2_pep = V4) %>%
    mutate(sp2_mid = (sp2_start + sp2_stop) /2)

  # add the genomic coordinates to the MBH table :
  temp_MBH_table <- merge(temp_MBH_table,species1_bed)
  MBH_table_to_return <- merge(temp_MBH_table,species2_bed)

  # Filter to keep only longest scaffolds if necessary
  if (sp1_keep_n_chr > 0) {
    temp_summary <- MBH_table_to_return %>% group_by(sp1_chr) %>% summarise( n = n()) %>% arrange(desc(n))
    temp_scaffolds_to_keep <- temp_summary$sp1_chr[1:sp1_keep_n_chr]
    MBH_table_to_return <- subset(MBH_table_to_return,sp1_chr %in% temp_scaffolds_to_keep)
  }
  if (sp2_keep_n_chr > 0) {
    temp_summary <- MBH_table_to_return %>% group_by(sp2_chr) %>% summarise( n = n()) %>% arrange(desc(n))
    temp_scaffolds_to_keep <- temp_summary$sp2_chr[1:sp2_keep_n_chr]
    MBH_table_to_return <- subset(MBH_table_to_return,sp2_chr %in% temp_scaffolds_to_keep)
  }
  # Filter to keep only the specified chromosome list
  if (! is.null(sp1_keep_chr_names)) {
    MBH_table_to_return <- subset(MBH_table_to_return,sp1_chr %in% sp1_keep_chr_names)
  }
  if (! is.null(sp2_keep_chr_names)) {
    MBH_table_to_return <- subset(MBH_table_to_return,sp2_chr %in% sp2_keep_chr_names)
  }


  # Calculate the indexes :
  MBH_table_to_return <- MBH_table_to_return %>%
    arrange(sp1_chr,sp1_mid) %>%
    mutate(sp1_index = row_number()) %>%
    select(-sp1_mid) %>%
    arrange(sp2_chr,sp2_mid) %>%
    mutate(sp2_index = row_number()) %>%
    select(-sp2_mid)

  return(MBH_table_to_return)
}
