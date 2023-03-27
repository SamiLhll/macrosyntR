# get_syntenic_genes
#
# This is a function to extract all the syntenic genes from an orthologs_df.
#' @title get the syntenic genes as a table
#' @description This is a function to extract all the syntenic genes from an orthologs_df. It requires as input an orthologs_df loaded by load_orthologs(). 
#'
#' @param orthologs_df dataframe. orthologs with genomic coordinates loaded by load_orthologs()
#' 
#' @importFrom dplyr select
#' 
#' @return dataframe composed of details for each detected syntenic block of genes. It contains the following columns : sp1.Chr, sp1.Start, sp1.End, sp2.Chr, sp2.Start, sp2.End, size, sp1.IDs, sp2.IDs
#' 
#' @seealso [load_orthologs()]
#' 
#' @examples 
#' # basic usage of get_syntenic_genes :
#' 
#' orthologs_file <- system.file("extdata","Bflo_vs_Pech.tab",package="macrosyntR")
#' bedfile_sp1 <- system.file("extdata","Bflo.protein_products.bed",package="macrosyntR")
#' bedfile_sp2 <- system.file("extdata","Pech.protein_products.bed",package="macrosyntR")
#' 
#' my_orthologs <- load_orthologs(orthologs_table = orthologs_file,
#'                                sp1_bed = bedfile_sp1,
#'                                sp2_bed = bedfile_sp2)
#'                                
#' my_syntenic_block_of_genes <- get_syntenic_genes(my_orthologs)
#' 
#' @export


get_syntenic_genes <- function(orthologs_df) {
  
  ##################################################
  ####### ERROR CHECK : 
  # Error check : format of orthologs_df 
  required_fields <- c("sp1.ID","sp1.Index","sp1.Chr","sp2.ID","sp2.Index","sp2.Chr")
  for (i in required_fields) {
    if (isFALSE(i %in% colnames(orthologs_df))) {
      required_fields_character <- paste(required_fields,sep=",")
      stop("Missing fields in the provided orthologs_df. All the following columns are required : sp1.ID,sp1.Index,sp1.Chr,sp2.ID,sp2.Index,sp2.Chr")
    }
  }
  # Error check : orthologs_df is empty
  if (length(orthologs_df$sp1.Chr) == 0) {stop("Table provided through the orthologs_df argument is empty")}
  ##################################################
  
  Syntenic_blocks <- NULL
  size <- sp1.Chr <- sp1.End <- sp1.IDs <- sp1.Index <- sp1.Start <- sp2.Chr <- sp2.End <- sp2.IDs <- sp2.Start <- NULL
  
  chromosomes_sp1 <- unique(orthologs_df$sp1.Chr)
  # Iterate on all chromosomes from 1st species
  for (i in (chromosomes_sp1)) {
    chr_table <- subset(orthologs_df,sp1.Chr == i)
    chromosome_sp2 <- unique(chr_table$sp2.Chr)
    # Iterate on all chromosomes from 2nd species
    for (ii in (chromosome_sp2)) {
      # get the subset of orthologs df and arrange the genes by the relative order in species 1 :
      chr_table <- subset(orthologs_df,sp1.Chr == i & sp2.Chr == ii) %>%
        arrange(sp1.Index)
      # Compute the distance (relative) of genes in species 2 following the order in species 1 : 
      distance_to_next_gene <- diff(chr_table$sp2.Index)
      distance_to_next_gene[abs(distance_to_next_gene) != 1] <- 0
      distance_to_next_gene <- abs(distance_to_next_gene)
      
      ### Start to walk through the list of distances to get all the syntenic genes
      Index_on_list <- 1
      while (Index_on_list < length(distance_to_next_gene)) {
        current_pair <- distance_to_next_gene[Index_on_list]
        if(current_pair == 1) {
          new_syntenic_block <- c(Index_on_list,Index_on_list + 1)
          Index_on_list <- Index_on_list + 1
          current_pair <- distance_to_next_gene[Index_on_list]
          while(current_pair == 1 & Index_on_list < length(distance_to_next_gene)) {
            new_syntenic_block <- c(new_syntenic_block, Index_on_list,Index_on_list + 1)
            Index_on_list <- Index_on_list + 1
            current_pair <- distance_to_next_gene[Index_on_list]
          }
          ## Build the piece of dataframe for the detected syntenic genes :
          new_syntenic_block <- unique(new_syntenic_block)
          genes_sp1_in_syntenic_block <- toString(chr_table$sp1.ID[new_syntenic_block],)
          genes_sp2_in_syntenic_block <- toString(chr_table$sp2.ID[new_syntenic_block],)
          start_sp1 <- min(chr_table$sp1.Start[new_syntenic_block])
          stop_sp1 <- max(chr_table$sp1.End[new_syntenic_block])
          start_sp2 <- min(chr_table$sp2.Start[new_syntenic_block])
          stop_sp2 <- max(chr_table$sp2.End[new_syntenic_block])
          temp_syntenic_df <- data.frame(sp1.IDs = genes_sp1_in_syntenic_block,
                                         sp2.IDs = genes_sp2_in_syntenic_block) %>% 
            mutate(sp1.Chr = i,
                   sp2.Chr = ii,
                   sp1.Start = start_sp1,
                   sp1.End = stop_sp1,
                   sp2.Start = start_sp2,
                   sp2.End = stop_sp2,
                   size = length(chr_table$sp1.ID[new_syntenic_block]))
          Syntenic_blocks <- rbind(Syntenic_blocks,temp_syntenic_df) 
        }
        else {
          Index_on_list <- Index_on_list +1
        }
      }
      
      Syntenic_blocks <- Syntenic_blocks %>%
        select(sp1.Chr,sp1.Start, sp1.End,sp2.Chr,sp2.Start,sp2.End,size,sp1.IDs,sp2.IDs)
    }
  }
  return(Syntenic_blocks)
}