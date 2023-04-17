# compute_macrosynteny
#
# This is a function to generate the contingency table of an orthologs dataframe and apply fischer test to calculate the significant associations.
#' @title Compute significant macrosynteny blocks
#' @description This is a function to generate the contingency table of an orthologs dataframe and 
#'   apply fischer test to calculate the significant associations. It outputs a dataframe shaped as following :
#'   sp1.Chr,sp2.Chr,a,pval,significant,pval_adj
#'
#' @param orthologs_df dataframe. orthologs with genomic coordinates loaded with load_orthologs()
#' @param pvalue_threshold numeric. threshold for significancy. (default equals 0.001)
#' 
#' @return A dataframe object
#'
#' @importFrom stats fisher.test p.adjust
#' @import tidyr
#' @importFrom dplyr rename mutate select
#' 
#' @examples 
#' # basic usage of compute_macrosynteny : 
#' 
#' orthologs_table <- system.file("extdata","my_orthologs.tab",package="macrosyntR")
#' 
#' my_orthologs <- read.table(orthologs_table,header=TRUE)
#'                                
#' my_macrosynteny <- compute_macrosynteny(my_orthologs)
#' 
#' @export


compute_macrosynteny <- function(orthologs_df,pvalue_threshold = 0.001) {
  
  ### construct the contingency table :
  final_i <- final_j <- final_a <- final_b <- final_c <- final_d <- NULL
  sp1.Chr <- sp2.Chr <- value <- variable <- pvalues_value <- pvalues_adj <- odds <- a <- significant <-NULL
  # Error check : proper format for arguments :
  if (!(is.numeric(pvalue_threshold) & length(pvalue_threshold) == 1)) {stop("Wrong format for 'pvalue_threshold' argument. Must be a single value of type numeric")}
  # Error check : format of orthologs_df 
  # required_fields <- c("sp1.ID","sp1.Index","sp1.Chr","sp2.ID","sp2.Index","sp2.Chr")
  required_fields <- c("sp1.ID","sp1.Chr","sp2.ID","sp2.Chr")
  for (i in required_fields) {
    if (isFALSE(i %in% colnames(orthologs_df))) {
      # stop("Missing fields in the provided orthologs_df. All the following columns are required : sp1.ID,sp1.Index,sp1.Chr,sp2.ID,sp2.Index,sp2.Chr")
      stop("Missing fields in the provided orthologs_df. All the following columns are required : sp1.ID,sp1.Chr,sp2.ID,sp2.Chr")
    }
  }
  # Error check : orthologs_df is empty
  if (length(orthologs_df$sp1.Chr) == 0) {stop("Table provided through the orthologs_df argument is empty")}
  # Warning check : when number of chromosomes is too high
  if (length(unique(orthologs_df$sp1.Chr)) >= 300) { 
    warning(paste0("The first species in the orthologs_df has ",length(unique(orthologs_df$sp1.Chr))," chromosomes. Computational time can be very long on fragmented genomes"))
  }
  if (length(unique(orthologs_df$sp2.Chr)) >= 300) { 
    warning(paste0("The second species in the orthologs_df has ",length(unique(orthologs_df$sp2.Chr))," chromosomes. Computational time can be very long on fragmented genomes"))
  }
  
  for (i in unique(orthologs_df$sp1.Chr)) {
    subsetted_orthologs_df_in <- subset(orthologs_df,sp1.Chr == i)
    subsetted_orthologs_df_out <- subset(orthologs_df,sp1.Chr != i)
    compared_specie_scaffs <- unique(subsetted_orthologs_df_in$sp2.Chr)
    for (j in compared_specie_scaffs) {
      cell_a <- length(subset(subsetted_orthologs_df_in, sp2.Chr == j)$sp1.Chr)
      final_a <- c(final_a,cell_a)
      cell_b <- length(subset(subsetted_orthologs_df_out, sp2.Chr == j)$sp1.Chr)
      final_b <- c(final_b,cell_b)
      cell_c <- length(subset(subsetted_orthologs_df_in, sp2.Chr != j)$sp1.Chr)
      final_c <- c(final_c,cell_c)
      cell_d <- length(subset(subsetted_orthologs_df_out, sp2.Chr != j)$sp1.Chr)
      final_d <- c(final_d,cell_d)
      final_j <- c(final_j,j)
    }
    final_i <- c(final_i,rep(i,times=length(compared_specie_scaffs)))
  }
  
  contingency_table <- data.frame(sp1 = final_i,
                                  sp2 = final_j,
                                  a = final_a,b = final_b,
                                  c = final_c,d = final_d)
  ### DONE
  
  rownames(contingency_table) <- paste(contingency_table$sp1,
                                       contingency_table$sp2,sep="-")
  contingency_table <- subset(contingency_table,select=c("a","b","c","d"))
  
  ### Calculate pvalues and odds :
  p_values <- apply(contingency_table, 1,
                    function(x) {
                      tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                      stats::fisher.test(tbl, alternative="greater")$p.value
                    })
  odds_values <- apply(contingency_table, 1,
                       function(x) {
                         tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                         stats::fisher.test(tbl, alternative="greater")$estimate[[1]]
                       })
  
  ### merge pvalues and odds with the contingency table :
  contingency_table$variable <- rownames(contingency_table)
  p_values.melted <- reshape2::melt(p_values) %>%
    dplyr::rename(pvalues_value = value)
  p_values.melted$variable <- rownames(p_values.melted)
  
  odds_values.melted <- reshape2::melt(odds_values) %>%
    dplyr::rename(odds_values = value)
  odds_values.melted$variable <- rownames(odds_values.melted)
  
  contingency_table <- merge(contingency_table,p_values.melted) %>%
    dplyr::mutate(pvalues_adj = stats::p.adjust(pvalues_value)) %>%
    dplyr::mutate(significant = "no",
                  significant = replace(significant,pvalues_adj <= pvalue_threshold,"yes"))
  
  contingency_table <- merge(contingency_table,odds_values.melted) %>%
    dplyr::mutate(odds = "> 1",
                  odds = replace(odds,odds_values < 1,"< 1")) %>%
    tidyr::separate(variable,c("sp1.Chr","sp2.Chr"),sep="-")
  
  ### reshape before returning :
  macrosynt_df <- contingency_table %>%
    dplyr::select(sp1.Chr,sp2.Chr,a,pvalues_value,significant,pvalues_adj) %>%
    dplyr::rename(orthologs = a,pval = pvalues_value,pval_adj = pvalues_adj)
  # copy the levels of orthologs_df to keep the same ordering when plotting :
  if (is.factor(orthologs_df$sp1.Chr)) {
    macrosynt_df$sp1.Chr <- factor(macrosynt_df$sp1.Chr,levels = levels(orthologs_df$sp1.Chr))
  }
  else {
    macrosynt_df$sp1.Chr <- factor(macrosynt_df$sp1.Chr)
  }
  if (is.factor(orthologs_df$sp1.Chr)) {
    macrosynt_df$sp2.Chr <- factor(macrosynt_df$sp2.Chr,levels = levels(orthologs_df$sp2.Chr))
  }
  else {
    macrosynt_df$sp2.Chr <- factor(macrosynt_df$sp2.Chr)
  }
  
  return(macrosynt_df)
  
}