# calculate_contingency_table
#
# This is a function to generate the contingency table of an MBH dataframe and apply fischer test to calculate the significant associations.
#' @title Calculate contingency table
#' @description This is a function to generate the contingency table of an MBH dataframe and apply fischer test to calculate the significant associations.
#'
#' @param MBH_table dataframe of mutual best hits with genomic coordinates loaded by the load_MBH_table function
#'
#' @return A gdataframe object
#'
#' @import dplyr
#' @export


calculate_contingency_table <- function(MBH_table) {

  # construct the contingency table :
  final_i <- NULL
  final_j <- NULL
  final_a <- NULL
  final_b <- NULL
  final_c <- NULL
  final_d <- NULL

  for (i in unique(MBH_table$sp1_chr)) {
    subsetted_MBH_table_in <- subset(MBH_table,sp1_chr == i)
    subsetted_MBH_table_out <- subset(MBH_table,sp1_chr != i)
    compared_specie_scaffs <- unique(subsetted_MBH_table_in$sp2_chr)
    for (j in compared_specie_scaffs) {
      cell_a <- length(subset(subsetted_MBH_table_in, sp2_chr == j)$sp1_chr)
      final_a <- c(final_a,cell_a)
      cell_b <- length(subset(subsetted_MBH_table_out, sp2_chr == j)$sp1_chr)
      final_b <- c(final_b,cell_b)
      cell_c <- length(subset(subsetted_MBH_table_in, sp2_chr != j)$sp1_chr)
      final_c <- c(final_c,cell_c)
      cell_d <- length(subset(subsetted_MBH_table_out, sp2_chr != j)$sp1_chr)
      final_d <- c(final_d,cell_d)
      final_j <- c(final_j,j)
    }
    final_i <- c(final_i,rep(i,times=length(compared_specie_scaffs)))
  }

  contingency_table <- data.frame(sp1 = final_i,
                                  sp2 = final_j,
                                  a = final_a,b = final_b,
                                  c = final_c,d = final_d)

  rownames(contingency_table) <- paste(contingency_table$sp1,contingency_table$sp2,sep="-")
  contingency_table <- subset(contingency_table,select=c("a","b","c","d"))

  ### Calculate pvalues and odds :
  p_values <- apply(contingency_table, 1,
                    function(x) {
                      tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                      fisher.test(tbl, alternative="greater")$p.value
                    })
  odds_values <- apply(contingency_table, 1,
                       function(x) {
                         tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
                         fisher.test(tbl, alternative="greater")$estimate[[1]]
                       })

  ### merge pvalues and odds with the contingency table :
  contingency_table$variable <- rownames(contingency_table)
  p_values.melted <- reshape2::melt(p_values) %>%
    rename(pvalues_value = value)
  p_values.melted$variable <- rownames(p_values.melted)

  odds_values.melted <- reshape2::melt(odds_values) %>%
    rename(odds_values = value)
  odds_values.melted$variable <- rownames(odds_values.melted)
  head(odds_values.melted)

  contingency_table <- merge(contingency_table,p_values.melted) %>%
    mutate(significant = "yes",
           significant = replace(significant,pvalues_value > 0.001,"no"))

  contingency_table <- merge(contingency_table,odds_values.melted) %>%
    mutate(odds = "> 1",
           odds = replace(odds,odds_values < 1,"< 1")) %>%
    separate(variable,c("sp1_chr","sp2_chr"),sep="-")

  return(contingency_table)

}
