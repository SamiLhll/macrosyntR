# reorder_syneny
#
# This is a function to reorder an mbh_df, that was generated with load_mbh_df(). It allows for automatic reordering, or manual filtering and reordering.
#' @title Reorder the mbh_df before plotting
#' @description his is a function to reorder an mbh_df, that was generated with load_mbh_df(). It allows for automatic reordering, or manual filtering and reordering.
#'
#' @param mbh_df dataframe. mutual best hits with genomic coordinates loaded with load_mbh_df()
#' @param pvalue_threshold numeric. threshold for significancy. (default equals 0.001)
#' @return A dataframe object
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph components
#' @importFrom igraph groups
#' @importFrom dplyr select rename arrange group_by summarise ungroup mutate setdiff
#' 
#' @export


reorder_synteny <- function(mbh_df,
                            pvalue_threshold = 0.001) {
  
  contingency_table <- significant_entries <- significant_entries_for_graph <- NULL
  significant_entries_for_graph <- significant_association_undirected_graph <- clusters <- cluster_df <- NULL
  sp1_amounts <- sp2_amounts <- NULL
  final_levels_sp1 <- final_levels_sp2 <- NULL
  significant <- pval <- orthologs <- sp1.Chr <- sp2.Chr <- amounts <- clust <-NULL
  
  ##### 1 - Build an undirected and unweighted graph of connected chromosomes (significant amount of orthologs)
  # Get a table with only significant association using calculate_macrosynt from this package :
  contingency_table <- calculate_macrosynt(mbh_df,pvalue_threshold)
  significant_entries <- subset(contingency_table,significant == "yes")
  
  # build an Undirected Graph between sp1.chr and sp2.chr containing all significant edges :
  significant_entries_for_graph <- significant_entries %>% dplyr::select(-significant,-pval,-orthologs) %>% dplyr::rename(from = sp1.Chr,to = sp2.Chr)
  # compute clusters of sp1.Chr and sp2.Chr that are directly or indirectly connected in the graph :
  significant_association_undirected_graph <- igraph::graph_from_data_frame(significant_entries_for_graph, directed = F)
  clusters <- igraph::components(significant_association_undirected_graph)
  ##### DONE : Built the graph
  
  ##### 2 - Compute amounts of orthologs in each cluster and order them by attributing them an index from 1 to n (1 being the larger cluster) :
  chroms_in_cluster_as_string <- NULL
  ortholog_amount_in_cluster <- NULL
  for (list_of_chrom_in_cluster in igraph::groups(clusters)) {
    temp_amounts <- subset(significant_entries,(sp1.Chr %in% list_of_chrom_in_cluster) & (sp2.Chr %in% list_of_chrom_in_cluster))
    chroms_in_cluster_as_string <- c(chroms_in_cluster_as_string,paste(list_of_chrom_in_cluster,collapse = ","))
    ortholog_amount_in_cluster <- c(ortholog_amount_in_cluster,sum(temp_amounts$orthologs))
  }
  cluster_df <- data.frame(clust = chroms_in_cluster_as_string,
                           amounts = ortholog_amount_in_cluster) %>%
    dplyr::arrange(desc(amounts))
  cluster_df$num_clust <- seq(1:length(cluster_df$amounts))
  ##### DONE : Computed order of clusters
  
  ##### 3 - Compute order of chromosomes by recomputing the levels of the factor. sorting by descending cluster size then descencding chromosome size.
  # get a table with number of orthologs for each chromosome
  sp1_amounts <- significant_entries %>% dplyr::group_by(sp1.Chr) %>% dplyr::summarise(n=sum(orthologs)) %>% dplyr::ungroup()
  sp2_amounts <- significant_entries %>% dplyr::group_by(sp2.Chr) %>% dplyr::summarise(n=sum(orthologs)) %>% dplyr::ungroup()
  chrom_clusters_reordered <- cluster_df$clust
  for(chrom_cluster in chrom_clusters_reordered) {
    chrom_cluster_list <- strsplit(chrom_cluster,",")[[1]]
    cluster_levels_sp1 <- NULL
    cluster_levels_sp2 <- NULL
    for (chrom in chrom_cluster_list) {
      if (chrom %in% mbh_df$sp1.Chr) {
        cluster_levels_sp1 <- c(cluster_levels_sp1,chrom)
      }
      else {
        cluster_levels_sp2 <- c(cluster_levels_sp2,chrom)
      }
    }
    # order the chromosomes from larger to smaller (in amount of orthologs) :
    cluster_levels_sp1_df <- subset(sp1_amounts,sp1.Chr %in% cluster_levels_sp1) %>%
      dplyr::arrange(desc(n))
    final_levels_sp1 <- c(final_levels_sp1,as.character(cluster_levels_sp1_df$sp1.Chr))
    cluster_levels_sp2_df <- subset(sp2_amounts,sp2.Chr %in% cluster_levels_sp2) %>%
      dplyr::arrange(desc(n))
    final_levels_sp2 <- c(final_levels_sp2,as.character(cluster_levels_sp2_df$sp2.Chr))
    
  }
  mbh_df_to_return <- subset(mbh_df, (sp1.Chr %in% final_levels_sp1) & (sp2.Chr %in% final_levels_sp2))
  mbh_df_to_return$sp1.Chr <- factor(mbh_df_to_return$sp1.Chr,levels = final_levels_sp1)
  mbh_df_to_return$sp2.Chr <- factor(mbh_df_to_return$sp2.Chr,levels = final_levels_sp2)
  ##### DONE : Computed order of chromosomes
  
  ##### 4 - Add an additional column with clusters number :
  final_df_with_groups <- NULL
  cluster_num <- 0
  # cluster_df calculated in part 2 of this function
  chrom_clusters_reordered <- cluster_df$clust
  for (i in chrom_clusters_reordered){
    chrom_cluster_list <- strsplit(i,",")[[1]]
    cluster_num <- cluster_num + 1
    temp <- subset(mbh_df_to_return,((sp1.Chr %in% chrom_cluster_list) & (sp2.Chr %in% chrom_cluster_list))) %>%
      dplyr::mutate(clust = cluster_num)
    final_df_with_groups <- rbind(final_df_with_groups,temp)
  }
  # convert to character for discrete values coloring :
  final_df_with_groups$clust <- as.character(final_df_with_groups$clust)
  # build a second dataframe with dots not on linkage groups to display in grey :
  non_linkage_df <- final_df_with_groups %>% dplyr::select(-clust)
  # check that the clust column doesn't exist :
  if ("clust" %in% colnames(mbh_df_to_return)) {
    mbh_df_to_return <- mbh_df_to_return %>%
      dplyr::select(-clust)
  }
  non_linkage_df <- dplyr::setdiff(mbh_df_to_return,non_linkage_df) %>%
    dplyr::mutate(clust = "0")
  final_df_with_groups <- rbind(final_df_with_groups,non_linkage_df)
  mbh_df_to_return <- merge(mbh_df_to_return,final_df_with_groups)
  
  return(mbh_df_to_return)
}