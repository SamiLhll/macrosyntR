# reorder_macrosynteny
#
# This is a function to reorder an orthologs_df, that was generated with load_orthologs(). It retrieves communities using igraph::cluster_fast_greedy.
#' @title Reorder the mbh_df before plotting
#' @description This is a function to reorder an orthologs_df, that was generated with load_orthologs(). It retrieves communities using igraph::cluster_fast_greedy.
#'
#' @param orthologs_df dataframe. mutual best hits with genomic coordinates loaded with load_orthologs()
#' @param pvalue_threshold numeric. threshold for significancy. (default equals 0.001)
#' @param keep_only_significant logical. (default equals FALSE) tells if the non significant linkage groups should be removed. It drastically speeds up the computation when using one highly fragmented genome.
#' @param keep_sp1_raw_order logical. (default equals FALSE) tells if the reordering should be constrained on the species1 and change just the order of the species2
#' 
#' @return A dataframe object
#' @seealso [load_orthologs()]
#' @seealso [compute_macrosynteny()]
#'
#' @importFrom igraph graph_from_data_frame cluster_fast_greedy groups
#' @importFrom dplyr select rename arrange group_by summarise ungroup mutate setdiff desc
#' @importFrom stringr str_subset str_replace
#' 
#'@examples 
#' # basic usage of reorder_macrosynteny : 
#' 
#' orthologs_table <- system.file("extdata","my_orthologs.tab",package="macrosyntR")
#' 
#' my_orthologs <- read.table(orthologs_table,header=TRUE)
#' 
#' my_orthologs_reordered <- reorder_macrosynteny(my_orthologs)
#' 
#' @export


reorder_macrosynteny <- function(orthologs_df,
                                 pvalue_threshold = 0.001,
                                 keep_only_significant = FALSE,
                                 keep_sp1_raw_order = FALSE) {
  
  contingency_table <- significant_entries <- significant_entries_for_graph <- NULL
  significant_entries_for_graph <- significant_association_undirected_graph <- clusters <- cluster_df <- NULL
  sp1_amounts <- sp2_amounts <- NULL
  final_levels_sp1 <- final_levels_sp2 <- NULL
  significant <- pval <- orthologs <- sp1.Chr <- sp2.Chr <- amounts <- clust <- n <- weight <- NULL
  
  # Error check : 
  # Error check : proper format for arguments :
  if (!(is.numeric(pvalue_threshold) & length(pvalue_threshold) == 1)) {stop("Wrong format for 'pvalue_threshold' argument. Must be a single value of type numeric")}
  if (!(is.logical(keep_only_significant) & length(keep_only_significant) == 1)) { stop("Wrong format for argument 'keep_only_significant'. Must be of type logical")}
  if (!(is.logical(keep_sp1_raw_order) & length(keep_sp1_raw_order) == 1)) { stop("Wrong format for argument 'keep_sp1_raw_order'. Must be of type logical")}
  # Error check : proper formatting of macrosynt_df
  # Error check : format of orthologs_df 
  required_fields <- c("sp1.ID","sp1.Index","sp1.Chr","sp2.ID","sp2.Index","sp2.Chr")
  for (i in required_fields) {
    if (isFALSE(i %in% colnames(orthologs_df))) {
      stop("Missing fields in the provided 'orthologs_df'. All the following columns are required : sp1.ID,sp1.Index,sp1.Chr,sp2.ID,sp2.Index,sp2.Chr")
    }
  }
  # Error check : orthologs_df is empty
  if (length(orthologs_df$sp1.Chr) == 0) {stop("Table provided through the 'orthologs_df' argument is empty")}
  # Warning check : when number of chromosomes is too high
  if (length(unique(orthologs_df$sp1.Chr)) >= 300) { 
    warning(paste0("The first species in the 'orthologs_df' has ",length(unique(orthologs_df$sp1.Chr))," chromosomes. Computational time can be very long on fragmented genomes"))
  }
  if (length(unique(orthologs_df$sp2.Chr)) >= 300) { 
    warning(paste0("The second species in the 'orthologs_df' has ",length(unique(orthologs_df$sp2.Chr))," chromosomes. Computational time can be very long on fragmented genomes"))
  }
  
  ##### 1 - Build an undirected and unweighted graph of connected chromosomes (significant amount of orthologs)
  # rename_chromosomes to ensure uniqueness of names in graph :
  levels(orthologs_df$sp1.Chr) <- paste0("sp1_Chr_",levels(orthologs_df$sp1.Chr))
  levels(orthologs_df$sp2.Chr) <- paste0("sp2_Chr_",levels(orthologs_df$sp2.Chr))
  # These prefix are removed at the end of the function
  
  # Get a table with only significant association using compute_macrosynteny from this package :
  contingency_table <- compute_macrosynteny(orthologs_df,pvalue_threshold)
  if (keep_only_significant) {
    significant_entries <- subset(contingency_table,significant == "yes")
  }
  else {
    significant_entries <- contingency_table
  }
  
  
  # build an Undirected weighted graph between sp1.chr and sp2.chr containing all significant edges :
  significant_entries_for_graph <- significant_entries %>% dplyr::select(-significant,-pval,) %>% dplyr::rename(from = sp1.Chr,to = sp2.Chr, weight = orthologs)
  # compute clusters of sp1.Chr and sp2.Chr that are directly or indirectly connected in the graph :
  significant_association_undirected_graph <- igraph::graph_from_data_frame(significant_entries_for_graph, directed = F)
  clusters <- igraph::cluster_fast_greedy(significant_association_undirected_graph)
  
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
    dplyr::arrange(dplyr::desc(amounts))
  cluster_df$num_clust <- letters[1:length(cluster_df$amounts)]
  ##### DONE : Computed order of clusters
  
  ##### 3 - Compute order of chromosomes by recomputing the levels of the factor. sorting by descending cluster size then descencding chromosome size.
  # get a table with number of orthologs for each chromosome
  if (keep_sp1_raw_order) {
    final_levels_sp1 <- levels(orthologs_df$sp1.Chr)
    final_levels_sp2 <- NULL
    sp1_amounts <- significant_entries %>% dplyr::group_by(sp1.Chr) %>% dplyr::summarise(n=sum(orthologs)) %>% dplyr::ungroup()
    sp2_amounts <- significant_entries %>% dplyr::group_by(sp2.Chr) %>% dplyr::summarise(n=sum(orthologs)) %>% dplyr::ungroup()
    chrom_clusters_reordered <- cluster_df$clust
    for (i in final_levels_sp1) {
      matching_clusters <- stringr::str_subset(chrom_clusters_reordered,paste0(i,","))
      chrom_cluster_list <- strsplit(matching_clusters,",")[[1]]
      for (chrom in chrom_cluster_list) {
        if (!(chrom %in% final_levels_sp1)) {
          if (!(chrom %in% final_levels_sp2)) {
            final_levels_sp2 <- c(final_levels_sp2,chrom)
          }
        }
      }
    }
  }
  else {
    sp1_amounts <- significant_entries %>% dplyr::group_by(sp1.Chr) %>% dplyr::summarise(n=sum(orthologs)) %>% dplyr::ungroup()
    sp2_amounts <- significant_entries %>% dplyr::group_by(sp2.Chr) %>% dplyr::summarise(n=sum(orthologs)) %>% dplyr::ungroup()
    chrom_clusters_reordered <- cluster_df$clust
    for(chrom_cluster in chrom_clusters_reordered) {
      chrom_cluster_list <- strsplit(chrom_cluster,",")[[1]]
      cluster_levels_sp1 <- NULL
      cluster_levels_sp2 <- NULL
      for (chrom in chrom_cluster_list) {
        if (chrom %in% orthologs_df$sp1.Chr) {
          cluster_levels_sp1 <- c(cluster_levels_sp1,chrom)
        }
        else {
          cluster_levels_sp2 <- c(cluster_levels_sp2,chrom)
        }
      }
      # order the chromosomes from larger to smaller (in amount of orthologs) :
      cluster_levels_sp1_df <- subset(sp1_amounts,sp1.Chr %in% cluster_levels_sp1) %>%
        dplyr::arrange(dplyr::desc(n))
      final_levels_sp1 <- c(final_levels_sp1,as.character(cluster_levels_sp1_df$sp1.Chr))
      cluster_levels_sp2_df <- subset(sp2_amounts,sp2.Chr %in% cluster_levels_sp2) %>%
        dplyr::arrange(dplyr::desc(n))
      final_levels_sp2 <- c(final_levels_sp2,as.character(cluster_levels_sp2_df$sp2.Chr))
      
    }}
  orthologs_df_to_return <- subset(orthologs_df, (sp1.Chr %in% final_levels_sp1) & (sp2.Chr %in% final_levels_sp2))
  orthologs_df_to_return$sp1.Chr <- factor(orthologs_df_to_return$sp1.Chr,levels = final_levels_sp1)
  orthologs_df_to_return$sp2.Chr <- factor(orthologs_df_to_return$sp2.Chr,levels = final_levels_sp2)
  ##### DONE : Computed order of chromosomes
  
  ##### 4 - Add an additional column with clusters number :
  final_df_with_groups <- NULL
  cluster_num <- 0
  # cluster_df calculated in part 2 of this function
  chrom_clusters_reordered <- cluster_df$clust
  for (i in chrom_clusters_reordered){
    chrom_cluster_list <- strsplit(i,",")[[1]]
    cluster_num <- cluster_num + 1
    temp <- subset(orthologs_df_to_return,((sp1.Chr %in% chrom_cluster_list) & (sp2.Chr %in% chrom_cluster_list))) %>%
      dplyr::mutate(clust = letters[cluster_num])
    final_df_with_groups <- rbind(final_df_with_groups,temp)
  }
  # convert to character for discrete values coloring :
  # final_df_with_groups$clust <- as.character(final_df_with_groups$clust)
  # build a second dataframe with dots not on linkage groups to display in grey :
  non_linkage_df <- final_df_with_groups %>% dplyr::select(-clust)
  # check that the clust column doesn't exist :
  if ("clust" %in% colnames(orthologs_df_to_return)) {
    orthologs_df_to_return <- orthologs_df_to_return %>%
      dplyr::select(-clust)
  }
  non_linkage_df <- dplyr::setdiff(orthologs_df_to_return,non_linkage_df) %>%
    dplyr::mutate(clust = "0")
  final_df_with_groups <- rbind(final_df_with_groups,non_linkage_df)
  orthologs_df_to_return <- merge(orthologs_df_to_return,final_df_with_groups)
  
  # Remove prefix that granted uniqueness of chromosome names :
  levels(orthologs_df_to_return$sp1.Chr) <- stringr::str_replace(levels(orthologs_df_to_return$sp1.Chr),"sp1_Chr_","")
  levels(orthologs_df_to_return$sp2.Chr) <- stringr::str_replace(levels(orthologs_df_to_return$sp2.Chr),"sp2_Chr_","")
  
  return(orthologs_df_to_return)
}