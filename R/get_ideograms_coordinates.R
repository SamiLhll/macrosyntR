# get_ideograms_coordinates (internal)
#
# This is an internal function that is called by plot_chord_diagram()
#' @title get the ideograms coordinates to plot chord diagrams (internal)
#' @description This is an internal function to compute the ideograms coordinates to plot chord diagrams.
#'
#' @param orthologs_df dataframe. orthologs with genomic coordinates loaded with load_orthologs()
#' @param species_number integer. the number of the species being processed. (default = 1)
#' @param amount_of_species integer. the total amount of species in the orthologs_df (default = 2)
#' @param ideogram_height integer. the height of the ideograms (default = 10)
#' @param gap_size integer. the size of the gaps between ideograms (default = 40)
#' 
#' 
#' @return A dataframe object
#'
#' @import tidyr
#' @importFrom rlang sym
#' @importFrom dplyr mutate select
#' 
#' 
#' @noRd

get_ideograms_coordinates <- function(orthologs_df,
                                      species_number = 1,
                                      amount_of_species = 2,
                                      ideogram_height = 10,
                                      gap_size = 40) {
  
  id_y <- starts <- stops <- NULL
  
  ### 1 - compute values based on number of orthologs per chromosome :
  species <- paste0("sp",species_number)
  Chr_column <- paste0(species,".Chr")
  Index_column <- paste0(species,".Index")
  ideograms_to_process <- orthologs_df %>%
    group_by(!!rlang::sym(Chr_column)) %>%
    arrange(!!rlang::sym(Index_column),.by_group = TRUE) %>%
    summarise(n = n()) %>% ungroup()
  # print(ideograms_to_process)
  ideograms_to_process$stops <- cumsum(ideograms_to_process$n) + 1
  ideograms_to_process$starts <- c(0,ideograms_to_process$stops[1:length(ideograms_to_process$stops) - 1])
  # print(ideograms_to_process)
  ### 2 - Compute new starts and stops with adding a gap between the chromosomes :
  genome_gap <- 0
  temp_ideograms <- NULL
  
  for (i in ideograms_to_process[[Chr_column]]) {
    temp_ideogram <- ideograms_to_process[ideograms_to_process[[Chr_column]] == i,]
    temp_ideogram <- temp_ideogram %>%
      mutate(starts = replace(starts, TRUE, starts + genome_gap),
             stops = replace(stops, TRUE, stops + genome_gap))
    temp_ideograms <- rbind(temp_ideograms,temp_ideogram)
    genome_gap <- genome_gap + gap_size
    
  }
  ### 3 - Compute mids that will be the x aes in the plot :
  ideograms_to_return <- temp_ideograms %>%
    mutate(mids = (starts + stops) / 2) %>%
    select(-n)
  
  ### 4 - Compute Y coordinates and add to the dataframe :
  ideograms_to_return <- ideograms_to_return %>%
    mutate(id_y = 100 - ((species_number -1) * (100 / amount_of_species)))%>%
    mutate(id_ymin = id_y - (ideogram_height / 2),
           id_ymax = id_y + (ideogram_height / 2))
  
  ### 5 - rename the columns with species number :
  colnames_to_return <- c(Chr_column,paste0(species,c(".id_stops",".id_starts",".id_mids",".id_y",".id_ymin",".id_ymax")))
  colnames(ideograms_to_return) <- colnames_to_return
  
  return(ideograms_to_return)
}