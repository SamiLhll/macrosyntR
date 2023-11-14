# plot_chord_diagram
#
# This is a function to plot the orthology as chord diagrams
#' @title plot the Macro-synteny as a chord diagram
#' @description This is a function to plot the chord diagrams to compare the macro synteny of two or more species. 
#' It requires as input an orthologs_df loaded by load_orthologs()
#'
#' @param orthologs_df dataframe. orthologs with genomic coordinates loaded by the load_orthologs()
#' @param species_labels list of characters. names of the species to display on the plot
#' @param species_labels_size integer. size of the labels (default = 2)
#' @param color_by string. name of the column in the orthologs_df to color the links by (default = "sp1.Chr")
#' @param custom_color_palette list of characters. palette to use for the coloring of the links following the argument color_by
#' @param reorder_chromosomes logical. (default = TRUE) tells whether to reorder the chromosomes in clusters as implemented in reorder_macrosynteny()
#' @param remove_non_linkage_orthologs logical. (default = TRUE) tells wether to remove the orthologs that are not within significant linkage groups as calculated by comput_linkage_groups().
#' @param species_labels_hpos (default =-400)
#' @param label_size integer. size of the labels to display on the ideograms (default = 2)
#' @param ideogram_fill character. name of the colors to fill the ideograms with (default = "white")
#' @param ideogram_color character. name of the colors to draw the borders of the ideograms with (default = "black")
#' @param ideogram_height integer. height of the ideograms (default = 4)
#' @param ribbons_curvature float. curvature of the ribbons (default = 0.1)
#' @param ribbons_alpha float. alpha of the ribbons (default = 0.5)
#' 
#' @seealso [load_orthologs()]
#' @seealso [reorder_macrosynteny()]
#' @seealso [compute_linkage_groups()]
#' 
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @importFrom rlang sym :=
#' @importFrom dplyr arrange group_by mutate row_number ungroup
#' 
#' @examples 
#' # basic usage of plot_oxford_grid : 
#' 
#' orthologs_table <- system.file("extdata","my_orthologs.tab",package="macrosyntR")
#' 
#' my_orthologs <- read.table(orthologs_table,header=TRUE)
#'
#' plot_chord_diagram(my_orthologs,species_labels = c("B. flo","P. ech"))
#' 
#' @export

plot_chord_diagram <- function(orthologs_df,
                               species_labels = NULL,
                               species_labels_size = 5,
                               color_by = "sp1.Chr",
                               custom_color_palette = NULL,
                               reorder_chromosomes = TRUE,
                               remove_non_linkage_orthologs = TRUE,
                               species_labels_hpos = -400,
                               label_size = 2, 
                               ideogram_fill = "white",
                               ideogram_color = "black",
                               ideogram_height = 4,
                               ribbons_curvature = 0.1,
                               ribbons_alpha = 0.5) {
  
  ### 1 - Compute the amount of species (max 60 species) and reorder Chromosomes :
  
  amount_of_species <- length(which(colnames(orthologs_df) %in% paste0("sp",c(1:60),".Chr")))
  
  if (reorder_chromosomes) {
    orthologs_df_to_plot <- reorder_multiple_macrosyntenies(orthologs_df)
  }
  else { orthologs_df_to_plot <- orthologs_df}
  
  if (remove_non_linkage_orthologs){
    orthologs_df_to_plot <- subset_linkage_orthologs(orthologs_df_to_plot)
  }
  
  ### 2 - Initialize plot :
  
  p <- ggplot() + ggplot2::theme_void()
  
  ################################################
  ### Start for loop to compute ideogram and link coordinates 
  
  for (i in c(1:amount_of_species)) {
    
    ### 3 - Compute Ideograms coordinates : 
    species <- paste0("sp",i)
    
    # Compute X and Y coordinates for the ideograms of the species :
    species_ideograms_coordinates <- get_ideograms_coordinates(orthologs_df_to_plot,
                                                               species_number = i,
                                                               amount_of_species = amount_of_species,
                                                               ideogram_height = ideogram_height)
    
    
    ### 4 - Compute links coordinates and plot :
    
    if (i+1 <= amount_of_species) {
      departure_species <- species
      arrival_species <- paste0("sp",i+1)
      # subset orthologs_df to have only departure and arrival species
      
      # merge with departure species :
      links_df <- merge(orthologs_df_to_plot,species_ideograms_coordinates, by = paste0(departure_species,".Chr"))
      # merge with arrival species :
      arrival_species_coordinates <- get_ideograms_coordinates(orthologs_df_to_plot,
                                                               species_number = i + 1,
                                                               amount_of_species = amount_of_species,
                                                               ideogram_height = ideogram_height)
      links_df <- merge(links_df,arrival_species_coordinates,by = paste0(arrival_species,".Chr"))
      
      # recompute indexes by adding the id_starts to it :
      Chr_column <- paste0(departure_species,".Chr")
      Start_column <- paste0(departure_species,".Start")
      Index_column <- paste0(departure_species,".Index")
      links_df <- links_df %>% group_by(!!rlang::sym(Chr_column)) %>% arrange(!!rlang::sym(Index_column),.by_group = TRUE) %>%
        mutate( "{departure_species}.Index" := dplyr::row_number()) %>% ungroup()
      
      links_df[[paste0(departure_species,".Index")]] <- links_df[[paste0(departure_species,".Index")]] + links_df[[paste0(departure_species,".id_starts")]]
      
      Chr_column <- paste0(arrival_species,".Chr")
      Start_column <- paste0(arrival_species,".Start")
      Index_column <- paste0(arrival_species,".Index")
      links_df <- links_df %>% group_by(!!rlang::sym(Chr_column)) %>% arrange(!!rlang::sym(Index_column),.by_group = TRUE) %>%
        mutate( "{arrival_species}.Index" := dplyr::row_number()) %>% ungroup()
      
      links_df[[paste0(arrival_species,".Index")]] <- links_df[[paste0(arrival_species,".Index")]] + links_df[[paste0(arrival_species,".id_starts")]]
      # add X_center :
      links_df[[paste0(departure_species,"_",arrival_species,".Xcenter")]] <- (links_df[[paste0(departure_species,".Index")]] + links_df[[paste0(arrival_species,".Index")]]) / 2
      # add Y_center :
      links_df[[paste0(departure_species,"_",arrival_species,".Ycenter")]] <- (links_df[[paste0(departure_species,".id_y")]] + links_df[[paste0(arrival_species,".id_y")]]) / 2
      
      ## shuffle links order :
      random_rows_order <- sample(nrow(links_df))
      links_df <- links_df[random_rows_order,]
      
      # PLOT links : 
      color_by_sym <- ggplot2::ensym(color_by)
      p <- p + 
        # First half of the links : 
        geom_curve(data = links_df,
                   aes(x = .data[[paste0(departure_species,".Index")]],
                       xend = .data[[paste0(departure_species,"_",arrival_species,".Xcenter")]],
                       y = .data[[paste0(departure_species,".id_ymin")]],
                       yend = .data[[paste0(departure_species,"_",arrival_species,".Ycenter")]],
                       color = !!color_by_sym),
                   curvature = -ribbons_curvature,
                   linewidth = .5, alpha = ribbons_alpha) +
        # Second half of the links :
        geom_curve(data = links_df,
                   aes(x = .data[[paste0(departure_species,"_",arrival_species,".Xcenter")]],
                       xend = .data[[paste0(arrival_species,".Index")]],
                       y = .data[[paste0(departure_species,"_",arrival_species,".Ycenter")]],
                       yend = .data[[paste0(arrival_species,".id_ymax")]],
                       color = !!color_by_sym),
                   curvature = ribbons_curvature,
                   linewidth = .5, alpha = ribbons_alpha)
    } # END IF
    
    ### Plot ideograms :
    # Plot ideograms :
    p <- p + 
      # draw ideograms :
      annotate("rect", 
               xmin = species_ideograms_coordinates[[paste0(species,".id_starts")]],
               xmax = species_ideograms_coordinates[[paste0(species,".id_stops")]],
               ymin = species_ideograms_coordinates[[paste0(species,".id_ymin")]], 
               ymax = species_ideograms_coordinates[[paste0(species,".id_ymax")]],
               alpha = 1,fill = ideogram_fill, color = ideogram_color) +
      # add chromosome labels :
      geom_point(data = species_ideograms_coordinates,
                 aes(x=.data[[paste0(species,".id_mids")]],
                     y=.data[[paste0(species,".id_y")]]),
                 size = label_size + 1,color = "white",shape = 15) +
      annotate(geom="text",
               x=species_ideograms_coordinates[[paste0(species,".id_mids")]],
               y=species_ideograms_coordinates[[paste0(species,".id_y")]],
               label=species_ideograms_coordinates[[paste0(species,".Chr")]],
               size = label_size,angle = 0)
    # add species labels :
    if (! is.null(species_labels)) {
      p <- p + 
        annotate(geom="text",
                 x=species_labels_hpos,
                 y=min(species_ideograms_coordinates[[paste0(species,".id_y")]]),
                 label=species_labels[i],
                 size = species_labels_size,angle = 0,
                 fontface = 'italic')
      
    }
  } # END_FOR
  
  
  if (is.null(custom_color_palette)) {
    Ref_palette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                     "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                     "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
                     "#8A7C64", "#599861")
    custom_color_palette <- Ref_palette[1:length(unique(links_df[[color_by]]))]
  }
  
  
  p <- p + scale_color_manual(values = custom_color_palette)
  
  
  return(p)
  
}