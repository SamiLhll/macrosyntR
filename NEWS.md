# macrosyntR 0.3.3

### Enhancements : 

* removed dependency on ggthemes package that is scheduled for archival.

# macrosyntR 0.3.1

### New features :


* Two species and more can now be visualized on a chord diagram using the function plot_chord_diagram()


### Enhancements :

* Redundant chromosome names from one species to another are now accepted   

* load_orthologs now handles more than two species. It uses the argument 'bedfiles' instead of sp1_bed and sp2_bed which are still working when used together. 

# macrosyntR 0.2.21

### Enhancements :

* added control when chromosome names are redundant in the 'orthologs_df' provided to 'reorder_macrosynteny()'.
It now raises an error with explanations about what's wrong.

# macrosyntR 0.2.19

### New features :

* created the function get_syntenic_genes(). It takes an orthologs_df as input and
outputs a table with the details of all detected blocks of two or more consecutive genes.

* created the function reverse_species_order(). It takes an orthologs_df as input and
outputs it with sp1 changed in sp2 and the other way around. It can be called in plot_oxford_grid().

* Added the argument keep_sp1_raw_order in reorder_macrosynteny(). When set to TRUE, it returns an orthologs_df with only the sp2.Chr reordered, and doesn't change the order of sp1.Chr compared with the input data. It can be also be called in plot_oxford_grid().

* 'plot_oxford_grid()' now features an option to (dis)able the coloring of orthologs depending on if they
occur on significant linkage groups or not. By default, when setting a color_by argument, the orthologs that are located on non-significant linkage groups are displayed in grey. It is possible to disable this behavior by calling the function with setting argument *shade_non_significant* to *TRUE*.


### Bug fixes :

* Corrected a bug that occurred when loading bed files with more than 3 fields in 'load_orthologs()'.
This function now handles bed files that have fields after the 4th column (seqName).   

* Corrected a bug that happened when trying to set a custom color palette through the color_palette argument in 'plot_oxford_grid()'. It is now possible to set a custom color palette as a list of color names.


### Enhancements :

* Added a `NEWS.md` file to track changes to the package.   
* Added documentation about how to customize the plots
