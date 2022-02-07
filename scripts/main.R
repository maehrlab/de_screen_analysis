# ========================================================================================
# ===================  Change this path to point to this folder. =========================
# ========================================================================================
proj_dir = "PATH/TO/THIS/FOLDER"
# ========================================================================================

setwd(proj_dir)
renv::activate()
library(cowplot)
library(magrittr)
library(freezr)
library(Seurat)
library(thymusatlastools2)
library(ggplot2)
library(Matrix)
flash_freeze = freezr::configure_flash_freeze( project_directory = proj_dir, 
                                              setup_scripts = 
                                                c( "tools/guide_handling_tools.Rmd", 
                                                   "tools/gAmp_merging_tools.Rmd",
                                                   "tools/diff_expr_tools.Rmd",
                                                   "tools/heatmapping_tools.Rmd",
                                                   "tools/guides_by_cluster_tools.Rmd" ))
REPLICATE_COLORS = c("Terarep1" = "goldenrod3", "Terarep2" = "navy", 
                     "he_perturb_rep1" = "goldenrod3", "he_perturb_rep2" = "navy")
HE_IDENT_COLORS = c("0" = "cadetblue4", "1" = "orange", "2" = "cadetblue2", 
                    "A" = "cadetblue4", "C" = "orange", "B" = "cadetblue2")
DE_IDENT_COLORS = c("0" = "goldenrod", "1" = "deeppink", "2" = "cadetblue", "3" = "brown")
UMAP_AXES = c("UMAP1", "UMAP2")


## ----------------------------------------------------------------------------------------------------------------------------------
inventory_make( inv_location = file.path(proj_dir, "results"))
getwd()
print(proj_dir)
flash_freeze()


## ----------------------------------------------------------------------------------------------------------------------------------

# Explore raw data and then integrate guides and do data cleaning.
flash_freeze( "clustering/setup.Rmd" )
flash_freeze( "clustering/integrate_gAmp.Rmd" ) 
flash_freeze( "clustering/recluster_no_dub.Rmd" ) 

# Re-cluster, characterize clusters, and study guide effects.
flash_freeze( "clustering/recluster_no_SOX17.Rmd" ) 
flash_freeze( "clustering/guides_by_cluster.Rmd" )
flash_freeze( "clustering/cluster_markers.Rmd" ) 
flash_freeze( "quantseq/quantseq.Rmd" ) 

# Reanalyze main cluster
flash_freeze( "C0/C0_recluster.Rmd" )
flash_freeze( "C0/C0_cluster_markers.Rmd" ) 

# Export data for use in MIMOSCA, then grab the results once it's done.

# Side effect: unlike the rest of my scripts, which only save to the freezr destination, 
# mimosca_export.Rmd script sends data to mimosca/data/ .
flash_freeze( "guide_effects/mimosca_export.Rmd" ) 
# Here, you need to move to mimosca/ and run `python tera_de_c0.py` . You'll need numpy and pandas.
# Then move the results you want to use into a folder and point the import script to that.
flash_freeze( "guide_effects/mimosca_import.Rmd" ) 

# Characterize guides' effects
flash_freeze( "guide_effects/guide_effects.Rmd") 
flash_freeze( "guide_effects/guide_effects_assemble.Rmd") 
flash_freeze( "guide_effects/guide_effects_plot.Rmd") 
flash_freeze( "guide_effects/guide_effects_plot_abundance.Rmd") 

# Analyze follow-up: a FOXA2 knockdown in a further-along lineage, hepatic endoderm.
flash_freeze( "hepatic/HE_set_up_data.Rmd" ) 
flash_freeze( "hepatic/HE_cluster.Rmd" ) 
flash_freeze( "hepatic/HE_plot.Rmd" )  
flash_freeze( "hepatic/HE_markers.Rmd" ) 


