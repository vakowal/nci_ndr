# Make predictions from random forest
library(ranger)
library(raster)

# extent bounding boxes of a regular grid of 60 tiles covering the extent of covariates
EXTENT_DF <- read.csv("C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/globe_tiles/extent_df.csv")

# noxn and covariate observations for surface water
NOXN_PREDICTOR_SURF_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_surf_2000_2015.csv"

# train the random forest model
train_df = read.csv(NOXN_PREDICTOR_SURF_DF_PATH)
ranger_train_df <- train_df[complete.cases(train_df), ]
ranger_model <- ranger(noxn~., keep.inbag=TRUE, data=ranger_train_df,
                       mtry=2, min.node.size=5)  # using mtry and min.node.size tuned by caret

# surface covariate rasters
surface_covar_path_list <- list(
  'n_export'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/sum_aggregate_to_0.084100_n_export_baseline_napp_rate_global_md5_b210146a5156422041eb7128c147512f.tif",
  'precip_variability'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/wc2.0_bio_5m_15.tif",
  'population'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/inhabitants_avg_1990_2015.tif",
  'pigs'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/5_Pg_2010_Da.tif",
  'cattle'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/5_Ct_2010_Da.tif",
  'flash_flow'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/mean_div_range_1990_2015.tif",
  'average_flow'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/average_flow_1990_2015.tif",
  'percent_no_sanitation'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/no_sanitation_provision_avg_2000-2015.tif",
  'proportion_urban'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/perc_urban_5min.tif"
)

# Function to use ranger.predict onto a raster
predfun <- function(covar_stack, model, save_as){
  out <- raster(covar_stack, layer=0)
  bs <- blockSize(out)
  out <- writeStart(out, save_as, overwrite=TRUE)
  for (i in 1:bs$n) {
    covar_df <- as.data.frame(getValues(covar_stack, row=bs$row[i], nrows=bs$nrows[i]))
    na_vals <- apply(covar_df, 1, function(x) sum(is.na(x)))
    p <- numeric(length=nrow(covar_df))
    p[na_vals > 0] <- NA
    if (length(p[na_vals==0]) > 0) {
      p[na_vals==0] <- predict(model, covar_df[na_vals==0, ], type='se')$se
    }
    out <- writeValues(out, p, bs$row[i])
  }
  out <- writeStop(out)
  return(out)
}

# extract covariates to small tiles
native_format_list <- list()
covar_tile_dir <- "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/raster_format_crop"
dir.create(covar_tile_dir, showWarnings=FALSE)
pred_tile_dir <- "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/baseline_surface_noxn_se_tiles"
dir.create(pred_tile_dir, showWarnings=FALSE)
merge_list <- list()
for (r in 1:NROW(EXTENT_DF)) {
  extent <- extent(as.numeric(EXTENT_DF[r, ]))
  for(cv in names(surface_covar_path_list)){
    save_as <- paste(covar_tile_dir, paste(
      tools::file_path_sans_ext(basename(surface_covar_path_list[[cv]])), '.grd', sep=''),
      sep='/')
    crop(raster(surface_covar_path_list[[cv]]), y=extent,
         filename=save_as, overwrite=TRUE)
    native_format_list[[cv]] <- save_as
  }
  covar_stack <- stack(lapply(native_format_list, FUN=raster))
  names(covar_stack) <- names(surface_covar_path_list)
  
  pred_tile_filename <- paste(pred_tile_dir, paste("se_pred_", r, ".tif"), sep='/')
  pred <- predfun(covar_stack, model=ranger_model, save_as=pred_tile_filename)
 merge_list[[r]] <- pred_tile_filename
}
# merge the tiles together
merged_filename <- "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/baseline_surface_noxn_se.tif"
merged <- Reduce(merge, lapply(merge_list, raster))
writeRaster(merged, filename=merged_filename)

