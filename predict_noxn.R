# Make predictions from random forest
library(ranger)
library(raster)

# extent bounding boxes of a regular grid of 60 tiles covering the extent of covariates
EXTENT_DF <- read.csv("C:/Users/ginge/Documents/Python/nci_ndr/extent_df.csv")

# aligned covariate rasters
COVAR_PATH_LIST <- list(
  'precip_variability'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/wc2.0_bio_5m_15.tif",
  'population'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/inhabitants_avg_1990_2015.tif",
  'pigs'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/5_Pg_2010_Da.tif",
  'cattle'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/5_Ct_2010_Da.tif",
  'flash_flow'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/mean_div_range_1990_2015.tif",
  'average_flow'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/average_flow_1990_2015.tif",
  'percent_no_sanitation'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/no_sanitation_provision_avg_2000-2015.tif",
  'proportion_urban'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/perc_urban_5min.tif",
  'sand_percent'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/SNDPPT_M_sl1_10km_ll.tif",
  'clay_percent'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/CLYPPT_M_sl1_10km_ll.tif",
  'depth_to_groundwater'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/gwt_cm_sav_level12.tif"
  )

# aligned n export rasters for a set of scenarios
N_EXPORT_PATH_LIST = list(
  'extensification_bmps_irrigated'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_extensification_bmps_irrigated_300.0_D8_export_reduce28sum_md5_f28fd3b8e444044165bf8356f16da2cb.tif',
  'extensification_bmps_rainfed'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_extensification_bmps_rainfed_300.0_D8_export_reduce28sum_md5_6bf450d1e48ce523a96ef39971d14dff.tif',
  'extensification_current_practices'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_extensification_current_practices_300.0_D8_export_reduce28sum_md5_a7f97c7ea643732dcb9a37a38261d9be.tif',
  'extensification_intensified_irrigated'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_extensification_intensified_irrigated_300.0_D8_export_reduce28sum_md5_91b8c02c962a3ded85d7437a04544518.tif',
  'extensification_intensified_rainfed'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_extensification_intensified_rainfed_300.0_D8_export_reduce28sum_md5_1160ef8f4fb9b51a92c45fa7c081da97.tif',
  'fixedarea_bmps_irrigated'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_fixedarea_bmps_irrigated_300.0_D8_export_reduce28sum_md5_ce57595bd7f1aa0ee0b5aba91c5e5417.tif',
  'fixedarea_bmps_rainfed'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_fixedarea_bmps_rainfed_300.0_D8_export_reduce28sum_md5_b7d9eda12a76ceb6c5ed72501304318f.tif',
  'fixedarea_intensified_irrigated'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_fixedarea_intensified_irrigated_300.0_D8_export_reduce28sum_md5_5ed83d7dd073df9dd2ff5e55f493d1c9.tif',
  'fixedarea_intensified_rainfed'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_fixedarea_intensified_rainfed_300.0_D8_export_reduce28sum_md5_41c8043f4399419fcc995117be7c2a8f.tif',
  'grazing_expansion'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_grazing_expansion_300.0_D8_export_reduce28sum_md5_ebc5c0c91bbabff356fa1a0c57bfc01c.tif',
  'restoration'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_restoration_300.0_D8_export_reduce28sum_md5_0f7080ca2c1cf5b8ee10d7b39e3cce85.tif',
  'sustainable_currentpractices'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/compressed_sustainable_currentpractices_300.0_D8_export_reduce28sum_md5_11cbca90d8c151b07336aba2e917ec05.tif',
  'baseline_currentpractices'='C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/baseline_currentpractices_300.0_D8_export.tif'
)

# Use random foresets model to predict se from a raster stack
predict_se <- function(covar_stack, model, save_as) {
  out <- raster(covar_stack, layer=0)
  bs <- blockSize(out)
  out <- writeStart(out, save_as, overwrite=TRUE)
  for (i in 1:bs$n) {
    covar_df <- as.data.frame(getValues(covar_stack, row=bs$row[i], nrows=bs$nrows[i]))
    na_vals <- apply(covar_df, 1, function(x) sum(is.na(x)))
    p <- numeric(length=nrow(covar_df))
    p[na_vals > 0] <- NA
    if (length(p[na_vals==0]) > 0) {
      p[na_vals==0] <- raster::predict(
        model, covar_df[na_vals==0, ], type='se')$se
    }
    out <- writeValues(out, p, bs$row[i])
  }
  out <- writeStop(out)
  return(out)
}

# Predict response from a raster stack
predict_response <- function(covar_stack, model, save_as) {
  out <- raster(covar_stack, layer=0)
  bs <- blockSize(out)
  out <- writeStart(out, save_as, overwrite=TRUE)
  for (i in 1:bs$n) {
    covar_df <- as.data.frame(getValues(covar_stack, row=bs$row[i], nrows=bs$nrows[i]))
    na_vals <- apply(covar_df, 1, function(x) sum(is.na(x)))
    p <- numeric(length=nrow(covar_df))
    p[na_vals > 0] <- NA
    if (length(p[na_vals==0]) > 0) {
      p[na_vals==0] <- raster::predict(
        model, covar_df[na_vals==0, ], type='response')$predictions
    }
    out <- writeValues(out, p, bs$row[i])
  }
  out <- writeStop(out)
  return(out)
}

# Make response and se predictions by tiling
# Parameters:
#   - intermediate_dir (string): directory where intermediate tiled results
#     should be writting
#   - EXTENT_DF (data frame): data frame containing bounding boxes for tiles
#     covering the globe
#   - covariate_name_list (named list): named list giving covariates. The order
#     of covariates in this list must match the order of predictors that were
#     used to train the random forests model
#   - covariate_path_list (named list): list indexed by covariate name and
#     giving paths to aligned covariate rasters for the globe
#   - ranger_model (model object): trained ranger model that should be used to
#     make predictions
#   - response_target_path (string): path to location on disk where noxn
#     global raster should be saved
#   - se_target_path (string): path to location on disk where standard error
#     global raster should be saved
tile_and_predict <- function(
    intermediate_dir, EXTENT_DF, covariate_name_list, covariate_path_list,
    ranger_model, response_target_path, se_target_path) {
  native_format_list <- list()
  covar_tile_dir <- paste(intermediate_dir, 'covar_tile', sep='/')
  dir.create(covar_tile_dir, showWarnings=FALSE)
  pred_tile_dir <- paste(intermediate_dir, 'pred_tile', sep='/')
  dir.create(pred_tile_dir, showWarnings=FALSE)
  merge_list_se <- list()
  merge_list_response <- list()
  for (r in 1:NROW(EXTENT_DF)) {
    extent <- extent(as.numeric(EXTENT_DF[r, ]))
    for(cv in covariate_name_list){
      save_as <- paste(covar_tile_dir, paste(
        tools::file_path_sans_ext(
          basename(covariate_path_list[[cv]])), '.grd', sep=''),
        sep='/')
      crop(raster(covariate_path_list[[cv]]), y=extent,
           filename=save_as, overwrite=TRUE)
      native_format_list[[cv]] <- save_as
    }
    covar_stack <- stack(lapply(native_format_list, FUN=raster))
    names(covar_stack) <- covariate_name_list

    # predict response
    response_pred_tile_filename <- paste(
      pred_tile_dir, paste("response_", r, ".tif", sep=''), sep='/')
    response_pred <- predict_response(
      covar_stack, model=ranger_model, save_as=response_pred_tile_filename)
    merge_list_response[[r]] <- response_pred_tile_filename

    # predict standard error of the response
    se_pred_tile_filename <- paste(
      pred_tile_dir, paste("se_pred_", r, ".tif", sep=''), sep='/')
    se_pred <- predict_se(
      covar_stack, model=ranger_model, save_as=se_pred_tile_filename)
    merge_list_se[[r]] <- se_pred_tile_filename
  }
  # merge the tiles together
  merged_response <- Reduce(merge, lapply(merge_list_response, raster))
  writeRaster(merged_response, filename=response_target_path)
  merged_se <- Reduce(merge, lapply(merge_list_se, raster))
  writeRaster(merged_se, filename=se_target_path)

  # clean up
  unlink(covar_tile_dir, recursive=TRUE)
  unlink(pred_tile_dir, recursive=TRUE)
}

# surface noxn
# noxn and covariate observations for surface water
NOXN_PREDICTOR_SURF_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_surf_2000_2015.csv"

# train the random forest model
train_df = read.csv(NOXN_PREDICTOR_SURF_DF_PATH)
ranger_train_df <- train_df[complete.cases(train_df), ]
ranger_model <- ranger(noxn~., keep.inbag=TRUE, data=ranger_train_df,
                       mtry=2, min.node.size=5, importance='impurity')  # using mtry and min.node.size tuned by caret

output_dir <- "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/R_ranger_pred"
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
intermediate_dir <- "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/R_ranger_pred/surf_intermediate"
dir.create(intermediate_dir)
covariate_df <- subset(train_df, select=-c(noxn))
covariate_name_list <- colnames(covariate_df)
for (scenario in names(N_EXPORT_PATH_LIST)) {
  COVAR_PATH_LIST[['n_export']] = N_EXPORT_PATH_LIST[[scenario]]
  response_target_path <- paste(output_dir, paste(
    "surface_noxn_", scenario, ".tif", sep=''), sep='/')
  se_target_path <- paste(output_dir, paste(
    "surface_noxn_se_", scenario, ".tif", sep=''), sep='/')
  tile_and_predict(
    intermediate_dir, EXTENT_DF, covariate_name_list, COVAR_PATH_LIST,
    ranger_model, response_target_path, se_target_path)
}

# ground noxn
# GEMStat, China, USGS, and Ouedraogo (Africa) combined
NOXN_PREDICTOR_GR_DF_PATH = "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/noxn_predictor_df_gr_Ouedraogo_USGS_China_GEMStat.csv"

# train the random forest model
train_df <- read.csv(NOXN_PREDICTOR_GR_DF_PATH)
ranger_train_df <- train_df[complete.cases(train_df), ]
ranger_ground_model <- ranger(
  noxn~., keep.inbag=TRUE, data=ranger_train_df, mtry=2, min.node.size=5,
  importance='impurity')  # using mtry and min.node.size tuned by caret

output_dir <- "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/Ouedraogo_USGS_China_GEMStat/R_ranger_pred"
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
intermediate_dir <- "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/Ouedraogo_USGS_China_GEMStat/R_ranger_pred/ground_intermediate"
dir.create(intermediate_dir, recursive=TRUE, showWarnings=FALSE)
covariate_df <- subset(train_df, select=-c(noxn))
covariate_name_list <- colnames(covariate_df)
for (scenario in names(N_EXPORT_PATH_LIST)) {  # c("fixedarea_currentpractices")) {
  COVAR_PATH_LIST[['n_export']] = N_EXPORT_PATH_LIST[[scenario]]
  response_target_path <- paste(output_dir, paste(
    "ground_noxn_", scenario, ".tif", sep=''), sep='/')
  se_target_path <- paste(output_dir, paste(
    "ground_noxn_se_", scenario, ".tif", sep=''), sep='/')
  if(!file.exists(response_target_path)){
    tile_and_predict(
      intermediate_dir, EXTENT_DF, covariate_name_list, COVAR_PATH_LIST,
      ranger_ground_model, response_target_path, se_target_path)
  }
}
