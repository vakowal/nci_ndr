# Make predictions from random forest
library(ranger)
library(raster)

# extent bounding boxes of a regular grid of 60 tiles covering the extent of covariates
EXTENT_DF <- read.csv("C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/extent_df.csv")

# aligned covariate rasters
COVAR_PATH_LIST <- list(
  'precip_variability'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/wc2.0_bio_5m_15.tif",
  'population'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/inhabitants_avg_1990_2015.tif",
  'pigs'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/5_Pg_2010_Da.tif",
  'cattle'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/5_Ct_2010_Da.tif",
  'flash_flow'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/mean_div_range_1990_2015.tif",
  'average_flow'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/average_flow_1990_2015.tif",
  'percent_no_sanitation'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/no_sanitation_provision_avg_2000-2015.tif",
  'proportion_urban'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/perc_urban_5min.tif"
)

# aligned n export rasters for a set of scenarios
N_EXPORT_PATH_LIST = list(
  'baseline'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/sum_aggregate_to_0.084100_n_export_baseline_napp_rate_global_md5_b210146a5156422041eb7128c147512f.tif",
  'ag_expansion'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/sum_aggregate_to_0.084100_n_export_ag_expansion_global_md5_ea15fb82df52d49a1d0c4ffe197cdd0d.tif",
  'ag_intensification'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/sum_aggregate_to_0.084100_n_export_ag_intensification_global_md5_2734116e8c452f4c484ebcb574aab665.tif",
  'restoration'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_ground/sum_aggregate_to_0.084100_n_export_restoration_napp_rate_global_md5_7f9ddf313e414a68cbb8ba204101b190.tif"
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
    repsonse_pred_tile_filename <- paste(
      pred_tile_dir, paste("response_", r, ".tif", sep=''), sep='/')
    response_pred <- predict_response(
      covar_stack, model=ranger_model, save_as=repsonse_pred_tile_filename)
    merge_list_response[[r]] <- repsonse_pred_tile_filename

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
                       mtry=2, min.node.size=5)  # using mtry and min.node.size tuned by caret

output_dir <- "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/R_ranger_pred"
dir.create(output_dir, showWarnings=FALSE)
intermediate_dir <- "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate"
covariate_df <- subset(train_df, select=-c(noxn))
covariate_name_list <- colnames(covariate_df)
for (scenario in names(N_EXPORT_PATH_LIST)) {
  COVAR_PATH_LIST[['n_export']] = N_EXPORT_PATH_LIST[[scenario]]
  response_target_path <- paste(output_dir, paste(
    "surface_noxn_", scenario, ".tif", sep=''), sep='/')
  response_target_path <- paste(output_dir, paste(
    "surface_noxn_se_", scenario, ".tif", sep=''), sep='/')
  tile_and_predict(
    intermediate_dir, EXTENT_DF, covariate_name_list, COVAR_PATH_LIST,
    ranger_model, response_target_path, se_target_path)
}

# ground noxn
# noxn and covariate observations for ground water
NOXN_PREDICTOR_GR_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_gr.csv"