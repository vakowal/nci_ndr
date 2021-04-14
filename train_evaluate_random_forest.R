# Train and evaluate a random forests model.
#
#  Parameters:
#   noxn_obs_path (string): path to csv file containing nitrate observations at
#     points. Must contain a column 'OBJECTID'
#   covar_df_path (string): path to csv file containing predictor values at
#     points. Must contain a column 'OBJECTID' that matches values in
#     `noxn_obs_path`.
#   out_dir (string): path to directory where results should be stored.
#
#   Side effects:
#     creates out_dir and two files inside it:
#       ranger_rf_summary.txt, a summary of the random forest model performance
#       ranger_rf_var_importance.txt, a summary of variable importance values
#         in the trained model
#
#   Returns:
#     None
train_evaluate_rf <- function(
    noxn_obs_path, covar_df_path, out_dir) {
  library(caret)
  noxn_obs <- read.csv(noxn_obs_path)
  noxn_obs <- noxn_obs[, c('OBJECTID', 'noxn')]
  covar_df <- read.csv(covar_df_path)
  rf_df <- merge(covar_df, noxn_obs, by='OBJECTID')
  colnames(rf_df)[dim(rf_df)[2]] <- 'noxn'
  if ('X' %in% colnames(rf_df)) {
    rf_df <- subset(rf_df, select=-c(X, OBJECTID))
  }
  else {
    rf_df <- subset(rf_df, select=-c(OBJECTID))
  }

  dir.create(out_dir, recursive=TRUE)
  
  # K-fold cross-validation
  set.seed(491)
  fitControl <- trainControl(method='repeatedcv', number=10, repeats=10)
  # fit the random forests model
  ranger_rf <- train(noxn ~ ., data=rf_df,
                     method='ranger', trControl=fitControl,
                     na.action=na.omit, importance='impurity')
  sink(paste(out_dir, "ranger_rf_summary.txt", sep='/'))
  print(ranger_rf)
  sink()
  rf_var_imp <- varImp(ranger_rf)
  sink(paste(out_dir, "ranger_rf_var_importance.txt", sep='/'))
  print(rf_var_imp)
  sink()
}

# directory containing cloned repository
repo_dir <- "C:/Users/ginge/Documents/Python/nci_ndr"

# where results of the random forest model should be stored
results_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_3.18.21"

# testing 4 options for NDR predictor, 3.19.21
outer_dir <- 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/covar_df_3.19.21'
# surface
noxn_obs_path_surface <- paste(repo_dir, "/noxn_obs_surface.csv", sep='/')

covar_df_path <- paste(outer_dir, 'surface_WB_station_orig_5min_pixel_N_export_via_R', 'covar_df_surface.csv', sep='/')
out_dir <- paste(results_dir, "surface_WB_station_orig_5min_pixel_N_export_via_R", sep='/')
train_evaluate_rf(noxn_obs_path_surface, covar_df_path, out_dir)

# ground
noxn_obs_path <- paste(repo_dir, "/noxn_obs_GEMStat_Gu_USGS_Ouedraogo_ground.csv", sep='/')
outer_dir <- 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/covar_df_3.19.21'

covar_df_path <- paste(outer_dir, 'groundwater_GEMStat_Gu_USGS_Ouedraogo_5min_pixel_N_export_via_R', 'covar_df_ground.csv', sep='/')
out_dir <- paste(results_dir, "groundwater_GEMStat_Gu_USGS_Ouedraogo_5min_pixel_N_export_via_R", sep='/')
train_evaluate_rf(noxn_obs_path, covar_df_path, out_dir)

# one set of predictors
# surface water nitrate observations and predictors, 5 min resolution
noxn_obs_path_surface <- paste(repo_dir, "/noxn_obs_surface.csv", sep='/')
covar_df_path_5min <- paste(repo_dir, "/covar_df_5min_surface.csv", sep='/')
out_dir <- paste(results_dir, "/surface/subset_2000_2015/WB_station_orig_5min_pixel", sep='/')
train_evaluate_rf(noxn_obs_path_surface, covar_df_path_5min, out_dir)

# groundwater nitrate observations and predictors, 5 min resolution
noxn_obs_path <- paste(repo_dir, "/noxn_obs_GEMStat_Gu_USGS_Ouedraogo_ground.csv", sep='/')
covar_df_path_5min <- paste(repo_dir, "/covar_df_5min_ground.csv", sep='/')
out_dir <- paste(results_dir, "/ground/GEMStat_China_USGS_Ouedraogo_5min", sep='/')
train_evaluate_rf(noxn_obs_path, covar_df_path_5min, out_dir)

## add accumulated and accumulated-normalized predictors
# surface
# noxn_obs_path_surface <- paste(repo_dir, "/noxn_obs_surface.csv", sep='/')
# covar_df_path <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/accumulated/WB_station_orig_5min_pixel/covar_df_combined_with_local.csv"
# out_dir <- paste(results_dir, "/surface/subset_2000_2015/WB_station_orig_5min_pixel/accumulated", sep='/')
# train_evaluate_rf(noxn_obs_path_surface, covar_df_path, out_dir)
# 
# # ground
# noxn_obs_path <- paste(repo_dir, "/noxn_obs_GEMStat_Gu_USGS_Ouedraogo_ground.csv", sep='/')
# covar_df_path <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/accumulated/groundwater_GEMStat_Gu_USGS_Ouedraogo/covar_df_combined_with_local.csv"
# out_dir <- paste(results_dir, "/ground/GEMStat_China_USGS_Ouedraogo_5min/accumulated", sep='/')
# train_evaluate_rf(noxn_obs_path, covar_df_path, out_dir)
