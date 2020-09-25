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
  rf_df <- merge(stn_covar_df, noxn_obs, by='OBJECTID')
  colnames(rf_df)[dim(rf_df)[2]] <- 'noxn'
  rf_df <- rf_df[, 2:dim(rf_df)[2]]  # drop OBJECTID column
  
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

# surface water nitrate observations and predictors, 5 min resolution
noxn_obs_path_surface <- "~/noxn_obs_surface.csv"
covar_df_path_5min <- "~/covar_df_5min_surface.csv"
out_dir <- "~/surface/subset_2000_2015/WB_station_orig_5min_pixel"
train_evaluate_rf(noxn_obs_path, covar_df_path_30s, out_dir)

# surface water nitrate observations and predictors, 30 s resolution
noxn_obs_path_surface <- "~/noxn_obs_surface.csv"
covar_df_path_5min <- "~/covar_df_30s_surface.csv"
out_dir <- "~/surface/subset_2000_2015/WB_station_orig_30s_pixel"
train_evaluate_rf(noxn_obs_path, covar_df_path_30s, out_dir)

# groundwater nitrate observations and predictors, 5 min resolution
noxn_obs_path <- "~/noxn_obs_GEMStat_Gu_USGS_Ouedraogo_ground.csv"
covar_df_path_5min <- "~/covar_df_5min_ground.csv"
out_dir <- "~/ground/GEMStat_China_USGS_Ouedraogo_5min"
train_evaluate_rf(noxn_obs_path, covar_df_path_30s, out_dir)

# groundwater nitrate observations and predictors, 30 s resolution
noxn_obs_path <- "~/noxn_obs_GEMStat_Gu_USGS_Ouedraogo_ground.csv"
covar_df_path_30s <- "~/covar_df_30s_ground.csv"
out_dir <- "~/ground/GEMStat_China_USGS_Ouedraogo_30s"
train_evaluate_rf(noxn_obs_path, covar_df_path_30s, out_dir)