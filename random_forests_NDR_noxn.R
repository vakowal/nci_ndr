# random forests model: NOxN ~ NDR for World Bank

# Train and evaluate a random forests model.
#
#  Parameters:
#   noxn_obs_path (string): path to csv file containing nitrate observations at
#     points. Must contain a column 'OBJECTID'
#   covar_df_path (string): path to csv file containing predictor values at
#     points. Must contain a column 'OBJECTID' that matches values in
#     `noxn_obs_path`.
#   stn_covar_cols (vector of strings): column names in `covar_df_path` that
#     should be used as predictors
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
    noxn_obs_path, covar_df_path, stn_covar_cols, out_dir) {
  library(caret)
  noxn_obs <- read.csv(noxn_obs_path)
  noxn_obs <- noxn_obs[, c('OBJECTID', 'noxn')]
  covar_df <- read.csv(covar_df_path)
  # covar_df[covar_df$average_flow == 0, 'flash_flow'] <- 0
  stn_covar_df <- covar_df[, stn_covar_cols]
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

# match table: objectid+GEMS station number, all surface stations with noxn observations
OBJECTID_MATCH_CSV_SURF <- "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_objectid_stnid_match_table.csv"
STN_OBJ_MATCH_DF_SURF <- read.csv(OBJECTID_MATCH_CSV_SURF)

# surface stations with noxn observations
SURFACE_NOXN_STATION_CSV <- "F:/NCI_NDR/Data worldbank/station_data/WB_surface_stations_noxn_obs.csv"

# dates by which to restrict surface NOxN observations
MIN_DATE = "2000-01-01"  # "1995-01-01"
MAX_DATE = "2015-12-31"  # 2010-12-31"

# NOxN observations, all stations
NOXN_OBSERVATIONS_CSV <- "F:/NCI_NDR/Data worldbank/station_data/noxn_obs_all.csv"
NOXN_OBS <- read.csv(NOXN_OBSERVATIONS_CSV, stringsAsFactors=FALSE)
NOXN_OBS$Sample.Date <- as.Date(NOXN_OBS$Sample.Date, format="%Y-%m-%d")

# surface covariate table: 5 min pixel extent (not snapped)
# PIXEL_ORIG_COVAR_CSV <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/WB_station_orig_5min_pixel/combined_covariates.csv"
# N application predictor
# PIXEL_ORIG_COVAR_CSV <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/N_application_predictor/WB_station_orig_5min_pixel/combined_covariates.csv"
# 30 s scale
PIXEL_ORIG_COVAR_CSV <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/WB_station_orig_30s_pixel/combined_covariates.csv"

# station metadata
STATION_METADATA_CSV <- "F:/NCI_NDR/Data worldbank/station_data/station_metadata.csv"
STN_DF <- read.csv(STATION_METADATA_CSV)
stn_cols <- colnames(STN_DF)[c(1, 4, 5, 8:10, 15:16)]
stn_subset <- STN_DF[, stn_cols]
stn_subset[stn_subset$Water.Type == 'Groundwater station', 'ground_v_surface'] <- 'groundwater'
stn_subset[stn_subset$Water.Type == 'Lake station', 'ground_v_surface'] <- 'surface'
stn_subset[stn_subset$Water.Type == 'River station', 'ground_v_surface'] <- 'surface'
NOXN_BY_STN <- merge(NOXN_OBS, stn_subset, all.x=TRUE)
NOXN_BY_STN$Sample.Date <- as.Date(NOXN_BY_STN$Sample.Date, format="%Y-%m-%d")

# surface stations only: define filtered list of stations
# restrict to surface stations
surface_df_stn_list <- unique(NOXN_BY_STN[NOXN_BY_STN$ground_v_surface == 'surface',
                                               'GEMS.Station.Number'])

# aggregate noxn observations: one mean observation per station per year, inside
# the selected time period
surface_stn <- NOXN_BY_STN[NOXN_BY_STN$ground_v_surface == 'surface', ]
noxn_obs_subset <- surface_stn[(surface_stn$Sample.Date >= MIN_DATE) &
                               (surface_stn$Sample.Date <= MAX_DATE), ]
noxn_obs_subset$year <- format(noxn_obs_subset$Sample.Date, format="%Y")
mean_noxn_by_year <- aggregate(Value~GEMS.Station.Number+year+Unit,
                               data=noxn_obs_subset, FUN=mean)
noxn_obs_restr <- mean_noxn_by_year[, c(1, 4)]
colnames(noxn_obs_restr)[2] <- 'noxn'

###### process covariate data: 5min pixel level (not snapped) ######
covar_df <- read.csv(PIXEL_ORIG_COVAR_CSV)
covar_df[covar_df$average_flow == 0, 'flash_flow'] <- 0
match_df <- read.csv(OBJECTID_MATCH_CSV_SURF)
covar_df <- merge(covar_df, match_df)
STN_DF <- read.csv(SURFACE_NOXN_STATION_CSV)
stn_covar_df <- merge(covar_df, STN_DF)

stn_covar_cols <- c('GEMS.Station.Number', 'n_export', 'precip_variability',
                    'population', 'pigs', 'cattle', 'flash_flow', 'average_flow',
                    'percent_no_sanitation', 'proportion_urban')  # pixel-level covariates
stn_covar_df <- stn_covar_df[, stn_covar_cols]

# merge covariates with NOxN observations
combined_df <- merge(stn_covar_df, noxn_obs_restr, by="GEMS.Station.Number")  # , all=TRUE)
rf_covar_df <- combined_df[, 2:dim(combined_df)[2]]
# save filtered data frame for analysis in Python
# write.csv(rf_covar_df, "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_surf.csv", row.names=FALSE)

# random forests model, surface stations
# out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/surface/subset_2000_2015/WB_station_orig_5min_pixel"
# N application predictor
# out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/N_app_predictor/surface/subset_2000_2015/WB_station_orig_5min_pixel"
# 30 sec scale
out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/surface/subset_2000_2015/WB_station_orig_30s_pixel"
dir.create(out_dir, recursive=TRUE)

# K-fold cross-validation
set.seed(491)
fitControl <- trainControl(method='repeatedcv', number=10, repeats=10)

# fit the random forests model
ranger_rf <- train(noxn ~ ., data=rf_covar_df,
                   method='ranger', trControl=fitControl,
                   na.action=na.omit, importance='impurity')
sink(paste(out_dir, "ranger_rf_summary.txt", sep='/'))
print(ranger_rf)
sink()
rf_var_imp <- varImp(ranger_rf)
sink(paste(out_dir, "ranger_rf_var_importance.txt", sep='/'))
print(rf_var_imp)
sink()

# Bar chart: variable importance (for an example model)
var_imp_csv <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_3.2.20/ground/all_stations/R_ranger_pred/var_imp.csv"
var_imp <- read.csv(var_imp_csv)
var_imp$Predictor <- reorder(var_imp$variable, var_imp$importance, descending=TRUE)
# var_imp_res <- var_imp[var_imp$Overall > 7, ]
library(ggplot2)
p <- ggplot(var_imp, aes(x=Predictor, y=importance))
p <- p + geom_bar(stat='identity') + coord_flip()
p <- p + ylab("Variable Importance") + xlab("")
pngname <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_3.2.20/ground/all_stations/R_ranger_pred/var_imp_plot.png"
png(file=pngname, units="in", res=300, width=3.5, height=4)
print(p)
dev.off()

##### random forests for groundwater #####
GROUNDWATER_PREDICTORS <- c(
  'OBJECTID', 'n_export', 'precip_variability',
  'population', 'depth_to_groundwater', 'cattle',
  'average_flow', 'percent_no_sanitation', 'clay_percent',
  'flash_flow', 'sand_percent', 'pigs', 'proportion_urban')

# covariates collected at GEMStat, USGS, China, and Ouedraogo sites combined
noxn_obs_path <- "F:/NCI_NDR/Data groundwater merged/combined_noxn_groundwater_GEMStat_Gu_USGS_Ouedraogo.csv"
# N application as predictor: covar_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/N_application_predictor/groundwater_GEMStat_Gu_USGS_Ouedraogo/combined_covariates.csv")
# covar_df_path <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/groundwater_GEMStat_Gu_USGS_Ouedraogo/combined_covariates.csv"
# out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/ground/GEMStat_China_USGS_Ouedraogo"
# out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/N_app_predictor/ground/GEMStat_China_USGS_Ouedraogo"
# 30 s scale
covar_df_path <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/groundwater_GEMStat_Gu_USGS_Ouedraogo_30s_pixel/combined_covariates.csv"
out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/ground/GEMStat_China_USGS_Ouedraogo_30s"
train_evaluate_rf(
  noxn_obs_path, covar_df_path, GROUNDWATER_PREDICTORS, out_dir)
