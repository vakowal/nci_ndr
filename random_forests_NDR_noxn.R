# random forests model: NOxN ~ NDR for World Bank
library(caret)

# match table: objectid+GEMS station number, all surface stations with noxn observations
OBJECTID_MATCH_CSV_SURF <- "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_objectid_stnid_match_table.csv"
STN_OBJ_MATCH_DF_SURF <- read.csv(OBJECTID_MATCH_CSV_SURF)

# surface stations with noxn observations
SURFACE_NOXN_STATION_CSV <- "F:/NCI_NDR/Data worldbank/station_data/WB_surface_stations_noxn_obs.csv"

# stations with change in N fertilizer application rates, 1980 to 2013
DELTA_FERT_CSV_SURF <- "F:/NCI_NDR/Data fertilizer Lu Tian/percent_change_1980_2013_surf.csv"
DELTA_FERT_DF_SURF <- read.csv(DELTA_FERT_CSV_SURF)

# match table: objectid+GEMS station number, groundwater stations with noxn observations
OBJECTID_MATCH_CSV_GR <- "F:/NCI_NDR/Data worldbank/station_data/WB_groundwater_stations_noxn_obs_objectid_match.csv"
STN_OBJ_MATCH_DF_GR <- read.csv(OBJECTID_MATCH_CSV_GR)

# ground stations with noxn observations
GR_NOXN_STATION_CSV <- "F:/NCI_NDR/Data worldbank/station_data/WB_ground_stations_noxn_obs.csv"

DELTA_FERTILIZER_1980_2013_CSV_GR <- "F:/NCI_NDR/Data fertilizer Lu Tian/perc_change_mean_N_application_1980_2013_point_ground.csv"
DELTA_FERT_DF_GR <- read.csv(DELTA_FERTILIZER_1980_2013_CSV_GR)

# dates by which to restrict NOxN observations
MIN_DATE = "2000-01-01"  # "1995-01-01"
MAX_DATE = "2015-12-31"  # 2010-12-31"

# NOxN observations, all stations
NOXN_OBSERVATIONS_CSV <- "F:/NCI_NDR/Data worldbank/station_data/noxn_obs_all.csv"
NOXN_OBS <- read.csv(NOXN_OBSERVATIONS_CSV, stringsAsFactors=FALSE)
NOXN_OBS$Sample.Date <- as.Date(NOXN_OBS$Sample.Date, format="%Y-%m-%d")

# covariate table: 5 min pixel extent (not snapped)
PIXEL_ORIG_COVAR_CSV <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/WB_station_orig_5min_pixel/combined_covariates.csv"

# covariate table: 5 min pixel, groundwater
PIXEL_ORIG_COVAR_GR_CSV <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/WB_groundwater_station_orig_5min_pixel/combined_covariates.csv"

library(ggplot2)

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

# throwaway
# calculate annual 95% perc value
noxn_copy <- NOXN_BY_STN
noxn_copy$year <- format(noxn_copy$Sample.Date, format="%Y")
noxn_95_perc <- aggregate(Value~GEMS.Station.Number+year+Unit+ground_v_surface,
                          data=noxn_copy, FUN=quantile, probs=0.95)
colnames(noxn_95_perc)[5] <- 'noxn_95_perc'
exceeded_95_perc <- noxn_95_perc[noxn_95_perc$noxn_95_perc >= 5.6, ]

noxn_yearly <- aggregate(Value~GEMS.Station.Number+year+Unit+ground_v_surface,
                               data=noxn_copy, FUN=mean)
colnames(noxn_yearly)[5] <- 'noxn_mean'
atest <- merge(noxn_95_perc, noxn_yearly)

exceeded_yearly <- noxn_yearly[noxn_yearly$Value >= 5.6, ]

# save the raw data in the relevant subset
subset_to_save <- NOXN_BY_STN[, c('Sample.Date', 'Value', 'Country.Name', 'ground_v_surface')]
colnames(subset_to_save)[2] <- 'noxn_mg/L'
write.csv(subset_to_save, paste(outdir, 'noxn_observations_full.csv', sep='/'), row.names=FALSE)

# calculate and save average yearly noxn for all stations, all years
copy_df <- NOXN_BY_STN
copy_df$year <- format(copy_df$Sample.Date, format="%Y")
avg_annual_noxn_by_stn <- aggregate(Value~GEMS.Station.Number+year+ground_v_surface+Latitude+Longitude,
                                    data=copy_df, FUN=mean)
colnames(avg_annual_noxn_by_stn)[6] <- 'noxn_mg/L'
write.csv(avg_annual_noxn_by_stn, "F:/NCI_NDR/Data worldbank/avg_annual_noxn_by_stn.csv",
          row.names=FALSE)

# surface stations only: define filtered list of stations
# restrict to surface stations
surface_df_stn_list <- unique(NOXN_BY_STN[NOXN_BY_STN$ground_v_surface == 'surface',
                                               'GEMS.Station.Number'])
# restrict by stability of N fertilizer application rates
delta_fert_by_stn <- merge(DELTA_FERT_DF_SURF, STN_OBJ_MATCH_DF_SURF)
fert_app_subset_stn_list <- delta_fert_by_stn[
  delta_fert_by_stn$percent_change <= 171.7, 'GEMS.Station.Number']

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
# save data frame identifying stations with observations in modeling period
noxn_obs_surf_stations_in_modeling_period <- as.data.frame(
  noxn_obs_restr[!duplicated(noxn_obs_restr$GEMS.Station.Number), 'GEMS.Station.Number'])
colnames(noxn_obs_surf_stations_in_modeling_period)[1] <- 'GEMS.Station.Number'
noxn_obs_surf_stations_in_modeling_period$model_per <- 1
noxn_obs_surf_stations_in_modeling_period$fert_app <- 0
for (r in 1:NROW(noxn_obs_surf_stations_in_modeling_period)) {
  if (noxn_obs_surf_stations_in_modeling_period[r, 'GEMS.Station.Number'] %in%
      fert_app_subset_stn_list) {
    noxn_obs_surf_stations_in_modeling_period[r, 'fert_app'] <- 1
  }
}
write.csv(noxn_obs_surf_stations_in_modeling_period,
          "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/training_set_station_shapefiles/surface_stn_modeling_subset.csv", row.names=FALSE)

###### process covariate data: 5min pixel level (not snapped) ######
covar_df <- read.csv(PIXEL_ORIG_COVAR_CSV)
covar_df[covar_df$average_flow == 0, 'flash_flow'] <- 0
match_df <- read.csv(OBJECTID_MATCH_CSV_SURF)
covar_df <- merge(covar_df, match_df)
STN_DF <- read.csv(SURFACE_NOXN_STATION_CSV)
stn_covar_df <- merge(covar_df, STN_DF)

# add dummy variables for each water body type
# stn_covar_df$lake <- 0
# stn_covar_df[stn_covar_df$Water.Type == 'Lake station', 'lake'] <- 1
# stn_covar_df$river <- 0
# stn_covar_df[stn_covar_df$Water.Type == 'River station', 'river'] <- 1
stn_covar_cols <- c('GEMS.Station.Number', 'n_export', 'precip_variability',
                    'population', 'pigs', 'cattle', 'flash_flow', 'average_flow',
                    'percent_no_sanitation', 'proportion_urban')  # pixel-level covariates
stn_covar_df <- stn_covar_df[, stn_covar_cols]

# merge covariates with NOxN observations
combined_df <- merge(stn_covar_df, noxn_obs_restr, by="GEMS.Station.Number")  # , all=TRUE)

# save records not filtered by fertilizer application trends
non_filtered_df <- combined_df[, 2:dim(combined_df)[2]]
write.csv(non_filtered_df, "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_surf_2000_2015.csv", row.names=FALSE)

# restrict by stability of N fert application
combined_df_restr <- combined_df[combined_df$GEMS.Station.Number %in% fert_app_subset_stn_list, ]
rf_covar_df <- combined_df_restr[, 2:dim(combined_df_restr)[2]]
# save filtered data frame for analysis in Python
write.csv(rf_covar_df, "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_surf.csv", row.names=FALSE)

# random forests model, surface stations
out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_3.2.20/surface/subset_2000_2015/WB_station_orig_5min_pixel"
dir.create(out_dir, recursive=TRUE)

# K-fold cross-validation
set.seed(491)
fitControl <- trainControl(method='repeatedcv', number=10, repeats=10)
# rf_covar_df$lake <- as.factor(rf_covar_df$lake)
# rf_covar_df$river <- as.factor(rf_covar_df$river)

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
# define filtered list of stations
# groundwater stations
ground_df_stn_list <- unique(
    NOXN_BY_STN[NOXN_BY_STN$ground_v_surface == 'groundwater',
    'GEMS.Station.Number'])
# restrict by date
time_subset_stn_list <- unique(
  NOXN_OBS[(NOXN_OBS$Sample.Date >= MIN_DATE) &
  (NOXN_OBS$Sample.Date <= MAX_DATE), 'GEMS.Station.Number'])
# restrict by stability of N fertilizer application rates
# delta_fert_by_stn <- merge(DELTA_FERT_DF_GR, STN_OBJ_MATCH_DF_GR)
# fert_app_subset_stn_list <- delta_fert_by_stn[
#   delta_fert_by_stn$perc_change_N_application_1980_2013 <=
#   GROUND_FERT_MEDIAN_VAL, 'GEMS.Stati']
# subset_stn_list_time <- intersect(ground_df_stn_list, time_subset_stn_list)
# subset_stn_list <- intersect(subset_stn_list_time, fert_app_subset_stn_list)

# aggregate noxn observations: one mean observation per station per year
gr_stn <- NOXN_BY_STN[NOXN_BY_STN$ground_v_surface == 'groundwater', ]
gr_stn$year <- format(gr_stn$Sample.Date, format="%Y")
mean_noxn_by_year <- aggregate(Value~GEMS.Station.Number+year+Unit,
                               data=gr_stn, FUN=mean)
noxn_obs_restr <- mean_noxn_by_year[, c(1, 4)]
colnames(noxn_obs_restr)[2] <- 'noxn'

# process covariate data
covar_df <- read.csv(PIXEL_ORIG_COVAR_GR_CSV)
covar_df[covar_df$average_flow == 0, 'flash_flow'] <- 0
match_df <- read.csv(OBJECTID_MATCH_CSV_GR)
covar_df <- merge(covar_df, match_df)
STN_DF <- read.csv(GR_NOXN_STATION_CSV)
stn_covar_df <- merge(covar_df, STN_DF, by.x='GEMS.Stati', by.y='GEMS.Station.Number')
stn_covar_cols <- c('GEMS.Stati', 'n_export', 'precip_variability',
                    'population', 'depth_to_groundwater', 'cattle',
                    'average_flow', 'percent_no_sanitation', 'clay_percent',
                    'flash_flow', 'sand_percent', 'pigs', 'proportion_urban')
stn_covar_df <- stn_covar_df[, stn_covar_cols]

# merge groundwater covariates with noxn obs
combined_df <- merge(stn_covar_df, noxn_obs_restr, by.x='GEMS.Stati', by.y='GEMS.Station.Number')

# write out groundwater raw data for random forests analysis in Python
combined_df_restr <- combined_df[, 2:dim(combined_df)[2]]
write.csv(combined_df_restr, "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_gr.csv", row.names=FALSE)

# restrict by time period and n fert application
# combined_df_restr <- combined_df[
#   combined_df$GEMS.Stati %in% subset_stn_list, ]
# # OR
# # restrict by n fert application
# combined_df_restr <- combined_df[
#   combined_df$GEMS.Stati %in% fert_app_subset_stn_list, ]
# # OR
# # restrict by time period
# combined_df_restr <- combined_df[
#   combined_df$GEMS.Stati %in% subset_stn_list_time, ]

# random forests model
library(caret)
out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_3.2.20/ground/all_stations/WB_station_orig_5min_pixel"
dir.create(out_dir, recursive=TRUE)
rf_covar_df <- combined_df_restr

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

### random forests model: groundwater data from India
noxn_obs <- read.csv("F:/NCI_NDR/Data worldbank/India_groundwater/wells_noxn_obs.csv")
noxn_obs <- noxn_obs[, c('stncode', 'nitrateN_mean')]
match_df <- read.csv("F:/NCI_NDR/Data worldbank/India_groundwater/wells_match_table.csv")
covar_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/India_groundwater_5min_pixel/combined_covariates.csv")

covar_df[covar_df$average_flow == 0, 'flash_flow'] <- 0
stn_covar_df <- merge(covar_df, match_df)
stn_covar_cols <- c('stncode', 'n_export', 'precip_variability',
                    'population', 'depth_to_groundwater', 'cattle',
                    'average_flow', 'percent_no_sanitation', 'clay_percent',
                    'flash_flow', 'sand_percent', 'pigs', 'proportion_urban')
stn_covar_df <- stn_covar_df[, stn_covar_cols]
india_obs <- merge(stn_covar_df, noxn_obs, by='stncode')
colnames(india_obs)[dim(india_obs)[2]] <- 'noxn'

# random forests model
library(caret)
out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_3.2.20/ground/India_wells"
dir.create(out_dir, recursive=TRUE)
india_rf_df <- india_obs[, 2:dim(india_obs)[2]]  # drop stncode column

# K-fold cross-validation
set.seed(491)
fitControl <- trainControl(method='repeatedcv', number=10, repeats=10)
# fit the random forests model
ranger_rf <- train(noxn ~ ., data=india_rf_df,
                   method='ranger', trControl=fitControl,
                   na.action=na.omit, importance='impurity')
sink(paste(out_dir, "ranger_rf_summary.txt", sep='/'))
print(ranger_rf)
sink()
rf_var_imp <- varImp(ranger_rf)
sink(paste(out_dir, "ranger_rf_var_importance.txt", sep='/'))
print(rf_var_imp)
sink()

### random forests model: groundwater data from China
noxn_obs <- read.csv("F:/NCI_NDR/Data groundwater Gu et al/noxn_observations.csv")
covar_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/China_groundwater_5min_pixel/combined_covariates.csv")
covar_df[covar_df$average_flow == 0, 'flash_flow'] <- 0
stn_covar_cols <- c('OBJECTID', 'n_export', 'precip_variability',
                    'population', 'depth_to_groundwater', 'cattle',
                    'average_flow', 'percent_no_sanitation', 'clay_percent',
                    'flash_flow', 'sand_percent', 'pigs', 'proportion_urban')
stn_covar_df <- covar_df[, stn_covar_cols]
china_obs <- merge(stn_covar_df, noxn_obs, by='OBJECTID')
colnames(china_obs)[14] <- 'noxn'

# random forests model
library(caret)
out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_3.2.20/ground/China_wells"
dir.create(out_dir, recursive=TRUE)
china_rf_df <- china_obs[, 2:dim(china_obs)[2]]  # drop OBJECTID column

# K-fold cross-validation
set.seed(491)
fitControl <- trainControl(method='repeatedcv', number=10, repeats=10)
# fit the random forests model
ranger_rf <- train(noxn ~ ., data=china_rf_df,
                   method='ranger', trControl=fitControl,
                   na.action=na.omit, importance='impurity')
sink(paste(out_dir, "ranger_rf_summary.txt", sep='/'))
print(ranger_rf)
sink()
rf_var_imp <- varImp(ranger_rf)
sink(paste(out_dir, "ranger_rf_var_importance.txt", sep='/'))
print(rf_var_imp)
sink()

# combine Gemstat, India, and China data frames into a single groundwater df
gemstat_gr_df <- read.csv("C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_gr.csv")
groundwater_df <- rbind(gemstat_gr_df, india_pixel_rf_df, china_rf_df)  # or india_rf_df
write.csv(groundwater_df, "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/noxn_predictor_df_gr_IndiaPixel_China_GEMStat.csv",
          row.names=FALSE)

out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_3.2.20/ground/GEMStat_IndiaPixel_China_combined"
dir.create(out_dir, recursive=TRUE)

# K-fold cross-validation
set.seed(491)
fitControl <- trainControl(method='repeatedcv', number=10, repeats=10)
# fit the random forests model
ranger_rf <- train(noxn ~ ., data=groundwater_df,
                   method='ranger', trControl=fitControl,
                   na.action=na.omit, importance='impurity')
sink(paste(out_dir, "ranger_rf_summary.txt", sep='/'))
print(ranger_rf)
sink()
rf_var_imp <- varImp(ranger_rf)
sink(paste(out_dir, "ranger_rf_var_importance.txt", sep='/'))
print(rf_var_imp)
sink()

# combine Gemstat and China only
groundwater_df <- rbind(gemstat_gr_df, china_rf_df)
write.csv(groundwater_df, "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/noxn_predictor_df_gr_China_GEMStat.csv",
          row.names=FALSE)
out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_3.2.20/ground/GEMStat_China_combined"
dir.create(out_dir, recursive=TRUE)

# K-fold cross-validation
set.seed(491)
fitControl <- trainControl(method='repeatedcv', number=10, repeats=10)
# fit the random forests model
ranger_rf <- train(noxn ~ ., data=groundwater_df,
                   method='ranger', trControl=fitControl,
                   na.action=na.omit, importance='impurity')
sink(paste(out_dir, "ranger_rf_summary.txt", sep='/'))
print(ranger_rf)
sink()
rf_var_imp <- varImp(ranger_rf)
sink(paste(out_dir, "ranger_rf_var_importance.txt", sep='/'))
print(rf_var_imp)
sink()

# India data aggregated within 5min pixels
noxn_obs <- read.csv("F:/NCI_NDR/Data worldbank/India_groundwater/mean_noxn_by_5min_pixel.csv")
noxn_obs <- noxn_obs[, c('OBJECTID', 'nitrateN_mean')]
covar_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/India_groundwater_pixel_centroid_5min_pixel/combined_covariates.csv")

covar_df[covar_df$average_flow == 0, 'flash_flow'] <- 0
stn_covar_cols <- c('OBJECTID', 'n_export', 'precip_variability',
                    'population', 'depth_to_groundwater', 'cattle',
                    'average_flow', 'percent_no_sanitation', 'clay_percent',
                    'flash_flow', 'sand_percent', 'pigs', 'proportion_urban')
stn_covar_df <- covar_df[, stn_covar_cols]
india_obs <- merge(stn_covar_df, noxn_obs, by='OBJECTID')
colnames(india_obs)[dim(india_obs)[2]] <- 'noxn'

# random forests model
library(caret)
out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_3.2.20/ground/India_wells_aggregated_5min_pixel"
dir.create(out_dir, recursive=TRUE)
india_pixel_rf_df <- india_obs[, 2:dim(india_obs)[2]]  # drop OBJECTID column

# K-fold cross-validation
set.seed(491)
fitControl <- trainControl(method='repeatedcv', number=10, repeats=10)
# fit the random forests model
ranger_rf <- train(noxn ~ ., data=india_pixel_rf_df,
                   method='ranger', trControl=fitControl,
                   na.action=na.omit, importance='impurity')
sink(paste(out_dir, "ranger_rf_summary.txt", sep='/'))
print(ranger_rf)
sink()
rf_var_imp <- varImp(ranger_rf)
sink(paste(out_dir, "ranger_rf_var_importance.txt", sep='/'))
print(rf_var_imp)
sink()

# go to line 335 to combine India pixel-level db with GEMStat and China

# groundwater data for the US from the decadalchange in groundwater quality viewer
noxn_obs <- read.csv("F:/NCI_NDR/Data groundwater USGS/Decadal_change_in_groundwater_NO3N_datav2.csv")
noxn_obs <- noxn_obs[, c('OBJECTID', 'Mean_02_13')]
covar_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/USGS_groundwater_5min_pixel/combined_covariates.csv")

covar_df[covar_df$average_flow == 0, 'flash_flow'] <- 0
stn_covar_cols <- c('OBJECTID', 'n_export', 'precip_variability',
                    'population', 'depth_to_groundwater', 'cattle',
                    'average_flow', 'percent_no_sanitation', 'clay_percent',
                    'flash_flow', 'sand_percent', 'pigs', 'proportion_urban')
stn_covar_df <- covar_df[, stn_covar_cols]
usgs_obs <- merge(stn_covar_df, noxn_obs, by='OBJECTID')
colnames(usgs_obs)[dim(usgs_obs)[2]] <- 'noxn'
usgs_rf_df <- usgs_obs[, 2:dim(usgs_obs)[2]]  # drop OBJECTID column

# combine Gemstat, USGS, and China data frames into a single groundwater df
gemstat_gr_df <- read.csv("C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_gr.csv")
# generate china_rf_df with lines 301-316
groundwater_df <- rbind(gemstat_gr_df, china_rf_df, usgs_rf_df)
write.csv(groundwater_df, "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/noxn_predictor_df_gr_USGS_China_GEMStat.csv",
          row.names=FALSE)

out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_3.2.20/ground/GEMStat_China_USGS_combined"
dir.create(out_dir, recursive=TRUE)

# K-fold cross-validation
set.seed(491)
fitControl <- trainControl(method='repeatedcv', number=10, repeats=10)
# fit the random forests model
ranger_rf <- train(noxn ~ ., data=groundwater_df,
                   method='ranger', trControl=fitControl,
                   na.action=na.omit, importance='impurity')
sink(paste(out_dir, "ranger_rf_summary.txt", sep='/'))
print(ranger_rf)
sink()
rf_var_imp <- varImp(ranger_rf)
sink(paste(out_dir, "ranger_rf_var_importance.txt", sep='/'))
print(rf_var_imp)
sink()
