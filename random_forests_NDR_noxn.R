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

# compare mean nitrate in surface vs groundwater, by country, within time period 2000-2015
noxn_time_subset <- NOXN_BY_STN[(NOXN_BY_STN$Sample.Date >= MIN_DATE) &
                                  (NOXN_BY_STN$Sample.Date <= MAX_DATE), 
                                c("Value", "Country.Name", "ground_v_surface")]
mean_by_country_source <- aggregate(Value~Country.Name+ground_v_surface,
                                    data=noxn_time_subset, FUN=mean)
colnames(mean_by_country_source)[3] <- "noxn_mean"
sd_by_country_source <- aggregate(Value~Country.Name+ground_v_surface,
                                  data=noxn_time_subset, FUN=sd)
colnames(sd_by_country_source)[3] <- "noxn_sd"
median_by_country_source <- aggregate(Value~Country.Name+ground_v_surface,
                                      data=noxn_time_subset, FUN=median)
colnames(median_by_country_source)[3] <- "noxn_median"
summary_df <- merge(mean_by_country_source, sd_by_country_source, all=TRUE)
summary_df <- merge(summary_df, median_by_country_source, all=TRUE)
summary_df_reshape <- reshape(summary_df, idvar='Country.Name', timevar='ground_v_surface',
                              direction='wide')
summary_df_reshape$diff_mean <- summary_df_reshape$noxn_mean.surface - summary_df_reshape$noxn_mean.groundwater
summary_df_reshape$diff_median <- summary_df_reshape$noxn_median.surface - summary_df_reshape$noxn_median.groundwater
summary_df_reshape <- summary_df_reshape[
  c('Country.Name', 'diff_median', 'diff_mean', 'noxn_median.surface', 'noxn_median.groundwater',
    'noxn_mean.surface', 'noxn_mean.groundwater', 'noxn_sd.surface', 'noxn_sd.groundwater')
]
outdir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/ground_v_surface_observed_summary"
write.csv(summary_df_reshape, paste(outdir, "noxn_mean_median_by_country_ground_v_surface_2000-2015.csv", sep="/"),
          row.names=FALSE)

# another way: ANOVA
NOXN_BY_STN$year <- as.factor(format(NOXN_BY_STN$Sample.Date, format="%Y"))
NOXN_BY_STN$Country.Name <- as.factor(NOXN_BY_STN$Country.Name)
NOXN_BY_STN$ground_v_surface <- as.factor(NOXN_BY_STN$ground_v_surface)
model <- lm(Value~Country.Name+year+ground_v_surface, data=NOXN_BY_STN)
sink(paste(outdir, "noxn_anova_summary.txt", sep='/'))
summary(model)
sink()

# save the raw data in the relevant subset
subset_to_save <- NOXN_BY_STN[, c('Sample.Date', 'Value', 'Country.Name', 'ground_v_surface')]
colnames(subset_to_save)[2] <- 'noxn_mg/L'
write.csv(subset_to_save, paste(outdir, 'noxn_observations_full.csv', sep='/'), row.names=FALSE)

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
stn_covar_cols <- colnames(stn_covar_df)[c(1, 3:5, 8:11)]  # pixel-level covariates
stn_covar_df <- stn_covar_df[, stn_covar_cols]

# merge covariates with NOxN observations
combined_df <- merge(stn_covar_df, noxn_obs_restr, by="GEMS.Station.Number")  # , all=TRUE)

# restrict by stability of N fert application
combined_df_restr <- combined_df[combined_df$GEMS.Station.Number %in% fert_app_subset_stn_list, ]
rf_covar_df <- combined_df_restr[, 2:9]
# save filtered data frame for analysis in Python
write.csv(rf_covar_df, "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_surf.csv", row.names=FALSE)

# random forests model, surface stations
out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR/surface/N_application_subset_2000_2015/WB_station_orig_5min_pixel"

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
var_imp_csv <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/N_application_subset_2000_2015/WB_station_orig_5min_pixel/no_climate_zones/ranger_rf_var_importance.csv"
var_imp <- read.csv(var_imp_csv)
var_imp$Predictor <- reorder(var_imp$X, var_imp$Overall, descending=TRUE)
var_imp_res <- var_imp[var_imp$Overall > 7, ]
library(ggplot2)
p <- ggplot(var_imp_res, aes(x=Predictor, y=Overall))
p <- p + geom_bar(stat='identity') + coord_flip()
p <- p + ylab("Variable Importance") + xlab("")
pngname <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/N_application_subset_2000_2015/WB_station_orig_5min_pixel/no_climate_zones/var_imp_plot.png"
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
delta_fert_by_stn <- merge(DELTA_FERT_DF_GR, STN_OBJ_MATCH_DF_GR)
fert_app_subset_stn_list <- delta_fert_by_stn[
  delta_fert_by_stn$perc_change_N_application_1980_2013 <=
  GROUND_FERT_MEDIAN_VAL, 'GEMS.Stati']
subset_stn_list_time <- intersect(ground_df_stn_list, time_subset_stn_list)
subset_stn_list <- intersect(subset_stn_list_time, fert_app_subset_stn_list)

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
stn_covar_cols <- colnames(stn_covar_df)[c(1, 3:5, 8:13, 15)]  # pixel covar cols for groundwater (no basin area covar)
stn_covar_df <- stn_covar_df[, stn_covar_cols]

# merge groundwater covariates with noxn obs
combined_df <- merge(stn_covar_df, noxn_obs_restr, by.x='GEMS.Stati', by.y='GEMS.Station.Number')

# write out groundwater raw data for random forests analysis in Python
combined_df_restr <- combined_df[, 2:12]
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
out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR/ground/all_stations/WB_station_orig_5min_pixel"
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

##### ARCHIVE CODE ######
# station data from World Bank
# parameter metadata: observation descriptions
parameter_metadata_csv <- "F:/NCI_NDR/Data worldbank/station_data/parameter_metadata.csv"
parameter_metadata_df <- read.csv(parameter_metadata_csv)

# observations data related to Nitrogen
N_files_dir <- "F:/NCI_NDR/Data worldbank/station_data/Nitrogen_sheets_exported"
N_files_list <- list.files(N_files_dir)
N_df_list <- list()
for (file in N_files_list) {
  df <- read.csv(paste(N_files_dir, file, sep='/'))
  N_df_list[[file]] <- df
}
N_measurement_df <- do.call(rbind, N_df_list)
N_measurement_df$Sample.Date <- as.Date(N_measurement_df$Sample.Date, format="%Y-%m-%d")
Noxn_cols <- colnames(N_measurement_df)[c(1:2, 8:9)]
Noxn_obs_df <- N_measurement_df[N_measurement_df$Parameter.Code == "NOxN", Noxn_cols]
write.csv(Noxn_obs_df, NOXN_OBSERVATIONS_CSV, row.names=FALSE)

# stations with surface NOxN observations, all dates
obs_subset_stn_df <- stn_subset[stn_subset$GEMS.Station.Number %in% surface_df_stn_list, ]
write.csv(obs_subset_stn_df, SURFACE_NOXN_STATION_CSV,
          row.names=FALSE)

gr_subset_stn_df <- stn_subset[
  stn_subset$GEMS.Station.Number %in% ground_df_stn_list, ]
write.csv(gr_subset_stn_df, GR_NOXN_STATION_CSV)
