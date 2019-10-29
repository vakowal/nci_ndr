# random forests model: NOxN ~ NDR for World Bank

# match table: objectid+GEMS station number, all surface stations with noxn observations
OBJECTID_MATCH_CSV <- "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_objectid_stnid_match_table.csv"
STN_OBJ_MATCH_DF <- read.csv(OBJECTID_MATCH_CSV)

# stations with noxn observations
SURFACE_NOXN_STATION_CSV <- "F:/NCI_NDR/Data worldbank/station_data/WB_surface_stations_noxn_obs.csv"

# stations with change in N fertilizer application rates, 1980-2013
DELTA_FERTILIZER_1980_2013_CSV <- "F:/NCI_NDR/Data fertilizer Lu Tian/perc_change_mean_N_application_1980_2013.csv"
DELTA_FERT_DF <- read.csv(DELTA_FERTILIZER_1980_2013_CSV)

# NOxN observations, all stations
NOXN_OBSERVATIONS_CSV <- "F:/NCI_NDR/Data worldbank/station_data/noxn_obs_all.csv"
Noxn_obs_df <- read.csv(NOXN_OBSERVATIONS_CSV, stringsAsFactors=FALSE)
Noxn_obs_df$Sample.Date <- as.Date(Noxn_obs_df$Sample.Date, format="%Y-%m-%d")

# dates by which to restrict NOxN observations
MIN_DATE = "1990-01-01"
MAX_DATE = "1999-12-31"

# covariate table: basin extent
BASIN_COVAR_CSV <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/Rafa_watersheds_v3/combined_covariates.csv"

# covariate table: 5 min pixel extent
PIXEL_COVAR_CSV <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/WB_station_5min_pixel/combined_covariates.csv"

library(ggplot2)

# summarize NOXN training data from World bank
training_data_csv <- "F:/NCI_NDR/Data worldbank/training_data/NOXN/noxn_04_03_19.csv"
training_data <- read.csv(training_data_csv)

# collect min and max values from training data by column
min_by_column <- as.data.frame(apply(training_data, 2, min))
min_by_column$variable <- row.names(min_by_column)
colnames(min_by_column) <- c('min_value', 'variable')
max_by_column <- as.data.frame(apply(training_data, 2, max))
max_by_column$variable <- row.names(max_by_column)
colnames(max_by_column) <- c('max_value', 'variable')
col_summary <- merge(min_by_column, max_by_column)
write.csv(col_summary, "F:/NCI_NDR/Data worldbank/training_data/NOXN/min_max_values.csv")

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

# station metadata
station_metadata_csv <- "F:/NCI_NDR/Data worldbank/station_data/station_metadata.csv"
stn_df <- read.csv(station_metadata_csv)

# how many stations of each type (groundwater vs surface water) have noxn measurements?
stn_cols <- colnames(stn_df)[c(1, 4, 5, 8:10, 15:16)]
stn_subset <- stn_df[, stn_cols]
stn_subset[stn_subset$Water.Type == 'Groundwater station', 'ground_v_surface'] <- 'groundwater'
stn_subset[stn_subset$Water.Type == 'Lake station', 'ground_v_surface'] <- 'surface'
stn_subset[stn_subset$Water.Type == 'River station', 'ground_v_surface'] <- 'surface'
NOxN_by_stn <- merge(Noxn_obs_df, stn_subset, all.x=TRUE)
NOxN_by_stn$Sample.Date <- as.Date(NOxN_by_stn$Sample.Date, format="%Y-%m-%d")

# groundwater stations with NOxN measurements
NOxN_groundwater <- NOxN_by_stn[NOxN_by_stn$ground_v_surface == 'groundwater', ]
groundwater_freq_table <- as.data.frame(table(NOxN_groundwater$GEMS.Station.Number))
groundwater_freq_table <- groundwater_freq_table[groundwater_freq_table$Freq > 0, ]
colnames(groundwater_freq_table) <- c('GEMS.Station.Number', 'Num_obs_NOxN')
stn_coords <- stn_subset[stn_subset$GEMS.Station.Number %in% groundwater_freq_table$GEMS.Station.Number,
                         c("GEMS.Station.Number", "Country.Name", "Latitude", "Longitude")]
groundwater_freq_table <- merge(groundwater_freq_table, stn_coords)
write.csv(groundwater_freq_table, "F:/NCI_NDR/Data worldbank/station_data/groundwater_noxn_obs.csv",
          row.names=FALSE)

# restrict to surface stations
surface_df_stn_list <- unique(NOxN_by_stn[NOxN_by_stn$ground_v_surface == 'surface',
                                               'GEMS.Station.Number'])
# restrict surface stations by date
time_subset_stn_list <- unique(Noxn_obs_df[(Noxn_obs_df$Sample.Date >= MIN_DATE) &
                                         (Noxn_obs_df$Sample.Date <= MAX_DATE), 'GEMS.Station.Number'])
# restrict by stability of N fertilizer application rates
delta_fert_by_stn <- merge(DELTA_FERT_DF, STN_OBJ_MATCH_DF)
fert_app_subset_stn_list <- delta_fert_by_stn[delta_fert_by_stn$perc_change_mean_N_application_1980_2013 <= 155.4,
                                          'GEMS.Station.Number']
# intersection of subsets:
  # - surface
  # - time period
  # - stable trend in fertilizer application rates
subset1_stn_list <- intersect(surface_df_stn_list, time_subset_stn_list)
subset_stn_list <- intersect(subset1_stn_list, fert_app_subset_stn_list)

# stations with surface NOxN observations, all dates
obs_subset_stn_df <- stn_subset[stn_subset$GEMS.Station.Number %in% surface_df_stn_list, ]
write.csv(obs_subset_stn_df, SURFACE_NOXN_STATION_CSV,
          row.names=FALSE)

# merge covariate table with station characteristics
covar_df <- read.csv(BASIN_COVAR_CSV)
match_df <- read.csv(OBJECTID_MATCH_CSV)
covar_df <- merge(covar_df, match_df)
stn_df <- read.csv(SURFACE_NOXN_STATION_CSV)
stn_covar_df <- merge(covar_df, stn_df)

# add dummy variables for each water body type
stn_covar_df$lake <- 0
stn_covar_df[stn_covar_df$Water.Type == 'Lake station', 'lake'] <- 1
stn_covar_df$river <- 0
stn_covar_df[stn_covar_df$Water.Type == 'River station', 'river'] <- 1
stn_covar_cols <- colnames(stn_covar_df)[c(1, 3:12, 21)]
stn_covar_df <- stn_covar_df[, stn_covar_cols]
# get one observation per station, according to latest date of measurement inside selected time period
noxn_obs_subset <- Noxn_obs_df[(Noxn_obs_df$Sample.Date >= MIN_DATE) &
                                 (Noxn_obs_df$Sample.Date <= MAX_DATE), ]
max_date_df <- aggregate(Sample.Date~GEMS.Station.Number, data=noxn_obs_subset, FUN=max)
max_date_obs <- merge(noxn_obs_subset, max_date_df)
max_date_avg <- aggregate(Value~GEMS.Station.Number+Sample.Date+Unit,
                          data=max_date_obs, FUN=mean)
noxn_obs_restr <- max_date_avg[, c(1, 4)]
colnames(noxn_obs_restr)[2] <- 'noxn'
# merge covariates with NOxN observations
combined_df <- merge(stn_covar_df, noxn_obs_restr, by="GEMS.Station.Number")  # , all=TRUE)

# restrict by stability of N fert application
combined_df_restr <- combined_df[combined_df$GEMS.Station.Number %in% subset_stn_list, ]
combined_df_restr$climate_zone <- as.factor(combined_df_restr$climate_zone)
combined_df_restr$lake <- as.factor(combined_df_restr$lake)

############ random forests model ############
library(lattice)
library(mlbench)
library(caret)

covariate_vals <- combined_df_restr[, c(2:12)]
noxn <- combined_df_restr$noxn

# visualize bivariate relationships of noxn to all covariates
theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
featurePlot(x = combined_df_restr[, c(2:12)], 
            y = combined_df_restr$noxn, 
            plot = "scatter",
            type=c('p', 'smooth'),
            span=0.5,
            layout = c(5, 3))

# check for near zero variance covariates
nzv <- nearZeroVar(covariate_vals, saveMetrics=TRUE)
nzv[nzv$nzv, ][1:10, ]
# check for high correlations between covariates
covar_cor <- cor(covariate_vals)
summary(covar_cor[upper.tri(covar_cor)])
# center and scale continuous covariates
preProcValues <- preProcess(covariate_vals, method=c("center", "scale"))
covarTransformed <- predict(preProcValues, newdata=covariate_vals)
# K-fold cross-validation
set.seed(491)
fitControl <- trainControl(method='repeatedcv', number=10, repeats=10)
mtry <- floor(sqrt(dim(covariate_vals)[2]))  # number randomly selected predictors
# fit the random forests model
processed_data <- cbind(covariate_vals, noxn)
ranger_rf <- train(noxn ~ ., data=processed_data,
                   method='ranger', trControl=fitControl,
                   na.action=na.omit, importance='impurity')
ranger_rf
rf_var_imp <- varImp(ranger_rf)
