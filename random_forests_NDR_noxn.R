# random forests model: NOxN ~ NDR for World Bank

# prepare nitrate and covariate data for modeling with random forest
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
PIXEL_ORIG_COVAR_CSV <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/WB_station_orig_5min_pixel/combined_covariates.csv"
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
colnames(noxn_obs_restr) <- c('OBJECTID', 'noxn')
noxn_obs_path <- "C:/Users/ginge/Documents/Python/nci_ndr/noxn_obs_surface.csv"
write.csv(noxn_obs_restr, noxn_obs_path)

###### process covariate data: 5min pixel level (not snapped) ######
OBJECTID_MATCH_CSV_SURF <- "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_objectid_stnid_match_table.csv"
SURFACE_NOXN_STATION_CSV <- "F:/NCI_NDR/Data worldbank/station_data/WB_surface_stations_noxn_obs.csv"
# 4 options for NDR outputs, 3.19.21
covar_dir <- 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates'
PIXEL_ORIG_COVAR_CSV <- paste(covar_dir, 'surface_WB_station_orig_5min_pixel_N_export_via_R', 'combined_covariates.csv', sep='/')
PIXEL_ORIG_COVAR_CSV <- paste(covar_dir, 'surface_WB_station_orig_5min_pixel_N_export_via_GDAL', 'combined_covariates.csv', sep='/')
PIXEL_ORIG_COVAR_CSV <- paste(covar_dir, 'surface_WB_station_orig_5min_pixel_N_load_via_R', 'combined_covariates.csv', sep='/')
PIXEL_ORIG_COVAR_CSV <- paste(covar_dir, 'surface_WB_station_orig_5min_pixel_N_load_via_GDAL', 'combined_covariates.csv', sep='/')
# ^^ 4 options for NDR outputs

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
colnames(stn_covar_df)[1] <- 'OBJECTID'
# 5 min resolution
# covar_df_path <- "C:/Users/ginge/Documents/Python/nci_ndr/covar_df_5min_surface.csv"
# 4 options for NDR outputs, 3.19.21
outer_dir <- 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/covar_df_3.19.21'
covar_df_path <- paste(outer_dir, 'surface_WB_station_orig_5min_pixel_N_export_via_R', 'covar_df_surface.csv', sep='/')
covar_df_path <- paste(outer_dir, 'surface_WB_station_orig_5min_pixel_N_export_via_GDAL', 'covar_df_surface.csv', sep='/')
covar_df_path <- paste(outer_dir, 'surface_WB_station_orig_5min_pixel_N_load_via_R', 'covar_df_surface.csv', sep='/')
covar_df_path <- paste(outer_dir, 'surface_WB_station_orig_5min_pixel_N_load_via_GDAL', 'covar_df_surface.csv', sep='/')

write.csv(stn_covar_df, covar_df_path)

# 30 s resolution
# covar_df_path <- "C:/Users/ginge/Documents/Python/nci_ndr/covar_df_30s_surface.csv"
# write.csv(stn_covar_df, covar_df_path)

## add accumulated and accumulated-normalized predictors
accum_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/accumulated/WB_station_orig_5min_pixel/combined_covariates.csv")
STN_DF <- read.csv(SURFACE_NOXN_STATION_CSV)
stn_accum_df <- merge(accum_df, STN_DF)
accum_df_objectid <- stn_accum_df[, c('GEMS.Stati', colnames(accum_df))]
accum_df_objectid <- subset(accum_df_objectid, select=-c(OBJECTID))
colnames(accum_df_objectid)[1] <- 'OBJECTID'
local_covar_df_path <- "C:/Users/ginge/Documents/Python/nci_ndr/covar_df_5min_surface.csv"
local_covar_df <- read.csv(local_covar_df_path)
nload_covar_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/N_application_predictor/WB_station_orig_5min_pixel/combined_covariates.csv")
nload_stn_covar <- merge(nload_covar_df, STN_DF)
nload_stn_covar <- nload_stn_covar[, c('GEMS.Stati', 'n_export')]
colnames(nload_stn_covar) <- c('OBJECTID', 'n_load')
combined_df <- merge(accum_df_objectid, local_covar_df)
combined_df <- merge(combined_df, nload_stn_covar)
combined_df <- subset(combined_df, select=-c(X))
write.csv(combined_df, "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/accumulated/WB_station_orig_5min_pixel/covar_df_combined_with_local.csv")

# examine correlations among predictor variables, including accumulated
library(corrplot)
library(caret)
combined_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/accumulated/WB_station_orig_5min_pixel/covar_df_combined_with_local.csv")
surface_predictors <- subset(combined_df, select=-c(X, OBJECTID))
x1_cor <- cor(surface_predictors, use="complete.obs")
corrplot(x1_cor, type='upper', diag=FALSE)
corr_vars <- findCorrelation(x1_cor, 0.6)
variable.names(surface_predictors[,corr_vars])
pngname <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/accumulated/WB_station_orig_5min_pixel/covariate_correlation_plot_surface.png"
png(file=pngname, units="in", res=300, width=6, height=6)
corrplot(x1_cor, type='upper', diag=FALSE)
dev.off()

# random forests model, surface stations
# do this in train_evaluate_random_forest.R
# out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/surface/subset_2000_2015/WB_station_orig_5min_pixel"
# N application predictor
# out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/N_app_predictor/surface/subset_2000_2015/WB_station_orig_5min_pixel"
# 30 sec scale
# out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/surface/subset_2000_2015/WB_station_orig_30s_pixel"
# train_evaluate_rf(noxn_obs_path, covar_df_path, out_dir)

# Bar chart: variable importance (for an example model)
var_imp_csv <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/ground/GEMStat_China_USGS_Ouedraogo/var_imp.csv"
var_imp <- read.csv(var_imp_csv)
colnames(var_imp) <- c('variable', 'importance')
var_imp$Predictor <- reorder(var_imp$variable, var_imp$importance, descending=TRUE)
# var_imp_res <- var_imp[var_imp$Overall > 7, ]
library(ggplot2)
p <- ggplot(var_imp, aes(x=Predictor, y=importance))
p <- p + geom_bar(stat='identity') + coord_flip()
p <- p + ylab("Variable Importance") + xlab("")
pngname <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/ground/GEMStat_China_USGS_Ouedraogo/var_imp_plot.png"
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
noxn_obs <- read.csv(noxn_obs_path)
write.csv(noxn_obs, "C:/Users/ginge/Documents/Python/nci_ndr/noxn_obs_GEMStat_Gu_USGS_Ouedraogo_ground.csv")
# N application as predictor: covar_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/N_application_predictor/groundwater_GEMStat_Gu_USGS_Ouedraogo/combined_covariates.csv")
# 5 min resolution
covar_df_path <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/groundwater_GEMStat_Gu_USGS_Ouedraogo/combined_covariates.csv"
covar_df <- read.csv(covar_df_path)
covar_df <- covar_df[, GROUNDWATER_PREDICTORS]
covar_df_path_5min <- "C:/Users/ginge/Documents/Python/nci_ndr/covar_df_5min_ground.csv"
write.csv(covar_df, covar_df_path_5min)

# 4 options for NDR outputs, 3.19.12
covar_dir <- 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates'
covar_df_path <- paste(covar_dir, 'groundwater_GEMStat_Gu_USGS_Ouedraogo_5min_pixel_N_export_via_R', 'combined_covariates.csv', sep='/')
covar_df_path <- paste(covar_dir, 'groundwater_GEMStat_Gu_USGS_Ouedraogo_5min_pixel_N_export_via_GDAL', 'combined_covariates.csv', sep='/')
covar_df_path <- paste(covar_dir, 'groundwater_GEMStat_Gu_USGS_Ouedraogo_5min_pixel_N_load_via_R', 'combined_covariates.csv', sep='/')
covar_df_path <- paste(covar_dir, 'groundwater_GEMStat_Gu_USGS_Ouedraogo_5min_pixel_N_load_via_GDAL', 'combined_covariates.csv', sep='/')

covar_df <- read.csv(covar_df_path)
covar_df <- covar_df[, GROUNDWATER_PREDICTORS]
outer_dir <- 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/covar_df_3.19.21'

out_path <- paste(outer_dir, 'groundwater_GEMStat_Gu_USGS_Ouedraogo_5min_pixel_N_export_via_R', 'covar_df_ground.csv', sep='/')
out_path <- paste(outer_dir, 'groundwater_GEMStat_Gu_USGS_Ouedraogo_5min_pixel_N_export_via_GDAL', 'covar_df_ground.csv', sep='/')
out_path <- paste(outer_dir, 'groundwater_GEMStat_Gu_USGS_Ouedraogo_5min_pixel_N_load_via_R', 'covar_df_ground.csv', sep='/')
out_path <- paste(outer_dir, 'groundwater_GEMStat_Gu_USGS_Ouedraogo_5min_pixel_N_load_via_GDAL', 'covar_df_ground.csv', sep='/')

write.csv(covar_df, out_path)
# ^^ 4 options

## add accumulated and accumulated-normalized predictors
accum_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/accumulated/groundwater_GEMStat_Gu_USGS_Ouedraogo/combined_covariates.csv")
nload_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/N_application_predictor/groundwater_GEMStat_Gu_USGS_Ouedraogo/combined_covariates.csv")
nload_df <- nload_df[, c('OBJECTID', 'n_export')]
colnames(nload_df) <- c('OBJECTID', 'n_load')
combined_df <- merge(accum_df, covar_df)
combined_df <- merge(combined_df, nload_df)
write.csv(combined_df, "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/accumulated/groundwater_GEMStat_Gu_USGS_Ouedraogo/covar_df_combined_with_local.csv")

combined_df <- read.csv("C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/accumulated/groundwater_GEMStat_Gu_USGS_Ouedraogo/covar_df_combined_with_local.csv")
ground_predictors <- subset(combined_df, select=-c(X, OBJECTID))
x1_cor <- cor(ground_predictors, use="complete.obs")
corrplot(x1_cor, type='upper', diag=FALSE)
corr_vars <- findCorrelation(x1_cor, 0.6)
variable.names(ground_predictors[,corr_vars])
pngname <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/accumulated/groundwater_GEMStat_Gu_USGS_Ouedraogocovariate_correlation_plot_ground.png"
png(file=pngname, units="in", res=300, width=6, height=6)
corrplot(x1_cor, type='upper', diag=FALSE)
dev.off()

# out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/ground/GEMStat_China_USGS_Ouedraogo"
# out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/N_app_predictor/ground/GEMStat_China_USGS_Ouedraogo"
# 30 s scale
covar_df_30s_full <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/groundwater_GEMStat_Gu_USGS_Ouedraogo_30s_pixel/combined_covariates.csv"
covar_df <- read.csv(covar_df_30s_full)
covar_df <- covar_df[, GROUNDWATER_PREDICTORS]
covar_df_path_30s <- "C:/Users/ginge/Documents/Python/nci_ndr/covar_df_30s_ground.csv"
write.csv(covar_df, covar_df_path_30s)
out_dir <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Analysis_results/Updated_NDR_5.18.20/ground/GEMStat_China_USGS_Ouedraogo_30s"
train_evaluate_rf(noxn_obs_path, covar_df_path, out_dir)
