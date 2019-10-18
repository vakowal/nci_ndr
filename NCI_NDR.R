# data prep for NCI Nitrate ~ NDR for World Bank

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
  
library(ggplot2)

# NOXN training data from World bank
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

Noxn_obs_df <- read.csv(NOXN_OBSERVATIONS_CSV)
Noxn_obs_df$Sample.Date <- as.Date(Noxn_obs_df$Sample.Date, format="%Y-%m-%d")

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
MIN_DATE = "1990-01-01"
MAX_DATE = "1999-12-31"
time_subset_stn_list <- unique(Noxn_obs_df[(Noxn_obs_df$Sample.Date >= MIN_DATE) &
                                         (Noxn_obs_df$Sample.Date >= MAX_DATE), 'GEMS.Station.Number'])
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
obs_subset_stn_df <- stn_subset[stn_subset$GEMS.Station.Number %in% noxn_surface_obs_station_list, ]
write.csv(obs_subset_stn_df, SURFACE_NOXN_STATION_CSV,
          row.names=FALSE)

## raster processing
# memory-intensive raster processing for NCI_NDR
library(raster)

# generate aligned raster of pixel area values, for calculating area of watersheds
# aligned with DEM that was used to delineate the watersheds
template_raster_path <- "F:/NCI_NDR/Data flow direction DRT/globe_dem_shedsandh1k.asc"
area_km2_raster <- area(raster(template_raster_path))
writeRaster(area_km2_raster, "F:/NCI_NDR/Data flow direction DRT/pixel_area_km2.tif")

# compare basin size calculated via zonal statistics from delineated basins to reported basin size
basin_area_square_km_csv <- "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_for_snapping_ContribArea.csv"
basin_area_square_km_df <- read.csv(basin_area_square_km_csv)
basin_area_square_km_df <- basin_area_square_km_df[, colnames(basin_area_square_km_df)[
  c(1, 5, 10, 34)]]
colnames(basin_area_square_km_df) <- c("OBJECTID", "GEMS.Station.Number", "Upstream.Basin.Area", "Basin_area_square_km")

basin_area_square_km_df$reported_div_calc <- basin_area_square_km_df$Upstream.Basin.Area / basin_area_square_km_df$Basin_area_square_km
basin_area_square_km_df$perc_diff <- (
  (basin_area_square_km_df$Basin_area_square_km - basin_area_square_km_df$Upstream.Basin.Area) / basin_area_square_km_df$Upstream.Basin.Area)
basin_area_square_km_df$abs_perc_diff <- (
  abs(basin_area_square_km_df$Basin_area_square_km - basin_area_square_km_df$Upstream.Basin.Area) / basin_area_square_km_df$Upstream.Basin.Area)
dim(basin_area_square_km_df[!is.na(basin_area_square_km_df$Upstream.Basin.Area), ]) # 325 stations have reported basin size
success_subs <- basin_area_square_km_df[
  !is.na(basin_area_square_km_df$Upstream.Basin.Area) & (
    basin_area_square_km_df$abs_perc_diff < 0.2), ]
dim(success_subs)  # 305 stations have area within 50% of reported size
write.csv(basin_area_square_km_df, "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/compare_basin_area_calc_v_reported.csv")

# raster processing
library(raster)
# resample NDR outputs with block statistics
target_pixel_size <- 0.08333333333333332871  # 5 min
n_export_10s_tif <- "F:/NCI_NDR/Data NDR/nutrient_deficit_10s_cur_compressed_md5_031d4bb444325835315a2cc825be3fd4.tif"
n_export_5min_tif <- "F:/NCI_NDR/Data NDR/nutrient_deficit_5min_cur_compressed_md5_031d4bb444325835315a2cc825be3fd4.tif"
in_ras <- raster(n_export_10s_tif)
in_pixel_size = xres(in_ras)
aggregate_factor = round(target_pixel_size / in_pixel_size)
agg_ras <- aggregate(
  in_ras, fact=aggregate_factor, fun=sum, expand=TRUE, na.rm=TRUE)
out_ras <- reclassify(agg_ras, cbind(NA, -9999))
writeRaster(out_ras, filename=n_export_5min_tif, NAflag=-9999)

# aggregate across bands inside FLO1K NetCDF
min_year <- 1990
min_band <- min_year - 1959
temporal_subset = c(min_band:56)

# average of average flow, 1990-2015
flo1k_av_ncdf <- "F:/NCI_NDR/Data streamflow FLO1K/FLO1K.5min.ts.1960.2015.qav.nc"
flo1k_stack <- stack(flo1k_av_ncdf)
flo1k_subs <- flo1k_stack[[temporal_subset]]
subs_mean <- calc(flo1k_subs, fun=mean)
writeRaster(subs_mean, "F:/NCI_NDR/Data streamflow FLO1K/average_flow_1990_2015.tif")

# average minimum flow across years, 1990-2015
flo1k_min_ncdf <- "F:/NCI_NDR/Data streamflow FLO1K/FLO1K.5min.ts.1960.2015.qmi.nc"
flo1k_min_stack <- stack(flo1k_min_ncdf)
flo1k_min_subs <- flo1k_min_stack[[temporal_subset]]
min_subs_mean <- calc(flo1k_min_subs, fun=mean)
writeRaster(min_subs_mean, "F:/NCI_NDR/Data streamflow FLO1K/average_min_flow_1990_2015.tif")

# average maximum flow across years, 1990-2015
flo1k_max_ncdf <- "F:/NCI_NDR/Data streamflow FLO1K/FLO1K.5min.ts.1960.2015.qma.nc"
flo1k_max_stack <- stack(flo1k_max_ncdf)
flo1k_max_subs <- flo1k_max_stack[[temporal_subset]]
max_subs_mean <- calc(flo1k_max_subs, fun=mean)
writeRaster(max_subs_mean, "F:/NCI_NDR/Data streamflow FLO1K/average_max_flow_1990_2015.tif")

# ratio of average min to average, 1990-2015
div_fun <- function(x, y) { return (x/y) }
min_avg_ratio <- overlay(min_subs_mean, subs_mean, fun=div_fun)
writeRaster(min_avg_ratio, "F:/NCI_NDR/Data streamflow FLO1K/min_div_average_flow_1990_2015.tif")

# calculate mean irrigation area across years from HYDE3.2 dataset
irri_raster_list = list()
irri_path_pattern <- "F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri<year>AD.asc"
irri_year_sequence <- c(1990, 2000:2015)
i <- 1
for (year in irri_year_sequence) {
  filename <- gsub('<year>', year, irri_path_pattern)
  irri_raster_list[[i]] <- raster(filename)
  i <- i + 1
}
irri_stack <- stack(irri_raster_list)
irri_mean <- calc(irri_stack, fun=mean)
writeRaster(irri_mean, "F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri_avg_1990_2015.tif")

# pixel area of irrigation rasters
irri_area_km2_raster <- area(irri_mean)
writeRaster(irri_area_km2_raster, "F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri_pixel_area_km2.tif")

# mean population across years, from HYDE3.2 dataset
pop_raster_list = list()
pop_path_pattern <- "F:/NCI_NDR/Data population HYDE3.2/popd_<year>AD.asc"
pop_year_sequence <- c(1990, 2000:2015)
i <- 1
for (year in pop_year_sequence) {
  filename <- gsub('<year>', year, pop_path_pattern)
  pop_raster_list[[i]] <- raster(filename)
  i <- i + 1
}
pop_stack <- stack(pop_raster_list)
pop_mean <- calc(pop_stack, fun=mean)
# multiply by pixel area to get total inhabitants per pixel
area_raster <- raster("F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri_pixel_area_km2.tif")
inhabitants <- pop_mean * area_raster
writeRaster(inhabitants, "F:/NCI_NDR/Data population HYDE3.2/inhabitants_avg_1990_2015.tif")

# change in fertilizer application across years within watersheds
# from Lu and Tian 2017
fert_by_basin_year_csv <- "F:/NCI_NDR/Data fertilizer Lu Tian/N_by_OBJECTID_WB_surface_stations_noxn_for_snapping_ContribArea.csv"
fert_by_basin_df <- read.csv(fert_by_basin_year_csv)
fert_by_basin_df$mean <- fert_by_basin_df$sum / fert_by_basin_df$count

library(ggplot2)
since_1980_subset <- fert_by_basin_df[fert_by_basin_df$year >= 1980, ]
p <- ggplot(since_1980_subset, aes(x=year, y=mean, group=OBJECTID))
p <- p + geom_line(aes(color=OBJECTID))
print(p)

# change in mean application rate, 2013-1980
start_end_subs <- fert_by_basin_df[(fert_by_basin_df$year==1980) |
                                     (fert_by_basin_df$year==2013), colnames(fert_by_basin_df)[c(1, 7:8)]]
delta_1980_2013_df <- reshape(start_end_subs, idvar='OBJECTID',
                              timevar='year', v.names='mean', direction='wide')
delta_1980_2013_df$percent_change <- ((delta_1980_2013_df$mean.2013 - delta_1980_2013_df$mean.1980)/
                                        delta_1980_2013_df$mean.1980) * 100
delta_1980_2013_df[!is.na(delta_1980_2013_df$mean.1980) & delta_1980_2013_df$mean.1980==0, ]$percent_change <- 100
delta_1980_2013_df[!is.na(delta_1980_2013_df$delta) & delta_1980_2013_df$delta==0, ]$percent_change <- 100
delta_1980_2013_df[is.na(delta_1980_2013_df$percent_change), ]$percent_change <- 0
delta_df_to_save <- delta_1980_2013_df[, c("OBJECTID", "percent_change")]
colnames(delta_df_to_save)[2] <- 'perc_change_mean_N_application_1980_2013'
write.csv(delta_df_to_save, DELTA_FERTILIZER_1980_2013_CSV,
          row.names=FALSE)

# generate a table of stations/basins that can be used to subset by delta N application and date of observation
n_app_stn_df <- merge(delta_df_to_save, STN_OBJ_MATCH_DF)

# summarize sanitation data
sanitation_csv <- "F:/NCI_NDR/Data sanitation/sanitation_data_raw.csv"
sanitation_df <- data.frame(read.csv(sanitation_csv))
sanitation_cols_to_keep <- c("year", "iso3", "san_lim_n", "san_unimp_n", "san_od_n")
san_subs_df <- sanitation_df[, sanitation_cols_to_keep]
colnames(san_subs_df) <- c('year', 'ISO3', 'limited', 'unimproved', 'open_defecation')
san_subs_df$no_san_provision <- rowSums(san_subs_df[, c('limited', 'unimproved', 'open_defecation')],
                                        na.rm=TRUE)
san_subs_df[(is.na(san_subs_df$unimproved) & is.na(san_subs_df$limited) & is.na(san_subs_df$open_defecation)),
            'no_san_provision'] <- NA
san_subs_df <- san_subs_df[san_subs_df$year < 2016, ]

san_avg_2000_2015 <- aggregate(no_san_provision~ISO3, data=san_subs_df, FUN=mean)
# match these values with countryid raster via ISO3 country code field
countryid_iso3_match_csv <- "F:/NCI_NDR/Data world borders/countryid_iso3_match_table_tm_world_borders_0.3.csv"
countryid_match_df <- read.csv(countryid_iso3_match_csv)
san_avg_2000_2015_countryid <- merge(san_avg_2000_2015, countryid_match_df, all=TRUE)
# remove countries from JMP database that don't appear in countries raster
san_avg_2000_2015_countryid <- san_avg_2000_2015_countryid[!(is.na(san_avg_2000_2015_countryid$countryid)), ]

write.csv(san_avg_2000_2015_countryid, "F:/NCI_NDR/Data sanitation/no_sanitation_provision_avg_2000-2015.csv",
          row.names=FALSE)

# merge covariate table with station characteristics
covariate_csv <- "C:/Users/ginge/Desktop/combined_covariates.csv"
covar_df <- read.csv(covariate_csv)
match_df <- read.csv(OBJECTID_MATCH_CSV)
covar_df <- merge(covar_df, match_df)
stn_df <- read.csv(SURFACE_NOXN_STATION_CSV)
stn_covar_df <- merge(covar_df, stn_df)

