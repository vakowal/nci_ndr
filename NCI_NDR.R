# Spatial data processing for NCI-NDR project
library(raster)

# table giving percent change in mean fertilizer application rate between modeled period and prior period
# surface
FERT_PERC_CHANGE_CSV_SURF <- "F:/NCI_NDR/Data fertilizer Lu Tian/perc_change_fert_surf.csv"
# ground
FERT_PERC_CHANGE_CSV_GR <- "F:/NCI_NDR/Data fertilizer Lu Tian/perc_change_fert_gr.csv"
  
# match table: objectid+GEMS station number, all surface stations with noxn observations
OBJECTID_MATCH_CSV <- "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_objectid_stnid_match_table.csv"
STN_OBJ_MATCH_DF <- read.csv(OBJECTID_MATCH_CSV)

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

# resample NDR outputs to 5min resolution with block statistics
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

# resample countries raster up to 5min resolution
covar_key_csv <- "C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/covariate_path_key.csv"
covar_key_df <- read.csv(covar_key_csv, stringsAsFactors=FALSE)
covar_key_df$native_resolution <- NA
for (r in c(1:NROW(covar_key_df))) {
  tryCatch({
    temp_ras <- raster(covar_key_df[r, 'path_native_resolution'])
    covar_key_df[r, 'native_resolution'] <- xres(temp_ras)
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
target_pixel_size <- 0.08333333333333332871  # 5 min
in_ras <- raster(covar_key_df[11, 'path_native_resolution'])  # countries raster
in_pixel_size = xres(in_ras)
aggregate_factor = round(target_pixel_size / in_pixel_size)
agg_ras <- aggregate(
  in_ras, fact=aggregate_factor, fun=modal, expand=TRUE, na.rm=TRUE)
out_ras <- reclassify(agg_ras, cbind(NA, -9999))
filename <- covar_key_df[r, 'path_5min_resolution']
writeRaster(out_ras, filename=filename, NAflag=-9999)

# reclassify urban extent raster, calculate % urban at 5min resolution
urban_ras <- raster(covar_key_df[9, 'path_native_resolution'])  # urban extent raster
# Reclassify from {1=rural and 2=urban} to {0 = rural and 1 = urban}.
rcl <- c(0.5, 1.5, 0, 1.6, 2.5, 1)
rcl_mat <- matrix(rcl, ncol=3, byrow=TRUE)
urban_rc <- reclassify(urban_ras, rcl_mat)
# calculate % urban at 5 min resolution
in_pixel_size = xres(urban_rc)
aggregate_factor = round(target_pixel_size / in_pixel_size)
agg_ras <- aggregate(
  urban_rc, fact=aggregate_factor, fun=mean, expand=TRUE, na.rm=TRUE)
filename <- "F:/NCI_NDR/Data urban extent GRUMP/perc_urban_5min.tif"
writeRaster(agg_ras, filename=filename, NAflag=-9999, overwrite=TRUE)

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
subs_mean <- raster("F:/NCI_NDR/Data streamflow FLO1K/average_flow_1990_2015.tif")

# average minimum flow across years, 1990-2015
flo1k_min_ncdf <- "F:/NCI_NDR/Data streamflow FLO1K/FLO1K.5min.ts.1960.2015.qmi.nc"
flo1k_min_stack <- stack(flo1k_min_ncdf)
flo1k_min_subs <- flo1k_min_stack[[temporal_subset]]
min_subs_mean <- calc(flo1k_min_subs, fun=mean)
writeRaster(min_subs_mean, "F:/NCI_NDR/Data streamflow FLO1K/average_min_flow_1990_2015.tif")
min_subs_mean <- raster("F:/NCI_NDR/Data streamflow FLO1K/average_min_flow_1990_2015.tif")

# average maximum flow across years, 1990-2015
flo1k_max_ncdf <- "F:/NCI_NDR/Data streamflow FLO1K/FLO1K.5min.ts.1960.2015.qma.nc"
flo1k_max_stack <- stack(flo1k_max_ncdf)
flo1k_max_subs <- flo1k_max_stack[[temporal_subset]]
max_subs_mean <- calc(flo1k_max_subs, fun=mean)
writeRaster(max_subs_mean, "F:/NCI_NDR/Data streamflow FLO1K/average_max_flow_1990_2015.tif")
max_subs_mean <- raster("F:/NCI_NDR/Data streamflow FLO1K/average_max_flow_1990_2015.tif")

# ratio of average min to average, 1990-2015
div_fun <- function(x, y) { return (x/y) }
min_avg_ratio <- overlay(min_subs_mean, subs_mean, fun=div_fun)
writeRaster(min_avg_ratio, "F:/NCI_NDR/Data streamflow FLO1K/min_div_average_flow_1990_2015.tif")

# ratio of mean flow to range (max - min), 'flashiness'
flash_fun <- function(mean, min, max) { return (mean / (max - min)) }
flash_ratio <- overlay(subs_mean, min_subs_mean, max_subs_mean, fun=flash_fun)
writeRaster(flash_ratio, "F:/NCI_NDR/Data streamflow FLO1K/mean_div_range_1990_2015.tif")

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

# change in fertilizer application across years at pixels containing stations
# from Lu and Tian 2017
MIN_YEAR = 2000
MAX_YEAR = 2015
library(ggplot2)
n_fert_surface_stn_csv = "F:/NCI_NDR/Data fertilizer Lu Tian/N_by_OBJECTID_WB_surface_stations_noxn_obs_objectid_adj.csv"
fert_surf_df <- read.csv(n_fert_surface_stn_csv)
since_1965 <- fert_surf_df[fert_surf_df$year >= 1965, ]
p <- ggplot(since_1965, aes(x=year, y=fertilizer, group=OBJECTID))
p <- p + geom_line(aes(colour=OBJECTID))
print(p)
# change, 1980-2013
start_end_subs <- fert_surf_df[(fert_surf_df$year==1980) |
                                 (fert_surf_df$year==2013), ]
delta_df <- reshape(start_end_subs, idvar='OBJECTID', timevar='year', v.names='fertilizer',
                    direction='wide')
delta_df$percent_change <- ((delta_df$fertilizer.2013 - delta_df$fertilizer.1980) /
                              delta_df$fertilizer.1980) * 100
delta_df[(delta_df$fertilizer.1980==0) & (delta_df$fertilizer.2013==0), 'percent_change'] <- 0
delta_df <- delta_df[, c('OBJECTID', 'percent_change')]
write.csv(delta_df, "F:/NCI_NDR/Data fertilizer Lu Tian/percent_change_1980_2013_surf.csv",
          row.names=FALSE)
  
model_period_surf_df <- fert_surf_df[
  (fert_surf_df$year >= MIN_YEAR) & (fert_surf_df$year <= MAX_YEAR), ]
mean_model_period <- aggregate(fertilizer~OBJECTID, data=model_period_surf_df, FUN=mean)
colnames(mean_model_period)[2] <- 'mean_fertilizer_model_period'
prior_period <- fert_surf_df[
  (fert_surf_df$year >= (MIN_YEAR - 5)) & (fert_surf_df$year <= MIN_YEAR), ]
mean_prior_period <- aggregate(fertilizer~OBJECTID, data=prior_period, FUN=mean)
colnames(mean_prior_period)[2] <- 'mean_fertilizer_prior_period'
change_df <- merge(mean_model_period, mean_prior_period)
change_df$perc_change_fert <- abs(
  (change_df$mean_fertilizer_model_period - change_df$mean_fertilizer_prior_period) / 
  change_df$mean_fertilizer_model_period) * 100
change_df[(change_df$mean_fertilizer_model_period == 0) &
            (change_df$mean_fertilizer_prior_period == 0), 'perc_change_fert'] <- 0
hist(change_df$perc_change_fert, breaks=100)
change_df <- change_df[, c('OBJECTID', 'perc_change_fert')]
write.csv(change_df, FERT_PERC_CHANGE_CSV_SURF, row.names=FALSE)

# fertilizer over time, groundwater stations
n_fert_groundwater_stn_csv = "F:/NCI_NDR/Data fertilizer Lu Tian/N_by_OBJECTID_WB_groundwater_stations_noxn_obs.csv"
fert_gr_df <- read.csv(n_fert_groundwater_stn_csv)
model_period <- fert_gr_df[
  (fert_gr_df$year >= 1995) & (fert_gr_df$year <= 2010), ]
mean_model_period <- aggregate(fertilizer~OBJECTID, data=model_period, FUN=mean)
colnames(mean_model_period)[2] <- 'mean_fertilizer_model_period'
prior_period <- fert_gr_df[
  (fert_gr_df$year >= 1965) & (fert_gr_df$year <= 1995), ]
mean_prior_period <- aggregate(fertilizer~OBJECTID, data=prior_period, FUN=mean)
colnames(mean_prior_period)[2] <- 'mean_fertilizer_prior_period'
change_df <- merge(mean_model_period, mean_prior_period)
change_df$perc_change_fert <- abs(
  (change_df$mean_fertilizer_model_period - change_df$mean_fertilizer_prior_period) / 
  change_df$mean_fertilizer_model_period) * 100
change_df[(change_df$mean_fertilizer_model_period == 0) &
            (change_df$mean_fertilizer_prior_period == 0), 'perc_change_fert'] <- 0
hist(change_df$perc_change_fert, breaks=100)
change_df <- change_df[, c('OBJECTID', 'perc_change_fert')]
write.csv(change_df, FERT_PERC_CHANGE_CSV_GR, row.names=FALSE)

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

# match country-level endpoints data with ISO3 country codes
country_codes <- read.csv("F:/NCI_NDR/Data national_boundaries/countries_iso3_match_table.csv")
cancer_rates <- read.csv("F:/NCI_NDR/Data endpoints/colorectal_cancer_risk_2018_per_person_by_country.csv",
                         stringsAsFactors=FALSE)
water_sources <- read.csv("F:/NCI_NDR/Data endpoints/Freshwater_groundwater_World_Bank_countries.csv")

cancer_rates <- cancer_rates[, c("Country.Name", "World.Bank.Country.Code", "Rate.per.person")]
water_sources <- water_sources[, c("Country.Name", "Country.Code.WB", "Country.Code.shapefile", "X..Surface.water", "X..Ground.water")]
colnames(water_sources) <- c("Country.Name", "World.Bank.Country.Code", "iso3", "perc_surface", "perc_ground")

water_join <- merge(water_sources, country_codes, by='iso3', all=TRUE)  # water source table contains ISO3 codes
water_join <- water_join[!(is.na(water_join$id)), ]  # drop 2 countries with values in water source table missing from shapefile
water_join[is.na(water_join$perc_surface), 'perc_surface'] <- 0.5  # give missing values the background value
water_source_by_country_id <- water_join[, c('iso3', 'id', 'perc_surface')]
write.csv(water_source_by_country_id, "F:/NCI_NDR/Data endpoints/water_source_by_country_id.csv",
          row.names=FALSE)

cancer_good <- cancer_rates[cancer_rates$World.Bank.Country.Code != 'na', ]
cancer_good[cancer_good$Country.Name == 'Namibia', 'World.Bank.Country.Code'] <- 'NAM'
cancer_good <- cancer_good[cancer_good$World.Bank.Country.Code != 'PSE', ]
cancer_good <- cancer_good[cancer_good$World.Bank.Country.Code != 'XKX', ]
cancer_join1 <- merge(cancer_good, water_join, by='World.Bank.Country.Code', all.x=TRUE)
cancer_join1 <- cancer_join1[, c('Rate.per.person', 'iso3', 'id')]

# figuring out codes this way
need_code <- cancer_rates[cancer_rates$World.Bank.Country.Code == 'na', ]
cancer_fill <- merge(need_code, country_codes, by.x='Country.Name', by.y='nev_name', all.x=TRUE)
# add a few country codes manually
cancer_fill[cancer_fill$Country.Name == 'Bahamas', 'iso3'] <- 'BHS'
cancer_fill <- cancer_fill[cancer_fill$Country.Name != 'Saint Lucia', ]
cancer_fill <- cancer_fill[!is.na(cancer_fill$iso3), c("Rate.per.person", "iso3")]
cancer_join2 <- merge(cancer_fill, country_codes, by='iso3', all.x=TRUE)
cancer_join2 <- cancer_join2[, c('Rate.per.person', 'iso3', 'id')]

cancer_corr <- rbind(cancer_join1, cancer_join2)
cancer_full <- merge(cancer_corr, country_codes, all.y=TRUE)
cancer_full[is.na(cancer_full$Rate.per.person), 'Rate.per.person'] <- 0.000242  # use background rate for countries with missing data
cancer_rate_by_country_id <- cancer_full[, c('iso3', 'id', 'Rate.per.person')]
write.csv(cancer_rate_by_country_id, "F:/NCI_NDR/Data endpoints/cancer_rate_by_country_id.csv",
          row.names=FALSE)

# resample population raster to 5min resolution with block statistics
library(raster)
target_pixel_size <- 0.08333333333333332871  # 5 arc min
population_raster <- "F:/NCI_NDR/Data population Landscan/LandScan Global 2018/LandScan2018_WGS84.tif"
population_5min <- "F:/NCI_NDR/Data population Landscan/LandScan Global 2018/LandScan2018_WGS84_5min.tif"
in_ras <- raster(population_raster)
in_pixel_size = xres(in_ras)
aggregate_factor = round(target_pixel_size / in_pixel_size)
agg_ras <- aggregate(
  in_ras, fact=aggregate_factor, fun=sum, expand=TRUE, na.rm=TRUE)
out_ras <- reclassify(agg_ras, cbind(NA, -1.))
writeRaster(out_ras, filename=population_5min, NAflag=-1.)
