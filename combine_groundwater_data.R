# Generate combined table of groundwater nitrate observations

# table containing original OBJECTID and source for each dataset
combined_match_df <- read.csv("F:/NCI_NDR/Data groundwater merged/groundwater_GEMStat_Gu_USGS_Ouedraogo.csv")

# GEMStat observations
NOXN_OBSERVATIONS_CSV <- "F:/NCI_NDR/Data worldbank/station_data/noxn_obs_all.csv"
NOXN_OBS <- read.csv(NOXN_OBSERVATIONS_CSV, stringsAsFactors=FALSE)
STN_DF <- read.csv("F:/NCI_NDR/Data worldbank/station_data/station_metadata.csv")
stn_cols <- colnames(STN_DF)[c(1, 4, 5, 8:10, 15:16)]
stn_subset <- STN_DF[, stn_cols]
stn_subset[stn_subset$Water.Type == 'Groundwater station', 'ground_v_surface'] <- 'groundwater'
NOXN_BY_STN <- merge(NOXN_OBS, stn_subset, all.x=TRUE)
NOXN_BY_STN$Sample.Date <- as.Date(NOXN_BY_STN$Sample.Date, format="%Y-%m-%d")
gr_stn <- NOXN_BY_STN[NOXN_BY_STN$ground_v_surface == 'groundwater', ]
gr_stn$year <- format(gr_stn$Sample.Date, format="%Y")
mean_noxn_by_year <- aggregate(Value~GEMS.Station.Number+year+Unit,
                               data=gr_stn, FUN=mean)
noxn_obs_restr <- mean_noxn_by_year[, c(1, 4)]
colnames(noxn_obs_restr) <- c('GEMS.Stati', 'noxn')

gems_objectid_match_df <- read.csv("F:/NCI_NDR/Data worldbank/station_data/WB_groundwater_stations_noxn_obs_objectid_match.csv")
gems_objectid_match_df <- gems_objectid_match_df[, c('GEMS.Stati', 'OBJECTID')]
colnames(gems_objectid_match_df)[2] <- 'OBJECTID_o'
gems_noxn_obs <- merge(noxn_obs_restr, gems_objectid_match_df, all.x=TRUE)
gems_noxn_obs <- gems_noxn_obs[, c('OBJECTID_o', 'noxn')]
gems_noxn_obs$source <- 'GEMStat'  # 607 rows

# Gu et al observations
gu_noxn_obs <- read.csv("F:/NCI_NDR/Data groundwater Gu et al/noxn_observations.csv")
colnames(gu_noxn_obs) <- c('OBJECTID_o', 'noxn')
gu_noxn_obs$source <- 'Gu'  # 600 rows

# USGS observations
usgs_noxn_obs <- read.csv("F:/NCI_NDR/Data groundwater USGS/Decadal_change_in_groundwater_NO3N_datav2.csv")
usgs_noxn_obs <- usgs_noxn_obs[, c('OBJECTID', 'Mean_02_13')]
colnames(usgs_noxn_obs) <- c('OBJECTID_o', 'noxn')
usgs_noxn_obs$source <- 'USGS'  # 909 rows

# Ouedraogo et al observations
oue_noxn_obs <- read.csv("F:/NCI_NDR/Data groundwater Ouedrago et al/Database_Issoufou_OUEDRAOGO_cor_mean_noxn.csv")
oue_noxn_obs <- oue_noxn_obs[, c('OBJECTID', 'Mean_NO3N')]
colnames(oue_noxn_obs) <- c('OBJECTID_o', 'noxn')
oue_noxn_obs$source <- 'Ouedraogo'  # 195 rows

# combine noxn observations
noxn_obs_combined <- do.call(rbind, list(
    gems_noxn_obs, gu_noxn_obs, usgs_noxn_obs, oue_noxn_obs))  # 2311 rows

# merge with new combined OBJECTID
objectid_merge <- merge(combined_match_df, noxn_obs_combined)
noxn_obs_to_save <- objectid_merge[, c('OBJECTID', 'source', 'noxn')]
noxn_obs_to_save <- noxn_obs_to_save[complete.cases(noxn_obs_to_save), ]
write.csv(
    noxn_obs_to_save, "F:/NCI_NDR/Data groundwater merged/combined_noxn_groundwater_GEMStat_Gu_USGS_Ouedraogo.csv",
    row.names=FALSE)
