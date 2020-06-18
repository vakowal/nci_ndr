# process groundwater data for India, from Esha
library(readstata13)
library(ggplot2)

cgwb_path <- "F:/NCI_NDR/Data worldbank/India_groundwater/cgwb_groundwater.dta"
wells_filled_path <- "F:/NCI_NDR/Data worldbank/India_groundwater/wq-wells-filled.dta"
cgwb_df <- read.dta13(cgwb_path)

# how many stations have valid lat/long values and nitrate?
valid_cgwb <- cgwb_df[!is.na(cgwb_df$latitude_final) &
                        !is.na(cgwb_df$longitude_final) &
                        !is.na(cgwb_df$nitrate), ]
# any duplicated station records within years?
stn_year_summary <- aggregate(nitrate~year+stncode, data=valid_cgwb, FUN=length)
dim(stn_year_summary[stn_year_summary$nitrate > 1, ])  # 16 cases where a station/year
                                                       # combination has >1 observation
p <- ggplot(valid_cgwb, aes(x=nitrate)) + geom_histogram(bins=100)
pngname <- "F:/NCI_NDR/Data worldbank/India_groundwater/cgwb_nitrate_hist.png"
png(file=pngname, units="in", res=300, width=3.5, height=4)
print(p)
dev.off()

wells_filled_df <- read.dta13(wells_filled_path)
cols_to_keep <- c('source', 'year', 'stncode', 'latitude_final', 'longitude_final',
                  'nitrateN_mean')
wells_filled_ndf <- wells_filled_df[, cols_to_keep]
wells_filled_valid <- wells_filled_ndf[!is.na(wells_filled_ndf$latitude_final) &
                                         !is.na(wells_filled_ndf$longitude_final) &
                                         !is.na(wells_filled_ndf$nitrateN_mean), ]
summary(aggregate(year~stncode, data=wells_filled_valid, FUN=min)$year)  # from 2003 to 2010
summary(aggregate(year~stncode, data=wells_filled_valid, FUN=max)$year)  # through 2016

# find non-duplicated observations, removing Esha's filled time series
wells_not_filled <- wells_filled_valid[!(duplicated(wells_filled_valid[c("stncode", "nitrateN_mean")])), ]

p <- ggplot(wells_not_filled, aes(x=nitrateN_mean)) + geom_histogram(bins=100)
pngname <- "F:/NCI_NDR/Data worldbank/India_groundwater/wells_filled_nitrateN_hist.png"
png(file=pngname, units="in", res=300, width=3.5, height=4)
print(p)
dev.off()

wells_noxn_obs <- wells_not_filled[, c('year', 'stncode', "nitrateN_mean")]
write.csv(wells_data_to_save, "F:/NCI_NDR/Data worldbank/India_groundwater/wells_noxn_obs.csv")
wells_coords <- wells_not_filled[!(duplicated(wells_not_filled$stncode)),
                                 c('stncode', "latitude_final", "longitude_final")]
wells_coords$OBJECTID <- seq(1, NROW(wells_coords), by=1)
write.csv(wells_coords, "F:/NCI_NDR/Data worldbank/India_groundwater/wells_coordinates.csv",
          row.names=FALSE)
wells_match_table <- wells_coords[, c('stncode', 'OBJECTID')]
write.csv(wells_match_table, "F:/NCI_NDR/Data worldbank/India_groundwater/wells_match_table.csv",
          row.names=FALSE)
# mean nitrate across years within station
well_noxn_by_stn <- aggregate(nitrateN_mean~stncode, data=wells_not_filled, FUN=mean)
write.csv(well_noxn_by_stn, "F:/NCI_NDR/Data worldbank/India_groundwater/mean_noxn_by_stn.csv",
          row.names=FALSE)
