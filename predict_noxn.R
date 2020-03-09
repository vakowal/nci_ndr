# Make predictions from random forest
library(ranger)
library(raster)

# noxn and covariate observations for surface water
NOXN_PREDICTOR_SURF_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_surf_2000_2015.csv"

# train the random forest model
train_df = read.csv(NOXN_PREDICTOR_SURF_DF_PATH)
ranger_train_df <- train_df[complete.cases(train_df), ]
ranger_model <- ranger(noxn~., keep.inbag=TRUE, data=ranger_train_df,
                       mtry=2, min.node.size=5)  # using mtry and min.node.size tuned by caret

# surface covariate rasters
surface_covar_path_list <- list(
  'n_export'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/sum_aggregate_to_0.084100_n_export_baseline_napp_rate_global_md5_b210146a5156422041eb7128c147512f.tif",
  'precip_variability'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/wc2.0_bio_5m_15.tif",
  'population'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/inhabitants_avg_1990_2015.tif",
  'pigs'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/5_Pg_2010_Da.tif",
  'cattle'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/5_Ct_2010_Da.tif",
  'flash_flow'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/mean_div_range_1990_2015.tif",
  'average_flow'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/average_flow_1990_2015.tif",
  'percent_no_sanitation'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/no_sanitation_provision_avg_2000-2015.tif",
  'proportion_urban'="C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/intermediate/aligned_covariates_surface/perc_urban_5min.tif"
)
covar_stack <- stack(lapply(surface_covar_path_list, FUN=raster))
names(covar_stack) <- names(surface_covar_path_list)

# Function to use ranger.predict onto a raster
predfun <- function(covar_stack, model, save_as)
{
  out <- raster(covar_stack, layer=0)
  bs <- blockSize(out)
  out <- writeStart(out, save_as, overwrite=TRUE)
  for (i in 1:bs$n) {
    covar_df <- as.data.frame(getValues(covar_stack, row=bs$row[i], nrows=bs$nrows[i]))
    na_vals <- apply(covar_df, 1, function(x) sum(is.na(x)))
    p <- numeric(length=nrow(covar_df))
    p[na_vals > 0] <- NA
    p[na_vals==0] <- predict(model, covar_df[na_vals==0, ], type='se')$se
    out <- writeValues(out, p, bs$row[i])
  }
  out <- writeStop(out)
  return(out)
}

# Run predictions
pred <- predfun(covar_stack, model=ranger_model, save_as="C:/Users/ginge/Desktop/ranger_surface_noxn_se.tif")

