"""Process global datasets related to Nitrate.

The purpose here is to aggregate spatial data from various global datasets into
scalar predictor variables to be included in an analysis of the relationship
between catchment characteristics like landcover, topography, and
precipitation, and nitrite/nitrate loadings in surface water.

Major tasks:
 - summarize covariate datasets within catchment basins

"""
import os

from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import pygeoprocessing
import pygeoprocessing.routing
import pandas
# import taskgraph

# shapefile containing locations of stations with NOxN observations
_STATION_SHP_PATH = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/WB_surface_stations_noxn_obs_sa_bas.shp"

# shapefile containing watershed area corresponding to WB stations
_WB_BASIN_SHP_PATH = "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v2/WB_surface_stations_noxn_for_snapping_ContribArea.shp"

# outer data directory
_DATA_DIR = "F:/NCI_NDR"


def map_FID_to_field(shp_path, field):
    """Map FID of each feature, according to GetFID(), to the given field.

    This allows for mapping of a dictionary of zonal statistics, where keys
    correspond to FID according to GetFID(), to another field that is preferred
    to identify features.

    Parameters:
        shp_path (string): path to shapefile
        field (string): the field to map to FID

    Returns:
        dictionary indexed by the FID of each feature retrieved with GetFID(),
            and values are the value of `field` for the feature

    """
    vector = gdal.OpenEx(shp_path, gdal.OF_VECTOR)
    layer = vector.GetLayer()
    FID_to_field = {
        feature.GetFID(): feature.GetField(field) for feature in layer}

    # clean up
    vector = None
    layer = None
    return FID_to_field


def zonal_stats_by_objectid(raster_path, band):
    """Calculate zonal stats inside watersheds by OBJECTID.

    Use the pygeoprocessing zonal_statistics() function to calculate zonal
    statistics from the given raster inside watershed features in the shapefile
    _WB_BASIN_SHP_PATH. Re-map the zonal statistics to be indexed by the field
    OBJECTID.

    Parameters:
        raster_path (string): path to the base raster to analyze with zonal
            stats
        band (int): band index of the raster to analyze

    Returns:
        nested dictionary indexed by aggregating feature OBJECTID, and then by
            one of 'min' 'max' 'sum' 'count' and 'nodata_count'.  Example:
            {0: {'min': 0, 'max': 1, 'sum': 1.7, count': 3, 'nodata_count': 1}}

    """
    fid_to_objectid = map_FID_to_field(_WB_BASIN_SHP_PATH, "OBJECTID")
    zonal_stats_dict = pygeoprocessing.geoprocessing.zonal_statistics(
        (raster_path, band), _WB_BASIN_SHP_PATH)
    objectid_zonal_stats_dict = {
        objectid: zonal_stats_dict[fid] for (fid, objectid) in
        fid_to_objectid.items()
    }
    return objectid_zonal_stats_dict


def summarize_n_fert_in_watersheds():
    """Summarize N fertilizer application rates in WB watersheds over time."""
    N_raster_path_pattern = "F:/NCI_NDR/Data fertilizer Lu Tian/Lu-Tian_2017/Nfer_ASCII/nfery<yyyy>.asc"
    year_sequence = range(1900, 2014)
    df_list = []
    for year in year_sequence:
        fert_raster_path = N_raster_path_pattern.replace('<yyyy>', str(year))
        zonal_stats_dict = zonal_stats_by_objectid(fert_raster_path, 1)
        zonal_stats_df = pandas.DataFrame.from_dict(
            zonal_stats_dict, orient='index')
        zonal_stats_df.index.name = 'OBJECTID'
        zonal_stats_df['year'] = year
        df_list.append(zonal_stats_df)
    combined_df = pandas.concat(df_list)
    save_dir = "F:/NCI_NDR/Data fertilizer Lu Tian"
    save_as = os.path.join(
        save_dir,
        "N_by_OBJECTID_WB_surface_stations_noxn_for_snapping_ContribArea.csv")
    combined_df.to_csv(save_as)


def main():
    """Program entry point."""
    summarize_n_fert_in_watersheds()


if __name__ == '__main__':
    main()
