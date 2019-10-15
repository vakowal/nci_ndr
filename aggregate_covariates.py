"""Process global datasets related to Nitrate.

The purpose here is to aggregate spatial data from various global datasets into
scalar predictor variables to be included in an analysis of the relationship
between catchment characteristics like landcover, topography, and
precipitation, and nitrite/nitrate loadings in surface water.

Major tasks:
 - summarize covariate datasets within catchment basins

"""
import os
import collections

from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import pygeoprocessing
import pygeoprocessing.routing
import pandas
import taskgraph

# shapefile containing locations of stations with NOxN observations
_STATION_SHP_PATH = "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_for_snapping_snapped.shp"

# shapefile containing watershed area corresponding to WB stations
_WB_BASIN_SHP_PATH = "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_for_snapping_ContribArea.shp"

# outer data directory
_DATA_DIR = "F:/NCI_NDR"

# covariate datasets
_COVARIATE_PATH_DICT = {
    'n_export': "F:/NCI_NDR/Data NDR/nutrient_deficit_10s_cur_compressed_md5_031d4bb444325835315a2cc825be3fd4.tif",
    'average_flow': "F:/NCI_NDR/Data streamflow FLO1K/average_flow_1990_2015.tif",
    'min_flow': "F:/NCI_NDR/Data streamflow FLO1K/average_max_flow_1990_2015.tif",
    'max_flow': "F:/NCI_NDR/Data streamflow FLO1K/average_min_flow_1990_2015.tif",
    'min_div_avg_flow': "F:/NCI_NDR/Data streamflow FLO1K/min_div_average_flow_1990_2015.tif",
    'precip_variability': "F:/NCI_NDR/Data precip Worldclim/wc2.0_bio_5m_15.tif",
    'climate_zones': "F:/NCI_NDR/Data climate zones Koeppen-Geiger/5min_updated/Map_KG-Global/KG_1986-2010.tif",
    'irrigated_area': "F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri_avg_1990_2015.tif",
    'area': "F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri_pixel_area_km2.tif",
    'population': "F:/NCI_NDR/Data population HYDE3.2/inhabitants_avg_1990_2015.tif",
    'urban_extent': "F:/NCI_NDR/Data urban extent GRUMP/glurextents.asc",
    'sanitation': "TODO",
}


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


def extract_raster_values_at_points(point_shp_path, raster_path, band):
    """Collect values from a raster intersecting points in a shapefile.

    Parameters:
        point_shp_path (string): path to shapefile containing point features
            where raster values should be extracted. Must be in geographic
            coordinates
        raster_path (string): path to raster containing values that should be
            extracted at points
        band (int): band index of the raster to analyze

    Returns:
        a data frame with one column 'OBJECTID' containing OBJECTID values of
            point features, and one column 'raster_value' containing values
            from the raster at the point location

    """
    point_vector = ogr.Open(_STATION_SHP_PATH)
    point_layer = point_vector.GetLayer()
    point_defn = point_layer.GetLayerDefn()

    # build up a list of the original field names so we can copy it to report
    point_field_name_list = ['OBJECTID']

    # this will hold (x,y) coordinates for each point in its iterator order
    point_coord_list = []
    # this maps fieldnames to a list of the values associated with that
    # fieldname in the order that the points are read in and written to
    # `point_coord_list`.
    feature_attributes_fieldname_map = collections.defaultdict(list)
    for point_feature in point_layer:
        sample_point_geometry = point_feature.GetGeometryRef()
        for field_name in point_field_name_list:
            feature_attributes_fieldname_map[field_name].append(
                point_feature.GetField(field_name))
        point_coord_list.append(
            (sample_point_geometry.GetX(), sample_point_geometry.GetY()))
    point_layer = None
    point_vector = None

    # each element will hold the point samples for each raster in the order of
    # `point_coord_list`
    sampled_precip_data_list = []
    for field_name in point_field_name_list:
        sampled_precip_data_list.append(
            pandas.Series(
                data=feature_attributes_fieldname_map[field_name],
                name=field_name))
    raster = gdal.Open(raster_path)
    band = raster.GetRasterBand(band)
    geotransform = raster.GetGeoTransform()
    sample_list = []
    for point_x, point_y in point_coord_list:
        raster_x = int((
            point_x - geotransform[0]) / geotransform[1])
        raster_y = int((
            point_y - geotransform[3]) / geotransform[5])
        sample_list.append(
            band.ReadAsArray(raster_x, raster_y, 1, 1)[0, 0])
    sampled_precip_data_list.append(
        pandas.Series(data=sample_list, name='raster_value'))

    raster = None
    band = None

    report_table = pandas.DataFrame(data=sampled_precip_data_list)
    report_table = report_table.transpose()
    return report_table


def aggregate_covariates():
    """Calculate covariate values for each surface monitoring station."""
    temp_dir = tempfile.mkdtemp()
    temp_val_dict = {
        'urban_extent_rc': os.path.join(temp_dir, 'urban_extent_rc.tif'),
    }

    # collect the following covariate values for each station:
    # N export from NDR, inside catchment area
    n_export_df = zonal_stats_by_objectid(
        _COVARIATE_PATH_DICT['n_export'], 1)

    # streamflow at the station location
    avg_flow_df = extract_raster_values_at_points(
        _STATION_SHP_PATH, _COVARIATE_PATH_DICT['average_flow'], 1)
    avg_min_flow_df = extract_raster_values_at_points(
        _STATION_SHP_PATH, _COVARIATE_PATH_DICT['min_flow'], 1)
    avg_max_flow_df = extract_raster_values_at_points(
        _STATION_SHP_PATH, _COVARIATE_PATH_DICT['max_flow'], 1)
    avg_min_div_avg_flow_df = extract_raster_values_at_points(
        _STATION_SHP_PATH, _COVARIATE_PATH_DICT['min_div_avg_flow'], 1)

    # variability of precipitation inside the catchment area
    precip_variability = zonal_stats_by_objectid(
        _COVARIATE_PATH_DICT['precip_variability'], 1)
    # calculate mean precip variability from zonal stats: sum / count

    # climate zone of the station location
    climate_zone_df = extract_raster_values_at_points(
        _STATION_SHP_PATH, _COVARIATE_PATH_DICT['climate_zones'], 1)

    # area of irrigated crops inside the catchment (square km)
    irrigated_area_df = zonal_stats_by_objectid(
        _COVARIATE_PATH_DICT['irrigated_area'], 1)

    # area of each catchment, from HYDE irrigation rasters
    catchment_area_df = zonal_stats_by_objectid(
        _COVARIATE_PATH_DICT['area'], 1)

    # calculate proportion of catchment that is irrigated, by area, by dividing
    # irrigated area by catchment area

    # population density (inhabitants per square km)
    population_density_df = zonal_stats_by_objectid(
        _COVARIATE_PATH_DICT['population'], 1)

    # reclassify urban extent df from {1 = rural and 2 = urban} to
    # {0 = rural and 1 = urban}
    urban_raster_info = pygeoprocessing.geoprocessing.get_raster_info(
        _COVARIATE_PATH_DICT['urban_extent'])
    urban_datatype = urban_raster_info['datatype']
    urban_nodata = urban_raster_info['nodata']
    value_map = {
        1: 0,
        2: 1,
    }
    pygeoprocessing.geoprocessing.reclassify_raster(
        _COVARIATE_PATH_DICT['urban_extent'], value_map,
        temp_val_dict['urban_extent_rc'], urban_datatype, urban_nodata)
    urban_extent_df = zonal_stats_by_objectid(
        temp_val_dict['urban_extent_rc'], 1)
    # process zonal stat df: percent urban area is sum / count


def main():
    """Program entry point."""
    aggregate_covariates()


if __name__ == '__main__':
    main()
