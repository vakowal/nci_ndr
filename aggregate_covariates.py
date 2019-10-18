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
import tempfile
import shutil

from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import pygeoprocessing
import pandas
import taskgraph

# shapefile containing locations of stations with NOxN observations
_STATION_SHP_PATH = "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_for_snapping_snapped.shp"

# shapefile containing watersheds
_BASIN_SHP_PATH = "F:/NCI_NDR/Data Hydrosheds/basin outlines/world_bas_15s_beta.shp"
# "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_for_snapping_ContribArea.shp"

# shapefile containing watershed centroids
_BASIN_CENTROID_PATH = "F:/NCI_NDR/Data Hydrosheds/basin outlines/centroid_world_bas_15s_beta.shp"

# outer data directory
_DATA_DIR = "F:/NCI_NDR"

# covariate datasets
_COVARIATE_PATH_DICT = {
    'n_export': "F:/NCI_NDR/Data NDR/nutrient_deficit_5min_cur_compressed_md5_031d4bb444325835315a2cc825be3fd4.tif",
    'average_flow': "F:/NCI_NDR/Data streamflow FLO1K/average_flow_1990_2015.tif",
    'min_flow': "F:/NCI_NDR/Data streamflow FLO1K/average_max_flow_1990_2015.tif",
    'max_flow': "F:/NCI_NDR/Data streamflow FLO1K/average_min_flow_1990_2015.tif",
    'min_div_avg_flow': "F:/NCI_NDR/Data streamflow FLO1K/min_div_average_flow_1990_2015.tif",
    'precip_variability': "F:/NCI_NDR/Data precip Worldclim/wc2.0_bio_5m_15.tif",
    'climate_zones': "F:/NCI_NDR/Data climate zones Koeppen-Geiger/5min_updated/Map_KG-Global/KG_1986-2010.tif",
    'irrigated_area': "F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri_avg_1990_2015.tif",
    'area': "F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri_pixel_area_km2.tif",
    'population': "F:/NCI_NDR/Data population HYDE3.2/inhabitants_avg_1990_2015.tif",
    'urban_extent_rc': "F:/NCI_NDR/Data urban extent GRUMP/glurextents_rc.tif",
    'sanitation_table': "F:/NCI_NDR/Data sanitation/no_sanitation_provision_avg_2000-2015.csv",
    'countries_raster': "F:/NCI_NDR/Data world borders/TM_WORLD_BORDERS-03_countryid.tif",
}


def zonal_sum_to_csv(input_raster_path, sum_field_name, save_as):
    """Calculate zonal sum and save the result as a data frame in a csv file.

    Use zonal statistics to calculate the sum of pixel values from
    `input_raster` lying inside watershed features. Save the zonal sum as a csv
    file called `save_as`, in a data frame with two columns: OBJECTID to
    identify watershed features, and `sum_field_name` to identify the zonal
    sum.

    Parameters:
        input_raster_path (string): path to the raster containing values that
            should be summarized by zonal sum
        sum_field_name (string): field name in the csv data frame for the
            column containing zonal sum
        save_as (string): path to location on disk where the csv data frame
            should be saved

    Returns:
        None

    """
    print("Calculating zonal sum from {}".format(input_raster_path))
    zonal_df = zonal_stats_by_objectid(input_raster_path, 1)
    zonal_df.rename(columns={'sum': sum_field_name}, inplace=True)
    zonal_df.drop(
        ['min', 'max', 'count', 'nodata_count'], axis='columns', inplace=True)
    zonal_df.to_csv(save_as, index=False)


def zonal_max_value_to_csv(input_raster_path, field_name, save_as):
    """Save zonal max value inside watersheds as a data frame in a csv file.

    Use zonal statistics to calculate the maximum pixel value from
    `input_raster` lying inside watershed features. Save the zonal max as a csv
    file called `save_as`, in a data frame with two columns: OBJECTID to
    identify watershed features, and `field_name` to identify the zonal
    max.

    Parameters:
        input_raster_path (string): path to the raster containing values that
            should be summarized by zonal max
        field_name (string): field name in the csv data frame for the
            column containing zonal max
        save_as (string): path to location on disk where the csv data frame
            should be saved

    Returns:
        None

    """
    print("Calculating zonal max from {}".format(input_raster_path))
    zonal_df = zonal_stats_by_objectid(input_raster_path, 1)
    zonal_df.rename(columns={'max': field_name}, inplace=True)
    zonal_df.drop(
        ['min', 'sum', 'count', 'nodata_count'], axis='columns', inplace=True)
    zonal_df.to_csv(save_as, index=False)


def aggregate_precip_variability(save_as):
    """Calculate average precipitation variability inside watersheds."""
    precip_var_df = zonal_stats_by_objectid(
        _COVARIATE_PATH_DICT['precip_variability'], 1)
    precip_var_df['mean_precip_variability'] = (
        precip_var_df['sum'] / precip_var_df['count'])
    precip_var_df.drop(
        ['min', 'max', 'count', 'nodata_count', 'sum'], axis='columns',
        inplace=True)
    precip_var_df.to_csv(save_as, index=False)


def aggregate_proportion_irrigation(area_df_path, save_as):
    """Calculate the proportion of area inside the basin that is irrigated."""
    irrigated_area_df = zonal_stats_by_objectid(
        _COVARIATE_PATH_DICT['irrigated_area'], 1)
    irrigated_area_df.rename(
        columns={'sum': 'irrigated_area_sum'}, inplace=True)
    total_area_df = pandas.read_csv(area_df_path)
    irrigated_area_df = irrigated_area_df.merge(
        total_area_df, on='OBJECTID', suffixes=(False, False))
    irrigated_area_df['proportion_irrigated_area'] = (
        irrigated_area_df['irrigated_area_sum'] /
        irrigated_area_df['basin_sum_area_km2'])
    irrigated_area_df.drop([
        'min', 'max', 'count', 'nodata_count', 'irrigated_area_sum',
        'basin_sum_area_km2'], axis='columns', inplace=True)
    irrigated_area_df.to_csv(save_as, index=False)


def reclassify_urban_extent(save_as):
    """Reclassify urban extent raster.

    Reclassify from {1=rural and 2=urban} to {0 = rural and 1 = urban}.

    Parameters:
        save_as (string): path to save the reclassified raster

    Returns:
        None

    """
    Print("Reclassifying urban extent ...")
    urban_raster_info = pygeoprocessing.get_raster_info(
        _COVARIATE_PATH_DICT['urban_extent'])
    urban_datatype = urban_raster_info['datatype']
    urban_nodata = -9999  # weirdly, these two are not compatible
    value_map = {
        1: 0,
        2: 1,
    }
    pygeoprocessing.reclassify_raster(
        (_COVARIATE_PATH_DICT['urban_extent'], 1), value_map, save_as,
        urban_datatype, urban_nodata)


def aggregate_proportion_urban(urban_extent_rc_path, save_as):
    """Calculate the proportion of the watershed that is urban."""
    Print("Aggregating proportion urban ...")
    urban_extent_df = zonal_stats_by_objectid(
        urban_extent_rc_path, 1)
    urban_extent_df['proportion_urban'] = (
        urban_extent_df['sum'] / urban_extent_df['count'])
    urban_extent_df.drop(
        ['min', 'max', 'count', 'nodata_count', 'sum'], axis='columns',
        inplace=True)
    urban_extent_df.to_csv(save_as, index=False)


def reclassify_countries_by_sanitation(save_as):
    """Reclassify countries raster by sanitation provision per country."""
    sanitation_df = pandas.read_csv(_COVARIATE_PATH_DICT['sanitation_table'])
    countryid_to_sanitation = pandas.Series(
        sanitation_df.no_san_provision.values,
        index=sanitation_df.countryid).to_dict()
    pygeoprocessing.reclassify_raster(
        (_COVARIATE_PATH_DICT['countries_raster'], 1),
        countryid_to_sanitation, save_as, gdal.GDT_Float32, -9999.)


def aggregate_proportion_no_sanitation(sanitation_rc_path, save_as):
    """Aggregate sanitation.

    Calculate the proportion of inhabitants in the watershed without access to
    sanitation.

    Parameters:
        sanitation_rc_path (string): path to raster containing proportion of
            population without access to sanitation, by country. Each pixel
            gives the country-level value of proportion of the population
            without sanitation.
        save_as (string): path to location on disk to save the result data
            frame

    Returns:
        None

    """
    sanitation_df = zonal_stats_by_objectid(sanitation_rc_path, 1)
    sanitation_df['proportion_no_sanitation'] = (
        sanitation_df['sum'] / sanitation_df['count'])
    sanitation_df.drop(
        ['min', 'max', 'count', 'nodata_count', 'sum'], axis='columns',
        inplace=True)
    sanitation_df.to_csv(save_as, index=False)


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
    _BASIN_SHP_PATH. Re-map the zonal statistics to be indexed by the field
    OBJECTID.  Convert the zonal stats nested dictionary to a pandas dataframe.

    Parameters:
        raster_path (string): path to the base raster to analyze with zonal
            stats
        band (int): band index of the raster to analyze

    Returns:
        data frame where the index is OBJECTID of the watershed layer,
            containing the columns 'min', 'max', 'count', 'nodata_count',
            and 'sum'

    """
    fid_to_objectid = map_FID_to_field(_BASIN_SHP_PATH, "OBJECTID")
    zonal_stats_dict = pygeoprocessing.zonal_statistics(
        (raster_path, band), _BASIN_SHP_PATH)
    objectid_zonal_stats_dict = {
        objectid: zonal_stats_dict[fid] for (fid, objectid) in
        fid_to_objectid.items()
    }
    objectid_df = pandas.DataFrame(objectid_zonal_stats_dict)
    objectid_df_t = objectid_df.transpose()
    objectid_df_t['OBJECTID'] = objectid_df_t.index
    return objectid_df_t


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
    combined_df.to_csv(save_as, index=False)


def raster_values_at_points(
        point_shp_path, raster_path, band, raster_field_name, save_as):
    """Collect values from a raster intersecting points in a shapefile.

    Create
    Parameters:
        point_shp_path (string): path to shapefile containing point features
            where raster values should be extracted. Must be in geographic
            coordinates
        raster_path (string): path to raster containing values that should be
            extracted at points
        band (int): band index of the raster to analyze
        raster_field_name (string): name to assign to the field in the data
            frame that contains values extracted from the raster
        save_as (string): path to location to save the data frame

    Side effects:
        creates or modifies the csv file indicated by `save_as`:
        a data frame with one column 'OBJECTID' containing OBJECTID values of
            point features, and one column raster_field_name containing values
            from the raster at the point location

    Returns:
        None

    """
    point_vector = ogr.Open(point_shp_path)
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
        pandas.Series(data=sample_list, name=raster_field_name))

    raster = None
    band = None

    report_table = pandas.DataFrame(data=sampled_precip_data_list)
    report_table = report_table.transpose()
    report_table.to_csv(save_as, index=False)


def merge_data_frame_list(df_path_list, save_as):
    """Merge the data frames in `df_path_list` and save as one data frame.

    Merge each data frame in `df_path_list` into a single data frame. Save this
    merged data frame as `save_as`.

    Parameters:
        df_path_list (list): list of file paths indicating locations of data
            frames that should be merged. Each must include a column 'OBJECTID'
            identifying unique watersheds or monitoring stations, and a column
            of covariate data
        save_as (string): path to location on disk where the result should be
            saved

    Returns:
        None

    """
    combined_df = pandas.read_csv(df_path_list[0])
    df_i = 1
    while df_i < len(df_path_list):
        combined_df = combined_df.merge(
            pandas.read_csv(df_path_list[df_i]), on='OBJECTID',
            suffixes=(False, False), validate="one_to_one")
        df_i = df_i + 1
    combined_df.to_csv(save_as, index=False)


def aggregate_covariates(intermediate_dir_path, combined_covariate_table_path):
    """Calculate covariate values for each surface monitoring station.

    Parameters:
        intermediate_dir_path (string): path to folder where persistent
            intermediate outputs should be written, namely a csv for each
            covariate
        combined_covariate_table_path (string): path to location where the
            combined covariate table should be written, i.e. the table
            containing all covariate values

    """
    # list of paths to data frames each containing aggregated values for one
    # covariate
    df_path_list = []

    temp_dir = tempfile.mkdtemp()
    temp_val_dict = {
        'urban_extent_rc': os.path.join(temp_dir, 'urban_extent_rc.tif'),
        'sanitation': os.path.join(temp_dir, 'sanitation.tif'),
    }

    # avg_flow_at_point_path = os.path.join(
    #     intermediate_dir_path, 'avg_flow_point.csv')
    # df_path_list.append(avg_flow_at_point_path)
    # if not os.path.exists(avg_flow_at_point_path):
    #     raster_values_at_points(
    #         _STATION_SHP_PATH, _COVARIATE_PATH_DICT['average_flow'], 1,
    #         'average_flow_at_point', avg_flow_at_point_path)

    avg_flow_zonal_max_path = os.path.join(
        intermediate_dir_path, 'avg_flow_zonal_max.csv')
    df_path_list.append(avg_flow_zonal_max_path)
    if not os.path.exists(avg_flow_zonal_max_path):
        zonal_max_value_to_csv(
            _COVARIATE_PATH_DICT['average_flow'], 'average_flow_zonal_max',
            avg_flow_zonal_max_path)

    # avg_min_flow_path = os.path.join(intermediate_dir_path, 'avg_min_flow.csv')
    # df_path_list.append(avg_min_flow_path)
    # if not os.path.exists(avg_min_flow_path):
    #     raster_values_at_points(
    #         _STATION_SHP_PATH, _COVARIATE_PATH_DICT['min_flow'], 1,
    #         'average_min_flow', avg_min_flow_path)

    # avg_max_flow_path = os.path.join(intermediate_dir_path, 'avg_max_flow.csv')
    # df_path_list.append(avg_max_flow_path)
    # if not os.path.exists(avg_max_flow_path):
    #     raster_values_at_points(
    #         _STATION_SHP_PATH, _COVARIATE_PATH_DICT['max_flow'], 1,
    #         'average_max_flow', avg_max_flow_path)

    # avg_min_div_avg_flow_path = os.path.join(
    #     intermediate_dir_path, 'avg_min_div_avg_flow.csv')
    # df_path_list.append(avg_min_div_avg_flow_path)
    # if not os.path.exists(avg_min_div_avg_flow_path):
    #     raster_values_at_points(
    #         _STATION_SHP_PATH, _COVARIATE_PATH_DICT['min_div_avg_flow'], 1,
    #         'average_min_div_avg_flow', avg_min_div_avg_flow_path)

    climate_zone_df_path = os.path.join(
        intermediate_dir_path, 'climate_zone.csv')
    df_path_list.append(climate_zone_df_path)
    if not os.path.exists(climate_zone_df_path):
        raster_values_at_points(
            _BASIN_CENTROID_PATH, _COVARIATE_PATH_DICT['climate_zones'], 1,
            'climate_zone', climate_zone_df_path)

    population_df_path = os.path.join(intermediate_dir_path, 'population.csv')
    df_path_list.append(population_df_path)
    if not os.path.exists(population_df_path):
        zonal_sum_to_csv(
            _COVARIATE_PATH_DICT['population'], 'basin_sum_population',
            population_df_path)

    proportion_urban_df_path = os.path.join(
        intermediate_dir_path, 'proportion_urban.csv')
    df_path_list.append(proportion_urban_df_path)
    if not os.path.exists(proportion_urban_df_path):
        aggregate_proportion_urban(
            _COVARIATE_PATH_DICT['urban_extent_rc'], proportion_urban_df_path)

    precip_var_df_path = os.path.join(
        intermediate_dir_path, 'precip_variability.csv')
    df_path_list.append(precip_var_df_path)
    if not os.path.exists(precip_var_df_path):
        aggregate_precip_variability(precip_var_df_path)

    basin_area_df_path = os.path.join(
        intermediate_dir_path, 'basin_area.csv')
    df_path_list.append(basin_area_df_path)
    if not os.path.exists(basin_area_df_path):
        zonal_sum_to_csv(
            _COVARIATE_PATH_DICT['area'], 'basin_sum_area_km2',
            basin_area_df_path)

    proportion_irrigation_df_path = os.path.join(
        intermediate_dir_path, 'proportion_irrigated_area.csv')
    df_path_list.append(proportion_irrigation_df_path)
    if not os.path.exists(proportion_irrigation_df_path):
        aggregate_proportion_irrigation(
            basin_area_df_path, proportion_irrigation_df_path)

    proportion_no_sanitation_df_path = os.path.join(
        intermediate_dir_path, 'proportion_no_sanitation.csv')
    df_path_list.append(proportion_no_sanitation_df_path)
    if not os.path.exists(proportion_no_sanitation_df_path):
        reclassify_countries_by_sanitation(temp_val_dict['sanitation'])
        aggregate_proportion_no_sanitation(
            temp_val_dict['sanitation'], proportion_no_sanitation_df_path)

    ndr_df_path = os.path.join(intermediate_dir_path, 'n_export.csv')
    df_path_list.append(ndr_df_path)
    if not os.path.exists(ndr_df_path):
        zonal_sum_to_csv(
            _COVARIATE_PATH_DICT['n_export'], 'basin_sum_n_export',
            ndr_df_path)

    merge_data_frame_list(df_path_list, combined_covariate_table_path)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def aggregate_covariates_taskgraph_mode(
            intermediate_dir_path, combined_covariate_table_path):
    """Aggregate covariates, with taskgraph for parallel processing."""
    df_path_list = []

    temp_dir = tempfile.mkdtemp()
    temp_val_dict = {
        'urban_extent_rc': os.path.join(temp_dir, 'urban_extent_rc.tif'),
        'sanitation': os.path.join(temp_dir, 'sanitation.tif'),
    }

    work_token_dir = os.path.join(intermediate_dir_path, '_work_tokens')
    graph = taskgraph.TaskGraph(work_token_dir, n_workers=4)  # n_workers=4)

    ndr_df_path = os.path.join(intermediate_dir_path, 'n_export.csv')
    df_path_list.append(ndr_df_path)
    aggregate_ndr_task = graph.add_task(
        zonal_sum_to_csv,
        args=(
            _COVARIATE_PATH_DICT['n_export'], 'basin_sum_n_export',
            ndr_df_path),
        target_path_list=[ndr_df_path],
        task_name='aggregate_ndr_task')

    avg_flow_path = os.path.join(intermediate_dir_path, 'avg_flow.csv')
    df_path_list.append(avg_flow_path)
    # avg_flow_task = graph.add_task(
    #     raster_values_at_points,
    #     args=[
    #         _STATION_SHP_PATH, _COVARIATE_PATH_DICT['average_flow'], 1,
    #         'average_flow', avg_flow_path],
    #     target_path_list=avg_flow_path,
    #     task_name='avg_flow_task')

    avg_min_flow_path = os.path.join(intermediate_dir_path, 'avg_min_flow.csv')
    df_path_list.append(avg_min_flow_path)
    # avg_min_flow_task = graph.add_task(
    #     raster_values_at_points,
    #     args=[
    #         _STATION_SHP_PATH, _COVARIATE_PATH_DICT['min_flow'], 1,
    #         'average_min_flow', avg_min_flow_path],
    #     target_path_list=avg_min_flow_path,
    #     task_name='avg_min_flow_task')

    avg_max_flow_path = os.path.join(intermediate_dir_path, 'avg_max_flow.csv')
    df_path_list.append(avg_max_flow_path)
    # avg_max_flow_task = graph.add_task(
    #     raster_values_at_points,
    #     args=[
    #         _STATION_SHP_PATH, _COVARIATE_PATH_DICT['max_flow'], 1,
    #         'average_max_flow', avg_max_flow_path],
    #     target_path_list=avg_max_flow_path,
    #     task_name='avg_max_flow_task')

    avg_min_div_avg_flow_path = os.path.join(
        intermediate_dir_path, 'avg_min_div_avg_flow.csv')
    df_path_list.append(avg_min_div_avg_flow_path)
    # avg_min_div_avg_flow_task = graph.add_task(
    #     raster_values_at_points,
    #     args=[
    #         _STATION_SHP_PATH, _COVARIATE_PATH_DICT['min_div_avg_flow'], 1,
    #         'average_min_div_avg_flow', avg_min_div_avg_flow_path],
    #     target_path_list=avg_min_div_avg_flow_path,
    #     task_name='avg_min_div_avg_flow_task')

    precip_var_df_path = os.path.join(
        intermediate_dir_path, 'precip_variability.csv')
    df_path_list.append(precip_var_df_path)
    # precip_var_task = graph.add_task(
    #     aggregate_precip_variability,
    #     args=(precip_var_df_path),
    #     target_path_list=[precip_var_df_path],
    #     task_name='aggregate_precip_variability_task')

    climate_zone_df_path = os.path.join(
        intermediate_dir_path, 'climate_zone.csv')
    df_path_list.append(climate_zone_df_path)
    # climate_zone_task = graph.add_task(
    #     raster_values_at_points,
    #     args=[
    #         _STATION_SHP_PATH, _COVARIATE_PATH_DICT['climate_zones'], 1,
    #         'climate_zone', climate_zone_df_path],
    #     target_path_list=climate_zone_df_path,
    #     task_name='climate_zone_task')

    basin_area_df_path = os.path.join(
        intermediate_dir_path, 'basin_area.csv')
    df_path_list.append(basin_area_df_path)
    # aggregate_basin_area_task = graph.add_task(
    #     zonal_sum_to_csv,
    #     args=(
    #         _COVARIATE_PATH_DICT['area'], 'basin_sum_area_km2',
    #         basin_area_df_path),
    #     target_path_list=[basin_area_df_path],
    #     task_name='aggregate_basin_area_task')

    proportion_irrigation_df_path = os.path.join(
        intermediate_dir_path, 'proportion_irrigated_area.csv')
    df_path_list.append(proportion_irrigation_df_path)
    # irrigated_area_task = graph.add_task(
    #     aggregate_proportion_irrigation,
    #     args=(basin_area_df_path, proportion_irrigation_df_path),
    #     target_path_list=[proportion_irrigation_df_path],
    #     dependent_task_list=[aggregate_basin_area_task],
    #     task_name='proportion_irrigation_task')

    population_df_path = os.path.join(intermediate_dir_path, 'population.csv')
    df_path_list.append(population_df_path)
    # population_task = graph.add_task(
    #     zonal_sum_to_csv,
    #     args=(
    #         _COVARIATE_PATH_DICT['population'], 'basin_sum_population',
    #         population_df_path),
    #     target_path_list=[population_df_path],
    #     task_name='aggregate_population_task')

    # reclassify_urban_task = graph.add_task(
    #     reclassify_urban_extent,
    #     args=[temp_val_dict['urban_extent_rc']],
    #     target_path_list=[temp_val_dict['urban_extent_rc']],
    #     task_name='reclassify_urban_task')

    proportion_urban_df_path = os.path.join(
        intermediate_dir_path, 'proportion_urban.csv')
    df_path_list.append(proportion_urban_df_path)
    # proportion_urban_task = graph.add_task(
    #     aggregate_proportion_urban,
    #     args=(temp_val_dict['urban_extent_rc'], proportion_urban_df_path),
    #     target_path_list=[proportion_urban_df_path],
    #     dependent_task_list=[reclassify_urban_task],
    #     task_name='aggregate_proportion_urban')

    # reclassify_sanitation_task = graph.add_task(
    #     reclassify_countries_by_sanitation,
    #     args=(temp_val_dict['sanitation']),
    #     target_path_list=[temp_val_dict['sanitation']],
    #     task_name='reclassify_sanitation_task')

    proportion_no_sanitation_df_path = os.path.join(
        intermediate_dir_path, 'proportion_no_sanitation.csv')
    df_path_list.append(proportion_no_sanitation_df_path)
    # proportion_no_sanitation_task = graph.add_task(
    #     aggregate_proportion_no_sanitation,
    #     args=(temp_val_dict['sanitation'], proportion_no_sanitation_df_path),
    #     target_path_list=[proportion_no_sanitation_df_path],
    #     dependent_task_list=[reclassify_sanitation_task],
    #     task_name='aggregate_sanitation_task')

    graph.close()
    graph.join()

    merge_data_frame_list(df_path_list, combined_covariate_table_path)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def main():
    """Program entry point."""
    out_dir = 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/HydroSHEDS_watersheds'
    intermediate_dir_path = os.path.join(
        out_dir, 'intermediate_df_dir')
    combined_covariate_table_path = os.path.join(
        out_dir, 'combined_covariates.csv')
    aggregate_covariates(
        intermediate_dir_path, combined_covariate_table_path)


if __name__ == '__main__':
    main()
