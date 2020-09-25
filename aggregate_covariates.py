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
import numpy
import pandas
import taskgraph

# shapefile containing locations of snapped stations with NOxN observations
_SNAPPED_STATION_SHP_PATH = "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_for_snapping_snapped.shp"

# shapefile containing original station locations (not snapped)
_ORIG_STATION_SHP_PATH = "F:/NCI_NDR/Data worldbank/station_data/WB_surface_stations_noxn_obs_objectid.shp"

# shapefile containing watersheds
_BASIN_SHP_PATH = "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v3/WB_surface_stations_noxn_for_snapping_ContribArea.shp"
# "F:/NCI_NDR/Data Hydrosheds/basin outlines/world_bas_15s_beta.shp"

# shapefile containing watershed centroids
_BASIN_CENTROID_PATH = "F:/NCI_NDR/Data Hydrosheds/basin outlines/centroid_world_bas_15s_beta.shp"

_COVARIATE_PATH_DICT = None

# covariate datasets, 5 min resolution
_COVARIATE_PATH_DICT_5min = {
    'n_export': "F:/NCI_NDR/Data NDR/updated_5.18.20/sum_aggregate_to_0.084100_n_export_fixedarea_currentpractices_global.tif",
    # 'n_export': "F:/NCI_NDR/Data NDR/updated_5.18.20/n_load_inputs/n_load_natveg_ag_5min.tif",
    'average_flow': "F:/NCI_NDR/Data streamflow FLO1K/average_flow_1990_2015.tif",
    'flash_flow': "F:/NCI_NDR/Data streamflow FLO1K/mean_div_range_1990_2015.tif",
    'precip_variability': "F:/NCI_NDR/Data precip Worldclim/wc2.0_bio_5m_15.tif",
    'climate_zones': "F:/NCI_NDR/Data climate zones Koeppen-Geiger/5min_updated/Map_KG-Global/KG_1986-2010.tif",
    'irrigated_area': "F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri_avg_1990_2015.tif",
    'area': "F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri_pixel_area_km2.tif",
    'population': "F:/NCI_NDR/Data population HYDE3.2/inhabitants_avg_1990_2015.tif",
    'urban_extent_rc': "F:/NCI_NDR/Data urban extent GRUMP/perc_urban_5min.tif",
    'sanitation_table': "F:/NCI_NDR/Data sanitation/no_sanitation_provision_avg_2000-2015.csv",
    'countries_raster': "F:/NCI_NDR/Data world borders/TM_WORLD_BORDERS-03_countryid_5min.tif",
    'depth_to_groundwater': "F:/NCI_NDR/Data HydroATLAS/gwt_cm_sav_level12.tif",
    'clay_percent': "F:/NCI_NDR/Data soil ISRIC/CLYPPT_M_sl1_10km_ll.tif",
    'sand_percent': "F:/NCI_NDR/Data soil ISRIC/SNDPPT_M_sl1_10km_ll.tif",
    'cattle': "F:/NCI_NDR/Data GLW/5_Ct_2010_Da.tif",
    'pigs': "F:/NCI_NDR/Data GLW/5_Pg_2010_Da.tif",
}

# covariate datasets, 30 sec resolution when possible
_COVARIATE_PATH_DICT_30s = {
    'n_export': "F:/NCI_NDR/Data NDR/updated_5.18.20/n_export_30s_fixedarea_currentpractices_global.tif",
    'average_flow': "F:/NCI_NDR/Data streamflow FLO1K/30s_resolution/average_flow_1990_2015.tif",
    'flash_flow': "F:/NCI_NDR/Data streamflow FLO1K/30s_resolution/mean_div_range_1990_2015.tif",
    'precip_variability': "F:/NCI_NDR/Data precip Worldclim/wc2.1_30s_bio_15.tif",
    'climate_zones': "F:/NCI_NDR/Data climate zones Koeppen-Geiger/5min_updated/Map_KG-Global/KG_1986-2010.tif",
    'irrigated_area': "F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri_avg_1990_2015.tif",
    'area': "F:/NCI_NDR/Data irrigation HYDE3.2/tot_irri_pixel_area_km2.tif",
    'population': "F:/NCI_NDR/Data population HYDE3.2/inhabitants_avg_1990_2015.tif",
    'urban_extent_rc': "F:/NCI_NDR/Data urban extent GRUMP/glurextents_rc.tif",
    'sanitation_table': "F:/NCI_NDR/Data sanitation/no_sanitation_provision_avg_2000-2015.csv",
    'countries_raster': "F:/NCI_NDR/Data world borders/TM_WORLD_BORDERS-03_countryid.tif",
    'depth_to_groundwater': "F:/NCI_NDR/Data HydroATLAS/gwt_cm_sav_level12_30s.tif",
    'clay_percent': "F:/NCI_NDR/Data soil ISRIC/CLYPPT_M_sl6_1km_ll.tif",
    'sand_percent': "F:/NCI_NDR/Data soil ISRIC/SNDPPT_M_sl6_1km_ll.tif",
    'cattle': "F:/NCI_NDR/Data GLW/5_Ct_2010_Da.tif",
    'pigs': "F:/NCI_NDR/Data GLW/5_Pg_2010_Da.tif",
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


def zonal_mean_value_to_csv(input_raster_path, field_name, save_as):
    """Save zonal mean value inside watersheds as a data frame in a csv file.

    Use zonal statistics to calculate the mean pixel value from
    `input_raster` lying inside watershed features. Save the zonal mean as a
    csv file called `save_as`, in a data frame with two columns: OBJECTID to
    identify watershed features, and `field_name` to identify the zonal
    mean.

    Parameters:
        input_raster_path (string): path to the raster containing values that
            should be summarized by zonal mean
        field_name (string): field name in the csv data frame for the
            column containing zonal mean
        save_as (string): path to location on disk where the csv data frame
            should be saved

    Returns:
        None

    """
    print("Calculating zonal mean from {}".format(input_raster_path))
    zonal_df = zonal_stats_by_objectid(input_raster_path, 1)
    zonal_df[field_name] = (
        zonal_df['sum'] / zonal_df['count'])
    zonal_df.drop(
        ['min', 'max', 'count', 'nodata_count', 'sum'], axis='columns',
        inplace=True)
    zonal_df.to_csv(save_as, index=False)


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


def proportion_irrigation_at_point(point_shp_path, save_as):
    """Calculate proportion irrigated area for pixels intersecting points.

    From a raster containing pixel areas, get the area of the pixel
    intersecting each point feature; from a raster containing irrigated area,
    get the irrigated area of the pixel intersecting each point feature;
    calculate proportion irrigated area for each point feature by dividing
    irrigated area by pixel area.

    Parameters:
        point_shp_path (string): path to shapefile containing points where
            values should be calculated
        save_as (string): path to save the proportion irrigated area as csv

    Returns:
        None

    """
    temp_dir = tempfile.mkdtemp()
    # get area of pixel intersecting station features
    area_df_path = os.path.join(temp_dir, 'area_by_point.csv')
    raster_values_at_points(
        point_shp_path, _COVARIATE_PATH_DICT['area'], 1, 'pixel_area',
        area_df_path)

    # get irrigated area in pixel intersecting station features
    irrigated_area_df_path = os.path.join(
        temp_dir, 'irrigated_area_by_point.csv')
    raster_values_at_points(
        point_shp_path, _COVARIATE_PATH_DICT['irrigated_area'], 1,
        'irrigated_area', irrigated_area_df_path)

    # calculate proportion of the pixel that is irrigated by dividing
    area_df = pandas.read_csv(area_df_path)
    area_df = area_df.merge(
        pandas.read_csv(irrigated_area_df_path), on='OBJECTID',
        suffixes=(False, False), validate="one_to_one")
    area_df['proportion_irrigated_area'] = (
        area_df['irrigated_area'] / area_df['pixel_area'])
    area_df.drop([
        'irrigated_area', 'pixel_area'], axis='columns', inplace=True)
    area_df.to_csv(save_as, index=False)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def reclassify_urban_extent(save_as):
    """Reclassify urban extent raster.

    Reclassify from {1=rural and 2=urban} to {0 = rural and 1 = urban}.

    Parameters:
        save_as (string): path to save the reclassified raster

    Returns:
        None

    """
    print("Reclassifying urban extent ...")
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


def reclassify_countries_by_sanitation(countries_raster, save_as):
    """Reclassify countries raster by sanitation provision per country.

    Parameters:
        countries_raster (string): path to raster identifying countries
        save_as (string): location to save raster where country id values have
            been reclassified to proportion without sanitation provision.

    """
    sanitation_df = pandas.read_csv(_COVARIATE_PATH_DICT['sanitation_table'])
    countryid_to_sanitation = pandas.Series(
        sanitation_df.no_san_provision.values,
        index=sanitation_df.countryid).to_dict()
    pygeoprocessing.reclassify_raster(
        (countries_raster, 1), countryid_to_sanitation, save_as,
        gdal.GDT_Float32, -9999.)


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


def summarize_n_fert_at_points(point_shp_path, save_as):
    """Summarize N fertilizer application rates at pixels over time."""
    temp_dir = tempfile.mkdtemp()
    N_raster_path_pattern = "F:/NCI_NDR/Data fertilizer Lu Tian/Lu-Tian_2017/Nfer_ASCII/nfery<yyyy>.asc"
    year_sequence = range(1900, 2014)
    df_list = []
    for year in year_sequence:
        fert_raster_path = N_raster_path_pattern.replace('<yyyy>', str(year))
        intermediate_df_path = os.path.join(
            temp_dir, 'fert_{}.csv'.format(year))
        raster_values_at_points(
            point_shp_path, fert_raster_path, 1, 'fertilizer',
            intermediate_df_path)
        points_df = pandas.read_csv(intermediate_df_path)
        points_df['year'] = year
        df_list.append(points_df)
    combined_df = pandas.concat(df_list)
    combined_df.to_csv(save_as, index=False)

    # clean up temporary files
    shutil.rmtree(temp_dir)


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
    try:
        raster_nodata = pygeoprocessing.get_raster_info(
            raster_path)['nodata'][0]
    except ValueError:
        print("Raster does not exist: {}".format(raster_path))
        return
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

    # set nodata values to NA
    report_table = pandas.DataFrame(data=sampled_precip_data_list)
    report_table = report_table.transpose()
    try:
        report_table.loc[
            numpy.isclose(report_table[raster_field_name], raster_nodata),
            raster_field_name] = None
    except TypeError:
        report_table[raster_field_name] = pandas.to_numeric(
            report_table[raster_field_name], errors='coerce')
        report_table.loc[
            numpy.isclose(report_table[raster_field_name], raster_nodata),
            raster_field_name] = None
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


def aggregate_covariates(
        point_shp_path, intermediate_dir_path, combined_covariate_table_path):
    """Calculate covariate values for each surface monitoring station.

    Parameters:
        point_shp_path (string): path to shapefile containing points that
            correspond to basins where values are aggregated
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

    if not os.path.exists(intermediate_dir_path):
        os.makedirs(intermediate_dir_path)

    avg_flow_zonal_mean_path = os.path.join(
        intermediate_dir_path, 'avg_flow_zonal_mean.csv')
    df_path_list.append(avg_flow_zonal_mean_path)
    if not os.path.exists(avg_flow_zonal_mean_path):
        zonal_mean_value_to_csv(
            _COVARIATE_PATH_DICT['average_flow'], 'average_flow_zonal_mean',
            avg_flow_zonal_mean_path)

    flash_flow_zonal_mean_path = os.path.join(
        intermediate_dir_path, 'flash_flow_zonal_mean.csv')
    df_path_list.append(flash_flow_zonal_mean_path)
    if not os.path.exists(flash_flow_zonal_mean_path):
        zonal_mean_value_to_csv(
            _COVARIATE_PATH_DICT['flash_flow'], 'flash_flow_zonal_mean',
            flash_flow_zonal_mean_path)

    climate_zone_df_path = os.path.join(
        intermediate_dir_path, 'climate_zone.csv')
    df_path_list.append(climate_zone_df_path)
    if not os.path.exists(climate_zone_df_path):
        raster_values_at_points(
            point_shp_path, _COVARIATE_PATH_DICT['climate_zones'],
            1, 'climate_zone', climate_zone_df_path)

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
        zonal_mean_value_to_csv(
            _COVARIATE_PATH_DICT['urban_extent_rc'],
            'proportion_urban', proportion_urban_df_path)

    precip_var_df_path = os.path.join(
        intermediate_dir_path, 'precip_variability.csv')
    df_path_list.append(precip_var_df_path)
    if not os.path.exists(precip_var_df_path):
        zonal_mean_value_to_csv(
            _COVARIATE_PATH_DICT['precip_variability'],
            'mean_precip_variability', precip_var_df_path)

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
        reclassify_countries_by_sanitation(
            _COVARIATE_PATH_DICT['countries_raster'],
            temp_val_dict['sanitation'])
        zonal_mean_value_to_csv(
            temp_val_dict['sanitation'], 'proportion_no_sanitation',
            proportion_no_sanitation_df_path)

    ndr_df_path = os.path.join(intermediate_dir_path, 'n_export.csv')
    df_path_list.append(ndr_df_path)
    if not os.path.exists(ndr_df_path):
        zonal_sum_to_csv(
            _COVARIATE_PATH_DICT['n_export'], 'basin_sum_n_export',
            ndr_df_path)

    merge_data_frame_list(df_path_list, combined_covariate_table_path)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def collect_covariates_from_rasters(
        point_shp_path, intermediate_dir_path, combined_covariate_table_path,
        groundwater=None):
    """Extract covariate values at station points from rasters.

    Parameters:
        point_shp_path (string): path to shapefile containing points features
            where intersecting pixel values should be collected
        intermediate_dir_path (string): path to folder where persistent
            intermediate outputs should be written, namely a csv for each
            covariate
        combined_covariate_table_path (string): path to location where the
            combined covariate table should be written, i.e. the table
            containing all covariate values
        groundwater (boolean): flag indicating whether the covariates are to be
            collected for variables relevant to groundwater

    Side effects:
        writes csv tables, one per covariate, inside `intermediate_dir_path`
        writes combined table of covariates at `combined_covariate_table_path`

    Returns:
        None

    """
    # list of paths to data frames each containing aggregated values for one
    # covariate
    df_path_list = []

    temp_dir = tempfile.mkdtemp()
    temp_val_dict = {
        'sanitation': os.path.join(temp_dir, 'sanitation.tif'),
    }

    if not os.path.exists(intermediate_dir_path):
        os.makedirs(intermediate_dir_path)

    n_export_path = os.path.join(intermediate_dir_path, 'n_export.csv')
    df_path_list.append(n_export_path)
    if not os.path.exists(n_export_path):
        raster_values_at_points(
            point_shp_path, _COVARIATE_PATH_DICT['n_export'], 1,
            'n_export', n_export_path)

    average_flow_path = os.path.join(intermediate_dir_path, 'average_flow.csv')
    df_path_list.append(average_flow_path)
    if not os.path.exists(average_flow_path):
        raster_values_at_points(
            point_shp_path, _COVARIATE_PATH_DICT['average_flow'], 1,
            'average_flow', average_flow_path)

    flash_flow_path = os.path.join(intermediate_dir_path, 'flash_flow.csv')
    df_path_list.append(flash_flow_path)
    if not os.path.exists(flash_flow_path):
        raster_values_at_points(
            point_shp_path, _COVARIATE_PATH_DICT['flash_flow'], 1,
            'flash_flow', flash_flow_path)

    precip_variability_path = os.path.join(
        intermediate_dir_path, 'precip_variability.csv')
    df_path_list.append(precip_variability_path)
    if not os.path.exists(precip_variability_path):
        raster_values_at_points(
            point_shp_path, _COVARIATE_PATH_DICT['precip_variability'], 1,
            'precip_variability', precip_variability_path)

    population_path = os.path.join(intermediate_dir_path, 'population.csv')
    df_path_list.append(population_path)
    if not os.path.exists(population_path):
        raster_values_at_points(
            point_shp_path, _COVARIATE_PATH_DICT['population'], 1,
            'population', population_path)

    urban_extent_path = os.path.join(intermediate_dir_path, 'urban_extent.csv')
    df_path_list.append(urban_extent_path)
    if not os.path.exists(urban_extent_path):
        raster_values_at_points(
            point_shp_path, _COVARIATE_PATH_DICT['urban_extent_rc'], 1,
            'proportion_urban', urban_extent_path)

    sanitation_path = os.path.join(intermediate_dir_path, 'sanitation.csv')
    df_path_list.append(sanitation_path)
    if not os.path.exists(sanitation_path):
        reclassify_countries_by_sanitation(
            _COVARIATE_PATH_DICT['countries_raster'],
            temp_val_dict['sanitation'])
        raster_values_at_points(
            point_shp_path, temp_val_dict['sanitation'], 1,
            'percent_no_sanitation', sanitation_path)

    cattle_path = os.path.join(intermediate_dir_path, 'cattle.csv')
    df_path_list.append(cattle_path)
    if not os.path.exists(cattle_path):
        raster_values_at_points(
            point_shp_path, _COVARIATE_PATH_DICT['cattle'],
            1, 'cattle', cattle_path)

    pigs_path = os.path.join(intermediate_dir_path, 'pigs.csv')
    df_path_list.append(pigs_path)
    if not os.path.exists(pigs_path):
        raster_values_at_points(
            point_shp_path, _COVARIATE_PATH_DICT['pigs'],
            1, 'pigs', pigs_path)


    if groundwater:
        depth_to_groundwater_path = os.path.join(
            intermediate_dir_path, 'depth_to_groundwater.csv')
        df_path_list.append(depth_to_groundwater_path)
        if not os.path.exists(depth_to_groundwater_path):
            raster_values_at_points(
                point_shp_path, _COVARIATE_PATH_DICT['depth_to_groundwater'],
                1, 'depth_to_groundwater', depth_to_groundwater_path)

        clay_percent_path = os.path.join(
            intermediate_dir_path, 'clay_percent.csv')
        df_path_list.append(clay_percent_path)
        if not os.path.exists(clay_percent_path):
            raster_values_at_points(
                point_shp_path, _COVARIATE_PATH_DICT['clay_percent'],
                1, 'clay_percent', clay_percent_path)

        sand_percent_path = os.path.join(
            intermediate_dir_path, 'sand_percent.csv')
        df_path_list.append(sand_percent_path)
        if not os.path.exists(sand_percent_path):
            raster_values_at_points(
                point_shp_path, _COVARIATE_PATH_DICT['sand_percent'],
                1, 'sand_percent', sand_percent_path)


    merge_data_frame_list(df_path_list, combined_covariate_table_path)

    # clean up temporary files
    shutil.rmtree(temp_dir)


def n_fert_at_points():
    """Generate tables summarizing n fertilizer trends at station points."""
    n_fert_surface_stn_csv = "F:/NCI_NDR/Data fertilizer Lu Tian/N_by_OBJECTID_WB_surface_stations_noxn_obs_objectid_adj.csv"
    # use surface station locations adjusted to intersect with Lu Tian data
    pt_shp_path = "F:/NCI_NDR/Data worldbank/station_data/WB_surface_stations_noxn_obs_objectid_shift_to_lu_tian.shp"
    summarize_n_fert_at_points(pt_shp_path, n_fert_surface_stn_csv)
    n_fert_groundwater_stn_csv = "F:/NCI_NDR/Data fertilizer Lu Tian/N_by_OBJECTID_WB_groundwater_stations_noxn_obs.csv"
    summarize_n_fert_at_points(
        _GROUNDWATER_STATION_SHP_PATH, n_fert_groundwater_stn_csv)


def collect_surface_covariates_points_snapped():
    """Aggregate covariates from the pixel containing the snapped station."""
    out_dir = 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/WB_station_snapped_5min_pixel'
    intermediate_dir_path = os.path.join(
        out_dir, 'intermediate_df_dir')
    combined_covariate_table_path = os.path.join(
        out_dir, 'combined_covariates.csv')
    collect_covariates_from_rasters(
        _SNAPPED_STATION_SHP_PATH, intermediate_dir_path,
        combined_covariate_table_path)


def collect_surface_covariates_points_orig(out_dir, resolution='5min'):
    """Aggregate covariates from the pixel containing the original station."""
    intermediate_dir_path = os.path.join(
        out_dir, 'intermediate_df_dir')
    combined_covariate_table_path = os.path.join(
        out_dir, 'combined_covariates.csv')
    global _COVARIATE_PATH_DICT
    if resolution == '30s':
        _COVARIATE_PATH_DICT = _COVARIATE_PATH_DICT_30s
    else:
        _COVARIATE_PATH_DICT = _COVARIATE_PATH_DICT_5min
    collect_covariates_from_rasters(
        _ORIG_STATION_SHP_PATH, intermediate_dir_path,
        combined_covariate_table_path)


def collect_surface_covariates_basin():
    """Aggregate covariates inside the basin delineated for each station."""
    out_dir = 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/Rafa_watersheds_v3'
    intermediate_dir_path = os.path.join(
        out_dir, 'intermediate_df_dir')
    combined_covariate_table_path = os.path.join(
        out_dir, 'combined_covariates.csv')
    aggregate_covariates(
        _SNAPPED_STATION_SHP_PATH, intermediate_dir_path,
        combined_covariate_table_path)


def collect_groundwater_covariates_orig(out_dir, resolution='5min'):
    """Aggregate covariates from piixel containing groundwater stations."""
    global _COVARIATE_PATH_DICT
    if resolution == '30s':
        _COVARIATE_PATH_DICT = _COVARIATE_PATH_DICT_30s
    else:
        _COVARIATE_PATH_DICT = _COVARIATE_PATH_DICT_5min

    _GROUNDWATER_STATION_SHP_PATH = "F:/NCI_NDR/Data groundwater merged/groundwater_GEMStat_Gu_USGS_Ouedraogo.shp"
    intermediate_dir_path = os.path.join(
        out_dir, 'intermediate_df_dir')
    combined_covariate_table_path = os.path.join(
        out_dir, 'combined_covariates.csv')
    collect_covariates_from_rasters(
        _GROUNDWATER_STATION_SHP_PATH, intermediate_dir_path,
        combined_covariate_table_path, groundwater=True)

def main():
    """Program entry point."""
    # surface_outdir = 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/WB_station_orig_5min_pixel'
    surface_outdir = 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/WB_station_orig_30s_pixel'
    # collect_surface_covariates_points_orig(surface_outdir, resolution='30s')
    # ground_outdir = 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/groundwater_GEMStat_Gu_USGS_Ouedraogo'
    ground_outdir = 'C:/Users/ginge/Dropbox/NatCap_backup/NCI WB/Aggregated_covariates/groundwater_GEMStat_Gu_USGS_Ouedraogo_30s_pixel'
    collect_groundwater_covariates_orig(ground_outdir, resolution='30s')


if __name__ == '__main__':
    main()
