"""Miscellaneous data processing tasks for NCI-NDR project."""
import os
import numpy
import pandas

import pygeoprocessing


def raster_list_mean(raster_band_tuple_list, target_path):
    """Calculate the pixel-level mean across rasters in a list.

    Calculate the mean value per pixel across rasters in `raster_list`.  Areas
    containing a nodata value in any of the rasters in `raster_list` will be
    nodata in the result.

    Parameters:
        raster_band_tuple_list (list): list of (str, ing) tuples indicating
            path and band of rasters to summarize
        target_path (float or int):path to location to store the result

    Side effects:
        modifies or creates the raster indicated by `target_path`

    Returns:
        None

    """
    def raster_mean_op(*raster_list):
        invalid_mask = numpy.any(
            numpy.isclose(numpy.array(raster_list), input_nodata), axis=0)
        raster_mean = numpy.mean(raster_list, axis=0)
        raster_mean[invalid_mask] = target_nodata
        return raster_mean

    input_nodata = pygeoprocessing.geoprocessing.get_raster_info(
        raster_band_tuple_list[0][0])['nodata'][0]
    datatype_target = pygeoprocessing.geoprocessing.get_raster_info(
        raster_band_tuple_list[0][0])['datatype']

    pygeoprocessing.raster_calculator(
        raster_band_tuple_list, raster_mean_op, target_path, datatype_target,
        input_nodata)


def reduce_FLO1K():
    """Reduce the FLO1K dataset by calculating summary statistics across years.

    The FLO1K dataset of Barbarossa et al. contains annual average, minimum,
    and maximum flow values over the period 1960-2015. Reduce the dataset by
    calculating summary statistics over the time period 1990-2015.

    """
    out_dir = "F:/NCI_NDR/Data streamflow FLO1K"
    # time period over which to calculate summary statistics
    min_year = 1990
    av_flow_path = "F:/NCI_NDR/Data streamflow FLO1K/FLO1K.5min.ts.1960.2015.qav.nc"

    # bands in the dataset are organized by increasing time:
    #    band 1 = 1960, band 56 = 2015
    num_bands = pygeoprocessing.geoprocessing.get_raster_info(
        av_flow)['n_bands']

    # band corresponding to first year in temporal subset
    min_band = min_year - 1959
    temporal_subset = [b for b in range(min_band, num_bands + 1)]

    # average flow across temporal subset
    raster_band_tuple_list = [(av_flow_path, band) for band in temporal_subset]
    target_path = os.path.join(
        out_dir, "average_flow_{}_2015.tif".format(min_year))
    raster_list_mean(raster_band_tuple_list, target_path)


def compare_basin_size():
    """Compare the size of catchment areas draining to monitoring points.

    Calculate the size of drainage areas via zonal statistics on a raster
    containing pixel size in square km.  Compare this to reported basin sizes
    in the World Bank dataset.

    """
    basin_shp_path = "F:/NCI_NDR/Watersheds_DRT/WB_surface_stations_noxn_obs_snapped_joined_to_GLORIC_unique_fromML_ContribArea.shp"
    pixel_area_raster_path = "F:/NCI_NDR/Data flow direction DRT/pixel_area_km2.tif"

    zonal_stats_dict = pygeoprocessing.geoprocessing.zonal_statistics(
        (pixel_area_raster_path, 1), basin_shp_path)
    zonal_sum_table = {
        feature: inner_dict['sum'] for feature, inner_dict in
        zonal_stats_dict.items()}
    zonal_sum_df = pandas.DataFrame.from_dict(zonal_sum_table, orient='index')

    save_as = "F:/NCI_NDR/Watersheds_DRT/basin_area_km2_zonal_stats.csv"
    zonal_sum_df.to_csv(save_as)


if __name__ == '__main__':
    compare_basin_size()
