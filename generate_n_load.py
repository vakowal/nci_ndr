"""Generate n_load raster, an intermediate output of NDR."""
import os
import tempfile

import pandas
import numpy

from osgeo import gdal

import pygeoprocessing

# this value in n_load raster indicates substitute value from ag_load
_USE_AG_LOAD_ID = -999

# DATA_DIR = "F:/NCI_NDR/Data NDR/updated_5.18.20/n_load_inputs"
DATA_DIR = "C:/Users/ginge/Desktop/n_load_inputs"

PROCESSING_DIR = os.path.join(DATA_DIR, 'processing')

base_data_dict = {
    'lulc': os.path.join(
        DATA_DIR, 'LULC_map',
        'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_md5_1254d25f937e6d9bdee5779d377c5aa4.tif'),
    'biophysical_table': os.path.join(
        DATA_DIR, 'nci-NDR-biophysical_table_ESA_ARIES_RS3_md5_74d69f7e7dc829c52518f46a5a655fb8.csv'),
    'fertilizer_map': os.path.join(
        DATA_DIR, 'fertilizer_map',
        'scenarios050420_ExtensificationNapp_allcrops_rainfedfootprint_gapfilled_observedNappRevB_md5_1185e457751b672c67cc8c6bf7016d03.tif'),
    'natveg_load': os.path.join(
        DATA_DIR, 'load_n', 'load_n_fixedarea_currentpractices_global.tif'),
}


def generate_n_load():
    """
    The plan:

        * resample fertilizer map to match resolution of load_n
        * where load_n is -999, substitute the value from the fertilizer map
    """
    def ag_load_op(base_load_n_array, ag_load_array):
        """raster calculator replace _USE_AG_LOAD_ID with ag loads."""
        result = numpy.copy(base_load_n_array)
        if load_nodata is not None:
            nodata_load_mask = numpy.isclose(ag_load_array, load_nodata)
        else:
            nodata_load_mask = numpy.zeros(
                ag_load_array.shape, dtype=numpy.bool)
        ag_mask = (base_load_n_array == _USE_AG_LOAD_ID)
        result[ag_mask & ~nodata_load_mask] = (
            ag_load_array[ag_mask & ~nodata_load_mask])
        result[ag_mask & nodata_load_mask] = 0.0
        return result

    if not os.path.exists(PROCESSING_DIR):
        os.makedirs(PROCESSING_DIR)

    rescaled_fertilizer_path = os.path.join(
        DATA_DIR, 'fertilizer_map', 'align_to_load_n.tif')
    if not os.path.exists(rescaled_fertilizer_path):
        # align fertilizer map with load_n on natural veg
        with tempfile.NamedTemporaryFile(
                prefix='aligned_natveg',
                dir=PROCESSING_DIR) as natveg_load_temp_file:
            natveg_aligned_path = natveg_load_temp_file.name
            source_path_list = [
                base_data_dict['natveg_load'],
                base_data_dict['fertilizer_map']]
            aligned_path_list = [natveg_aligned_path, rescaled_fertilizer_path]
            target_pixel_size = pygeoprocessing.get_raster_info(
                base_data_dict['natveg_load'])['pixel_size']
            bounding_box = pygeoprocessing.get_raster_info(
                base_data_dict['natveg_load'])['bounding_box']
            pygeoprocessing.align_and_resize_raster_stack(
                source_path_list, aligned_path_list,
                ['near'] * len(source_path_list),
                target_pixel_size, bounding_box, raster_align_index=0)
    ag_load_raster_path = rescaled_fertilizer_path
    load_n_per_ha_raster_path = base_data_dict['natveg_load']

    load_nodata = pygeoprocessing.get_raster_info(
        ag_load_raster_path)['nodata'][0]

    nodata = pygeoprocessing.get_raster_info(
        load_n_per_ha_raster_path)['nodata'][0]

    target_ag_load_path = os.path.join(DATA_DIR, 'n_load_natveg_ag.tif')
    pygeoprocessing.raster_calculator(
        [(load_n_per_ha_raster_path, 1), (ag_load_raster_path, 1)],
        ag_load_op, target_ag_load_path,
        gdal.GDT_Float32, nodata)


if __name__ == "__main__":
    generate_n_load()
