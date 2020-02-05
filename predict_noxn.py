"""Use random forests to predict nitrate concentrations globally."""
import os
import tempfile
import shutil

from osgeo import gdal
import numpy
import pandas

from sklearn.ensemble import RandomForestRegressor
from joblib import dump, load
import pygeoprocessing

# full dict of covariate datasets at native resolution
_GLOBAL_COVARIATE_PATH_DICT = {
    'average_flow': "F:/NCI_NDR/Data streamflow FLO1K/average_flow_1990_2015.tif",
    'flash_flow': "F:/NCI_NDR/Data streamflow FLO1K/mean_div_range_1990_2015.tif",
    'precip_variability': "F:/NCI_NDR/Data precip Worldclim/wc2.0_bio_5m_15.tif",
    'population': "F:/NCI_NDR/Data population HYDE3.2/inhabitants_avg_1990_2015.tif",
    'proportion_urban': "F:/NCI_NDR/Data urban extent GRUMP/perc_urban_5min.tif",
    'percent_no_sanitation': "F:/NCI_NDR/Data sanitation/no_sanitation_provision_avg_2000-2015.tif",
    'depth_to_groundwater': "F:/NCI_NDR/Data HydroATLAS/gwt_cm_sav_level12.tif",
    'clay_percent': "F:/NCI_NDR/Data HydroATLAS/cly_pc_sav_level12.tif",
    'silt_percent': "F:/NCI_NDR/Data HydroATLAS/slt_pc_sav_level12.tif",
}

# n export rasters for a set of scenarios: each must be aggregated up to
# approximately 5 arc min resolution, each must have a unique basename
_N_EXPORT_PATH_DICT = {
    'baseline': "F:/NCI_NDR/Data NDR/nutrient_deficit_5min_cur_compressed_md5_031d4bb444325835315a2cc825be3fd4.tif",
    # continue with other scenarios here
}

# directory to hold temporary outputs
_PROCESSING_DIR = "F:/NCI_NDR/rf_processing"

# noxn and covariate observations for groundwater
_NOXN_PREDICTOR_GR_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_gr.csv"

# noxn and covariate observations for surface water
_NOXN_PREDICTOR_SURF_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_surf.csv"

# nodata value for inputs and result
_TARGET_NODATA = -1.0

def reclassify_nodata(target_path):
    """Reclassify the nodata value of a raster to _TARGET_NODATA.

    Check to see if the current nodata value for the raster indicated by
    `target_path` is _TARGET_NODATA. If not, convert all areas of nodata in the
    target raster to _TARGET_NODATA and set the nodata value for the taret
    raster to _TARGET_NODATA.

    Parameters:
        target_path (string): path to target raster

    Side effects:
        modifies the raster indicated by `target_path`

    Returns:
        None

    """
    def reclassify_op(target_raster):
        reclassified_raster = numpy.copy(target_raster)
        reclassify_mask = (target_raster == previous_nodata_value)
        reclassified_raster[reclassify_mask] = _TARGET_NODATA
        return reclassified_raster

    previous_nodata_value = pygeoprocessing.get_raster_info(
        target_path)['nodata'][0]
    if numpy.isclose(previous_nodata_value, _TARGET_NODATA):
        return

    file_handle, temp_path = tempfile.mkstemp(dir=_PROCESSING_DIR)
    shutil.copyfile(target_path, temp_path)
    pygeoprocessing.raster_calculator(
        [(temp_path, 1)], reclassify_op, target_path, gdal.GDT_Float32,
        _TARGET_NODATA)

    # clean up
    os.close(file_handle)
    os.remove(temp_path)


def prepare_covariates(predictor_names, aligned_covariate_dir):
    """Align covariate rasters and ensure that they all share one nodata value.

    Parameters:
        predictor_names (list): list of predictor names giving the order of
            covariates used to fit the random forests model.
        aligned_covariate_dir (string): path to directory where aligned
            covariate rasters should be stored

    Side effects:
        creates or modifies the directory `aligned_covariate_dir`
        creates or modifies rasters in `aligned_covariate_dir`, one for each
            predictor variable in `predictor_names`

    Returns:
        a dictionary whose keys are the covariate predictors in
            `predictor_names` and whose values contain paths to aligned
            covariate rasters

    """
    if not os.path.exists(aligned_covariate_dir):
        os.makedirs(aligned_covariate_dir)
    target_pixel_size = pygeoprocessing.get_raster_info(
        _GLOBAL_COVARIATE_PATH_DICT['average_flow'])['pixel_size']

    input_path_list = [
        _GLOBAL_COVARIATE_PATH_DICT[covar] for covar in predictor_names]
    aligned_path_list = [
        os.path.join(aligned_covariate_dir, os.path.basename(
            _GLOBAL_COVARIATE_PATH_DICT[covar])) for covar in predictor_names]
    if not all([os.path.isfile(r) for r in aligned_path_list]):
        pygeoprocessing.align_and_resize_raster_stack(
            input_path_list, aligned_path_list,
            ['near'] * len(input_path_list),
            target_pixel_size, 'intersection')
    # set nodata value for all aligned covariates to _TARGET_NODATA
    for path in aligned_path_list:
        reclassify_nodata(path)
    aligned_covariate_dict = {
        covar: os.path.join(aligned_covariate_dir, os.path.basename(
            _GLOBAL_COVARIATE_PATH_DICT[covar]))
        for covar in predictor_names
    }
    return aligned_covariate_dict


def predict_noxn(
        aligned_covariate_dir, n_export_path, rf_pickle_filename,
        predictor_names, output_path):
    """Generate global nitrate concentrations from a set of predictors.

    Parameters:
        aligned_covariate_dir (string): path to directory where aligned
            covariate rasters should be stored
        n_export_path (string): path to raster containing N export for a single
            scenario
        rf_pickle_filename (string): path to file on disk containing trained
            random forests model. This file should have the extension '.joblib'
        predictor_names (list): list of predictor names giving the order of
            covariates used to fit the random forests model. Any subsequent use
            of the model to make predictions on new covariate data must use
            covariates in this order.
        output_path (string): path to location on disk where output should be
            saved. This output will be a tif containing predicted nitrate
            concentration in mg/L

    Side effects:
        creates a tif giving predicted nitrate concentration in mg/L at the
            location given by `output_path`

    """
    def rf_predict_op(*covariate_list):
        """Use a trained random forest model to predict noxn."""
        invalid_mask = numpy.any(
            numpy.isclose(numpy.array(covariate_list), _TARGET_NODATA), axis=0)
        covar_arr = numpy.stack([r.ravel() for r in covariate_list], axis=1)
        noxn_arr = rf_model.predict(covar_arr)
        result = noxn_arr.reshape(covariate_list[0].shape)
        result[invalid_mask] = _TARGET_NODATA
        return result

    _GLOBAL_COVARIATE_PATH_DICT['n_export'] = n_export_path
    aligned_covariate_dict = prepare_covariates(
        predictor_names, aligned_covariate_dir)
    covariate_path_list = [
        aligned_covariate_dict[key] for key in predictor_names]

    rf_model = load(rf_pickle_filename)
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in covariate_path_list],
        rf_predict_op, output_path, gdal.GDT_Float32, _TARGET_NODATA)


def train_rf_model(
        noxn_predictor_df_path, max_features, rf_pickle_filename):
    """Train random forest model using observed data.

    Parameters:
        noxn_predictor_df_path (string): path to csv containing noxn
            obserations and predictor values for monitoring points. Must
            contain a column named 'noxn', containing observed nitrate
            concentrations, which is the target of the random forests model
        max_features (int): maximum number of variables to use in each tree
        rf_pickle_filename (string): location on disk where trained random
            forests model should be saved. This file should have the extension
            'joblib'

    Side effects:
        saves an instance of sklearn.ensemble.RandomForestRegressor, which has
            been trained on the observed data provided, and may be used to make
            predictions on new data, at the location `rf_pickle_filename`

    Returns:
        a list of predictor names giving the order of covariates used to fit
            the model. Any subsequent use of the model to make predictions on
            new covariate data must use covariates in this order.

    """
     # read data that was filtered and subsetted in R
    combined_df = pandas.read_csv(noxn_predictor_df_path)
    # drop rows containing missing data
    combined_df.dropna(inplace=True)

    # Separate data frame into response and predictors
    noxn_arr = numpy.array(combined_df['noxn'])
    predictor_df = combined_df.drop('noxn', axis=1)
    predictor_arr = numpy.array(predictor_df)

    # train random forests model
    # Instantiate model, using parameters chosen to match those used in R
    rf_model = RandomForestRegressor(
        random_state=42, n_estimators=500, criterion="mse",
        max_features=max_features, min_samples_leaf=5, oob_score=True,
        bootstrap=True)
    rf_model.fit(predictor_arr, noxn_arr)
    dump(rf_model, rf_pickle_filename)

    predictor_names = list(predictor_df.columns)
    return predictor_names


def surface_noxn_workflow():
    """Workflow to predict noxn in surface water."""
    surface_rf_path = os.path.join(_PROCESSING_DIR, 'surface_model.joblib')
    aligned_covariate_dir = os.path.join(
        _PROCESSING_DIR, 'aligned_covariates_surface')
    surface_max_features = 7
    surface_predictors = train_rf_model(
        _NOXN_PREDICTOR_SURF_DF_PATH, surface_max_features, surface_rf_path)
    for scenario_key in _N_EXPORT_PATH_DICT:
        n_export_path = _N_EXPORT_PATH_DICT[scenario_key]
        output_path = os.path.join(
            _PROCESSING_DIR, 'surface_noxn_{}.tif'.format(scenario_key))
        predict_noxn(
            aligned_covariate_dir, n_export_path, surface_rf_path,
            surface_predictors, output_path)


def groundwater_noxn_workflow():
    """Workflow to predict noxn in groundwater."""
    ground_rf_path = os.path.join(_PROCESSING_DIR, 'ground_model.joblib')
    aligned_covariate_dir = os.path.join(
        _PROCESSING_DIR, 'aligned_covariates_ground')
    ground_max_features = 2
    ground_predictors = train_rf_model(
        _NOXN_PREDICTOR_GR_DF_PATH, ground_max_features, ground_rf_path)
    for scenario_key in _N_EXPORT_PATH_DICT:
        n_export_path = _N_EXPORT_PATH_DICT[scenario_key]
        output_path = os.path.join(
            _PROCESSING_DIR, 'ground_noxn_{}.tif'.format(scenario_key))
        predict_noxn(
            aligned_covariate_dir, n_export_path, ground_rf_path,
            ground_predictors, output_path)


def main():
    """Program entry point."""
    surface_noxn_workflow()
    groundwater_noxn_workflow()


if __name__ == '__main__':
    __spec__ = None
    main()
