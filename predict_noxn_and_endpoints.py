"""Predict nitrate concentrations and health/economic endpoints globally."""
import os
import tempfile
import shutil
import psutil

from osgeo import gdal
import numpy
import pandas

from sklearn.ensemble import RandomForestRegressor
import forestci
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
    'clay_percent': "F:/NCI_NDR/Data soil ISRIC/CLYPPT_M_sl1_10km_ll.tif",
    'sand_percent': "F:/NCI_NDR/Data soil ISRIC/SNDPPT_M_sl1_10km_ll.tif",
    'cattle': "F:/NCI_NDR/Data GLW/5_Ct_2010_Da.tif",
    'pigs': "F:/NCI_NDR/Data GLW/5_Pg_2010_Da.tif",
}

# base data dictionary
_BASE_DATA_PATH_DICT = {
    'water_source_table': "F:/NCI_NDR/Data endpoints/water_source_by_country_id.csv",
    'cancer_rate_table': "F:/NCI_NDR/Data endpoints/cancer_rate_by_country_id.csv",
    'countries_raster': "F:/NCI_NDR/Data national_boundaries/countries_iso3.tif",
    'population_raster': "F:/NCI_NDR/Data population Landscan/LandScan Global 2018/LandScan2018_WGS84_5min.tif",
}

# n export rasters for a set of scenarios: each must be aggregated up to
# approximately 5 arc min resolution, each must have a unique basename
_N_EXPORT_PATH_DICT = {
    'baseline': "F:/NCI_NDR/Data NDR/updated_3.2.20/sum_aggregate_to_0.084100_n_export_baseline_napp_rate_global_md5_b210146a5156422041eb7128c147512f.tif",
    'ag_expansion': "F:/NCI_NDR/Data NDR/updated_3.2.20/sum_aggregate_to_0.084100_n_export_ag_expansion_global_md5_ea15fb82df52d49a1d0c4ffe197cdd0d.tif",
    'ag_intensification': "F:/NCI_NDR/Data NDR/updated_3.2.20/sum_aggregate_to_0.084100_n_export_ag_intensification_global_md5_2734116e8c452f4c484ebcb574aab665.tif",
    'restoration': "F:/NCI_NDR/Data NDR/updated_3.2.20/sum_aggregate_to_0.084100_n_export_restoration_napp_rate_global_md5_7f9ddf313e414a68cbb8ba204101b190.tif",
}

# directory to hold temporary outputs
_PROCESSING_DIR = None

# noxn and covariate observations for groundwater
_NOXN_PREDICTOR_GR_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_gr.csv"

# noxn and covariate observations for surface water
_NOXN_PREDICTOR_SURF_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_surf_2000_2015.csv"

# nodata value for inputs and result
_TARGET_NODATA = -1.0


def subtract_op(raster1, raster2):
    """Subtract raster2 from raster1 element-wise."""
    valid_mask = (
        (raster1 != _TARGET_NODATA) &
        (raster2 != _TARGET_NODATA))
    result = numpy.empty(raster1.shape, dtype=numpy.float32)
    result[:] = _TARGET_NODATA
    result[valid_mask] = raster1[valid_mask] - raster2[valid_mask]
    return result


def add_op(raster1, raster2):
    """Add two rasters."""
    valid_mask = (
        (raster1 != _TARGET_NODATA) &
        (raster2 != _TARGET_NODATA))
    result = numpy.empty(raster1.shape, dtype=numpy.float32)
    result[:] = _TARGET_NODATA
    result[valid_mask] = raster1[valid_mask] + raster2[valid_mask]
    return result

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
        reclassify_mask = numpy.isclose(target_raster, previous_nodata_value)
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


def calc_drinking_water_source_raster(save_as):
    """Generate a raster giving percent of drinking water from surface water.

    This raster is generated by assigning country-level values in
    _BASE_DATA_PATH_DICT['water_source_table'] to countries identified in
    _BASE_DATA_PATH_DICT['countries_raster'].

    Parameters:
        save_as (string): location to save raster where country id values have
            been reclassified to proportion of drinking water supplied by
            surface water

    Side effects:
        creates or modifies a global raster at the location `save_as`

    Returns:
        None

    """
    drinking_df = pandas.read_csv(_BASE_DATA_PATH_DICT['water_source_table'])
    countryid_to_drinking = pandas.Series(
        drinking_df.perc_surface.values, index=drinking_df.id).to_dict()
    pygeoprocessing.reclassify_raster(
        (_BASE_DATA_PATH_DICT['countries_raster'], 1), countryid_to_drinking,
        save_as, gdal.GDT_Float32, _TARGET_NODATA)


def calc_background_cancer_rate_raster(save_as):
    """Generate a raster giving background cancer rate by country.

    This raster is generated by assigning country-level values in
    _BASE_DATA_PATH_DICT['cancer_rate_table'] to countries identified in
    _BASE_DATA_PATH_DICT['countries_raster'].

    Parameters:
        save_as (string): location to save raster where country id values have
            been reclassified to background cancer rate

    Side effects:
        creates or modifies a global raster at the location `save_as`

    Returns:
        None

    """
    rate_df = pandas.read_csv(_BASE_DATA_PATH_DICT['cancer_rate_table'])
    countryid_to_rate = pandas.Series(
        rate_df['Rate.per.person'].values, index=rate_df.id).to_dict()
    pygeoprocessing.reclassify_raster(
        (_BASE_DATA_PATH_DICT['countries_raster'], 1), countryid_to_rate,
        save_as, gdal.GDT_Float32, _TARGET_NODATA)


def calc_noxn_in_drinking_water(
        surface_noxn, ground_noxn, drinking_water_source):
    """Calculate noxn concentration in drinking water.

    Noxn in drinking water is calculated as the weighted average of noxn in
    surface and ground water, accounting for the fraction of drinking water
    obtained from surface vs groundwater sources.

    Parameters:
        surface_noxn (numpy.ndarray): noxn concentration in surface water, in
            mg/L
        ground_noxn (numpy.ndarray): noxn concentration in groundwater, in mg/L
        drinking_water_source (numpy.ndarray): fraction of drinking water
            obtained from surface water. The remainder is obtained from ground
            water

    Returns:
        noxn_in_drinking_water, a numpy array of noxn concentration in drinking
            water, in mg/L

    """
    noxn_in_drinking_water = numpy.empty(
        surface_noxn.shape, dtype=numpy.float32)
    noxn_in_drinking_water[:] = _TARGET_NODATA
    valid_mask = (
        (surface_noxn != _TARGET_NODATA) &
        (ground_noxn != _TARGET_NODATA) &
        (drinking_water_source != _TARGET_NODATA))
    noxn_in_drinking_water[valid_mask] = (
        surface_noxn[valid_mask] * drinking_water_source[valid_mask] +
        ground_noxn[valid_mask] * (1. - drinking_water_source[valid_mask]))
    return noxn_in_drinking_water


def calc_cancer_cases(noxn_in_drinking_water, population, background_rate):
    """Calculate noxn-attributable colorectal cancer cases per year.

    Convert nitrogen concentrations in drinking water into nitrate-
    attributable colorectal cancer cases per year, using a cancer slope factor
    of 0.04 cases per year per unit noxn from Temkin et al. (2019).

    Parameters:
        noxn_in_drinking_water (numpy.ndarray): estimated nitrate concentration
            in drinking water, in mg/L
        population (numpy.ndarray): inhabitants per grid cell
        background_rate (numpy.ndarray): country-specific background rate of
            colorectal cancer, in new cases per person per year

    Returns:
        cancer_cases, a numpy array of estimated nitrate-attributable cancer
            cases per year

    """
    cancer_cases = numpy.empty(population.shape, dtype=numpy.float32)
    cancer_cases[:] = _TARGET_NODATA
    valid_mask = (
        (noxn_in_drinking_water != _TARGET_NODATA) &
        (population != _TARGET_NODATA) &
        (background_rate != _TARGET_NODATA))
    cancer_cases[valid_mask] = (
        noxn_in_drinking_water[valid_mask] * population[valid_mask] *
        background_rate[valid_mask] * 0.04)
    return cancer_cases


def calc_treatment_costs(noxn_in_drinking_water, population):
    """Calculate nitrate abatement costs in USD per year.

    Convert nitrogen concentration in drinking water into nitrogen abatement
    costs. Where noxn is below the WHO drinking water limits for nitrates
    (11.3 mg/L noxn), it's assumed that no treatment is used and therefore the
    cost is zero. For areas exceeding the limits, the abatement cost is
    estimated by multiplying a constant abatement cost of 200 USD per person
    per year by the population.

    Parameters:
        noxn_in_drinking_water (numpy.ndarray): estimated nitrate concentration
            in drinking water, in mg/L
        population (numpy.ndarray): inhabitants per grid cell

    Returns:
        abatement_costs, a numpy array of estimated abatement costs for nitrate
            in drinking water in USD per year

    """
    abatement_costs = numpy.empty(population.shape, dtype=numpy.float32)
    abatement_costs[:] = _TARGET_NODATA
    valid_mask = (
        (noxn_in_drinking_water != _TARGET_NODATA) &
        (population != _TARGET_NODATA))
    threshold_conc = 11.3
    cost_per_person = 200
    abatement_costs[valid_mask] = 0
    exceeded_mask = (valid_mask & (noxn_in_drinking_water > threshold_conc))
    abatement_costs[exceeded_mask] = (
        population[exceeded_mask] * cost_per_person)
    return abatement_costs


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


def generate_confidence_intervals(
        noxn_predictor_df_path, noxn_path, aligned_covariate_dir,
        n_export_path, rf_pickle_filename, output_dir, basename):
    """Generate confidence intervals from random forest model.

    Parameters:
        noxn_path (string): path to raster containing predicted nitrate
            concentration for which the confidence intervals should be
            estimated
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
        output_dir (string): location on disk where outputs should be saved
        basename (string): basename for error, upper and lower bound rasters
            that should be created

    Side effects:
        creates a raster named 'noxn_error_<basename>.tif' in output_dir
        creates a raster named 'noxn_upper_bound_<basename>.tif' in output_dir
        creates a raster named 'noxn_lower_bound_<basename>.tif' in output_dir

    Returns:
        None

    """
    def calc_error_bound(*covariate_list):
        """Calculate +/- symmetric error bound from random forest model.

        The upper and lower error bound is calculated as the square root of the
        unbiased sampling variance calculated from the trained model. Upper and
        lower bounds can be calculated from this by adding this error to
        predicted noxn to get the upper bound; and by subtracting this error
        from predicted noxn to get the lower bound.

        Parameters:
            covariate_list (list of numpy.ndarrays): list of covariate
                predictors used to predict nitrate concentration

        Returns:
            error_bound, the square root of unbiased sampling variance

        """
        for covar_arr in covariate_list:
            numpy.place(covar_arr, numpy.isnan(covar_arr), [_TARGET_NODATA])
        invalid_mask = numpy.any(
            numpy.isclose(numpy.array(covariate_list), _TARGET_NODATA), axis=0)
        covar_arr = numpy.stack([r.ravel() for r in covariate_list], axis=1)
        mem_avail_mb = (psutil.virtual_memory().available >> 20) * 0.5
        unbiased_error = forestci.random_forest_error(
            rf_model, predictor_arr, covar_arr, memory_constrained=True,
            memory_limit=mem_avail_mb)
        error_arr = numpy.sqrt(unbiased_error)
        error_bound = error_arr.reshape(covariate_list[0].shape)
        error_bound[invalid_mask] = _TARGET_NODATA
        return error_bound

    # get predictor covariate array
    combined_df = pandas.read_csv(noxn_predictor_df_path)
    combined_df.dropna(inplace=True)
    predictor_df = combined_df.drop('noxn', axis=1)
    predictor_arr = numpy.array(predictor_df)
    predictor_names = list(predictor_df.columns)

    _GLOBAL_COVARIATE_PATH_DICT['n_export'] = n_export_path
    aligned_covariate_dict = prepare_covariates(
        predictor_names, aligned_covariate_dir)
    covariate_path_list = [
        aligned_covariate_dict[key] for key in predictor_names]

    rf_model = load(rf_pickle_filename)

    # with tempfile.NamedTemporaryFile(
    #         prefix='error_bound') as error_temp_file:
    #     error_path = error_temp_file.name
    error_path = os.path.join(
        output_dir, 'noxn_error_{}.tif'.format(basename))
    print("Calculating error bound for {} ...".format(basename))
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in covariate_path_list],
        calc_error_bound, error_path, gdal.GDT_Float32, _TARGET_NODATA,
        largest_block=0, raster_driver_creation_tuple=(('GTIFF', (
            'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=LZW',
            'BLOCKXSIZE=128', 'BLOCKYSIZE=128'))))

    print("Calculating lower bound ...")
    lower_bound_path = os.path.join(
        output_dir, 'noxn_lower_bound_{}.tif'.format(basename))
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [noxn_path, error_path]],
        subtract_op, lower_bound_path, gdal.GDT_Float32, _TARGET_NODATA)
    print("Calculating upper bound ...")
    upper_bound_path = os.path.join(
        output_dir, 'noxn_upper_bound_{}.tif'.format(basename))
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [noxn_path, error_path]],
        add_op, upper_bound_path, gdal.GDT_Float32, _TARGET_NODATA)

    # clean up
    # os.remove(error_path)


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
        # for a reason I don't understand, NaN values are cropping up here.
        # mask them out.
        for covar_arr in covariate_list:
            numpy.place(covar_arr, numpy.isnan(covar_arr), [_TARGET_NODATA])
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
    training_rsq = rf_model.score(predictor_arr, noxn_arr)
    print('R^2 Training Score: {}'.format(training_rsq))
    print('OOB Score: {}'.format(rf_model.oob_score_))
    dump(rf_model, rf_pickle_filename)

    predictor_names = list(predictor_df.columns)
    return predictor_names


def predict_surface_noxn(output_dir):
    """Workflow to predict noxn in surface water.

    Parameters:
        output_dir (string): path to directory where outputs should be stored

    Side effects:
        creates or modifies an output raster in `output_dir` for each scenario

    Returns:
        None

    """
    surface_rf_path = os.path.join(_PROCESSING_DIR, 'surface_model.joblib')
    aligned_covariate_dir = os.path.join(
        _PROCESSING_DIR, 'aligned_covariates_surface')
    surface_max_features = 2
    surface_predictors = train_rf_model(
        _NOXN_PREDICTOR_SURF_DF_PATH, surface_max_features, surface_rf_path)
    for scenario_key in _N_EXPORT_PATH_DICT:
        n_export_path = _N_EXPORT_PATH_DICT[scenario_key]
        noxn_path = os.path.join(
            output_dir, 'surface_noxn_{}.tif'.format(scenario_key))
        predict_noxn(
            aligned_covariate_dir, n_export_path, surface_rf_path,
            surface_predictors, noxn_path)
        generate_confidence_intervals(
            _NOXN_PREDICTOR_SURF_DF_PATH, noxn_path, aligned_covariate_dir,
            n_export_path, surface_rf_path, output_dir,
            'surface_{}'.format(scenario_key))


def predict_groundwater_noxn(output_dir):
    """Workflow to predict noxn in groundwater.

    Parameters:
        output_dir (string): path to directory where outputs should be stored

    Side effects:
        creates or modifies an output raster in `output_dir` for each scenario

    Returns:
        None

    """
    ground_rf_path = os.path.join(_PROCESSING_DIR, 'ground_model.joblib')
    aligned_covariate_dir = os.path.join(
        _PROCESSING_DIR, 'aligned_covariates_ground')
    ground_max_features = 2
    ground_predictors = train_rf_model(
        _NOXN_PREDICTOR_GR_DF_PATH, ground_max_features, ground_rf_path)
    for scenario_key in _N_EXPORT_PATH_DICT:
        n_export_path = _N_EXPORT_PATH_DICT[scenario_key]
        noxn_path = os.path.join(
            output_dir, 'ground_noxn_{}.tif'.format(scenario_key))
        predict_noxn(
            aligned_covariate_dir, n_export_path, ground_rf_path,
            ground_predictors, noxn_path)
        generate_confidence_intervals(
            _NOXN_PREDICTOR_GR_DF_PATH, noxn_path, aligned_covariate_dir,
            n_export_path, ground_rf_path, output_dir,
            'ground_{}'.format(scenario_key))


def calc_endpoints(output_dir):
    """Calculate health and economic endpoints from nitrate concentrations.

    Parameters:
        output_dir (string): path to directory where outputs should be stored

    Side effects:
        creates or modifies an output raster in `output_dir` for each scenario

    Returns:
        None

    """
    target_pixel_size = pygeoprocessing.get_raster_info(
        _BASE_DATA_PATH_DICT['countries_raster'])['pixel_size']

    drinking_water_source_path = os.path.join(
        output_dir, 'frac_surface.tif')
    if not os.path.isfile(drinking_water_source_path):
        calc_drinking_water_source_raster(drinking_water_source_path)

    background_cancer_rate_path = os.path.join(
        output_dir, 'bg_cancer_rate.tif')
    if not os.path.isfile(background_cancer_rate_path):
        calc_background_cancer_rate_raster(background_cancer_rate_path)

    for scenario_key in _N_EXPORT_PATH_DICT:
        surface_noxn_path = os.path.join(
            output_dir, 'surface_noxn_{}.tif'.format(scenario_key))
        ground_noxn_path = os.path.join(
            output_dir, 'ground_noxn_{}.tif'.format(scenario_key))

        # align all input rasters
        aligned_dir = tempfile.mkdtemp(prefix='aligned_', dir=_PROCESSING_DIR)
        input_path_dict = {
            'surface_noxn': surface_noxn_path,
            'ground_noxn': ground_noxn_path,
            'fraction_surface': drinking_water_source_path,
            'background_rate': background_cancer_rate_path,
            'population': _BASE_DATA_PATH_DICT['population_raster'],
        }
        aligned_input_dict = {
            key: os.path.join(aligned_dir, os.path.basename(path))
            for key, path in input_path_dict.items()
        }
        input_path_list = [
            input_path_dict[k] for k in sorted(input_path_dict.keys())]
        aligned_path_list = [
            aligned_input_dict[k] for k in sorted(input_path_dict.keys())]
        pygeoprocessing.align_and_resize_raster_stack(
            input_path_list, aligned_path_list,
            ['near'] * len(input_path_list), target_pixel_size, 'intersection')

        noxn_in_drinking_water_path = os.path.join(
            output_dir,
            'noxn_in_drinking_water_{}.tif'.format(scenario_key))
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                aligned_input_dict['surface_noxn'],
                aligned_input_dict['ground_noxn'],
                aligned_input_dict['fraction_surface']]],
            calc_noxn_in_drinking_water, noxn_in_drinking_water_path,
            gdal.GDT_Float32, _TARGET_NODATA)

        cancer_cases_path = os.path.join(
            output_dir, 'cancer_cases_{}.tif'.format(scenario_key))
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                noxn_in_drinking_water_path,
                aligned_input_dict['population'],
                aligned_input_dict['background_rate']]],
            calc_cancer_cases, cancer_cases_path, gdal.GDT_Float32,
            _TARGET_NODATA)

        abatement_cost_path = os.path.join(
            output_dir, 'abatement_costs_{}.tif'.format(scenario_key))
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in [
                noxn_in_drinking_water_path,
                aligned_input_dict['population']]],
            calc_treatment_costs, abatement_cost_path, gdal.GDT_Float32,
            _TARGET_NODATA)

        # clean up
        shutil.rmtree(aligned_dir)


def main():
    """Program entry point."""
    outer_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015"
    output_dir = os.path.join(outer_dir, 'output')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    global _PROCESSING_DIR
    _PROCESSING_DIR = os.path.join(outer_dir, 'intermediate')
    if not os.path.exists(_PROCESSING_DIR):
        os.makedirs(_PROCESSING_DIR)
    predict_surface_noxn(output_dir)
    predict_groundwater_noxn(output_dir)
    calc_endpoints(output_dir)


def test_confidence_intervals():
    """Throwaway code to test calculation of confidence intervals."""
    outer_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015"
    output_dir = os.path.join(outer_dir, 'output')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    global _PROCESSING_DIR
    _PROCESSING_DIR = os.path.join(outer_dir, 'intermediate')
    if not os.path.exists(_PROCESSING_DIR):
        os.makedirs(_PROCESSING_DIR)
    # for testing, work with baseline only
    global _N_EXPORT_PATH_DICT
    _N_EXPORT_PATH_DICT.pop('ag_expansion', None)
    _N_EXPORT_PATH_DICT.pop('ag_intensification', None)
    _N_EXPORT_PATH_DICT.pop('restoration', None)
    predict_surface_noxn(output_dir)


if __name__ == '__main__':
    __spec__ = None  # for running with pdb
    # main()
    test_confidence_intervals()
