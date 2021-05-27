"""Predict nitrate concentrations and health/economic endpoints globally."""
import os
import tempfile
import shutil

from osgeo import gdal
import numpy
import pandas

# from sklearn.ensemble import RandomForestRegressor
# from joblib import dump, load
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
# _N_EXPORT_PATH_DICT = {
#     'baseline': "F:/NCI_NDR/Data NDR/updated_3.2.20/sum_aggregate_to_0.084100_n_export_baseline_napp_rate_global_md5_b210146a5156422041eb7128c147512f.tif",
#     'ag_expansion': "F:/NCI_NDR/Data NDR/updated_3.2.20/sum_aggregate_to_0.084100_n_export_ag_expansion_global_md5_ea15fb82df52d49a1d0c4ffe197cdd0d.tif",
#     'ag_intensification': "F:/NCI_NDR/Data NDR/updated_3.2.20/sum_aggregate_to_0.084100_n_export_ag_intensification_global_md5_2734116e8c452f4c484ebcb574aab665.tif",
#     'restoration': "F:/NCI_NDR/Data NDR/updated_3.2.20/sum_aggregate_to_0.084100_n_export_restoration_napp_rate_global_md5_7f9ddf313e414a68cbb8ba204101b190.tif",
# }
# file prefix identifying the baseline scenario
_N_EXPORT_BASELINE_KEY = 'baseline_currentpractices'
_N_EXPORT_PATH_LIST = [
    'extensification_bmps_irrigated',
    'extensification_bmps_rainfed',
    'extensification_current_practices',
    'extensification_intensified_irrigated',
    'extensification_intensified_rainfed',
    'fixedarea_bmps_irrigated',
    'fixedarea_bmps_rainfed',
    'fixedarea_intensified_irrigated',
    'fixedarea_intensified_rainfed',
    'grazing_expansion',
    'restoration',
    'sustainable_currentpractices'
]

# directory to hold temporary outputs
_PROCESSING_DIR = None

# noxn and covariate observations for groundwater
_NOXN_PREDICTOR_GR_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_gr.csv"

# noxn and covariate observations for surface water
_NOXN_PREDICTOR_SURF_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_surf_2000_2015.csv"

# nodata value for inputs and result
_TARGET_NODATA = -1.0

# nodata value for noxn rasters calculated outside this script
_NOXN_NODATA = None


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
        (~numpy.isclose(surface_noxn, _NOXN_NODATA)) &
        (~numpy.isclose(ground_noxn, _NOXN_NODATA)) &
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


def prepare_covariates_March2021_NDR():
    """Align new NDR outputs with existing aligned covariates."""
    nexport_dir = "F:/NCI_NDR/Data NDR/updated_3.27.21/resampled_by_Becky"
    raw_nexport_path_list = [
        os.path.join(nexport_dir, f) for f in os.listdir(nexport_dir) if
        f.endswith('.tif')]

    aligned_covariate_dir = "C:/aligned_NDR"  # "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/"
    if not os.path.exists(aligned_covariate_dir):
        os.makedirs(aligned_covariate_dir)
    template_raster = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/template.tif"
    target_pixel_size = pygeoprocessing.get_raster_info(
        template_raster)['pixel_size']
    target_bb = pygeoprocessing.get_raster_info(
        template_raster)['bounding_box']

    input_path_list = ([template_raster] + raw_nexport_path_list)
    aligned_path_list = [
        os.path.join(aligned_covariate_dir, os.path.basename(
            f)) for f in input_path_list]
    pygeoprocessing.align_and_resize_raster_stack(
        input_path_list, aligned_path_list, ['near'] * len(input_path_list),
        target_pixel_size, target_bb, raster_align_index=0)

    # set nodata value for all aligned covariates to _TARGET_NODATA
    for path in aligned_path_list:
        reclassify_nodata(path)

    # copy aligned rasters to "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/"


def prepare_covariates_May2020_NDR():
    """Align covariate rasters and ensure that they all share one nodata value.

    Side effects:
        Creates aligned versions of covariate and N export rasters.

    Returns:
        None

    """
    nexport_dir = "F:/NCI_NDR/Data NDR/updated_5.18.20"
    raw_nexport_path_list = [
        os.path.join(nexport_dir, f) for f in os.listdir(nexport_dir) if
        f.endswith('.tif')]

    aligned_covariate_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/aligned_covariates_ground/"
    if not os.path.exists(aligned_covariate_dir):
        os.makedirs(aligned_covariate_dir)
    target_pixel_size = pygeoprocessing.get_raster_info(
        _GLOBAL_COVARIATE_PATH_DICT['average_flow'])['pixel_size']

    input_path_list = (
        raw_nexport_path_list +
        [_GLOBAL_COVARIATE_PATH_DICT[covar] for covar in
        _GLOBAL_COVARIATE_PATH_DICT])
    aligned_path_list = [
        os.path.join(aligned_covariate_dir, os.path.basename(
            f)) for f in input_path_list]
    pygeoprocessing.align_and_resize_raster_stack(
        input_path_list, aligned_path_list, ['near'] * len(input_path_list),
        target_pixel_size, 'intersection')
    # set nodata value for all aligned covariates to _TARGET_NODATA
    for path in aligned_path_list:
        reclassify_nodata(path)


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


def generate_confidence_intervals(noxn_path, se_path, output_dir, basename):
    """Generate confidence intervals from predicted noxn and standard error.

    Parameters:
        noxn_path (string): path to raster containing predicted nitrate
            concentration for which the confidence intervals should be
            estimated
        se_path (string): path to raster containing standard error estimate for
            the predicted nitrate concentrations
        output_dir (string): location on disk where outputs should be saved
        basename (string): basename for upper and lower bound rasters
            that should be created

    Side effects:
        creates a raster named 'noxn_95%_upper_bound_<basename>.tif' in
            output_dir
        creates a raster named 'noxn_95%_lower_bound_<basename>.tif' in
            output_dir

    Returns:
        None

    """
    def lower_ci_op(mean, standard_error):
        """Calculate lower bound of 95% confidence interval from mean and se."""
        valid_mask = (
            (~numpy.isclose(mean, noxn_nodata)) &
            (~numpy.isclose(standard_error, se_nodata)))
        result = numpy.empty(mean.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        result[valid_mask] = (
            mean[valid_mask] - 1.96 * standard_error[valid_mask])
        return result

    def upper_ci_op(mean, standard_error):
        """Calculate upper bound of 95% confidence interval from mean and se."""
        valid_mask = (
            (~numpy.isclose(mean, noxn_nodata)) &
            (~numpy.isclose(standard_error, se_nodata)))
        result = numpy.empty(mean.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        result[valid_mask] = (
            mean[valid_mask] + 1.96 * standard_error[valid_mask])
        return result

    noxn_nodata = pygeoprocessing.get_raster_info(noxn_path)['nodata'][0]
    se_nodata = pygeoprocessing.get_raster_info(se_path)['nodata'][0]
    lower_bound_path = os.path.join(
        output_dir, 'noxn_95%_lower_bound_{}.tif'.format(basename))
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [noxn_path, se_path]],
        lower_ci_op, lower_bound_path, gdal.GDT_Float32, _TARGET_NODATA)
    upper_bound_path = os.path.join(
        output_dir, 'noxn_95%_upper_bound_{}.tif'.format(basename))
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [noxn_path, se_path]],
        upper_ci_op, upper_bound_path, gdal.GDT_Float32, _TARGET_NODATA)


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


def calc_endpoints(noxn_dir, output_dir):
    """Calculate health and economic endpoints from nitrate concentrations.

    Parameters:
        noxn_dir (string): path to directory where rasters containing noxn in
            surface and groundwater are stored. These rasters may have been
            generated by this script or by another tool.
        output_dir (string): path to directory where outputs should be stored

    Side effects:
        creates or modifies an output raster in `output_dir` for each scenario

    Returns:
        None

    """
    target_pixel_size = pygeoprocessing.get_raster_info(
        _BASE_DATA_PATH_DICT['countries_raster'])['pixel_size']

    drinking_water_source_path = os.path.join(
        _PROCESSING_DIR, 'frac_surface.tif')
    if not os.path.isfile(drinking_water_source_path):
        calc_drinking_water_source_raster(drinking_water_source_path)

    background_cancer_rate_path = os.path.join(
        _PROCESSING_DIR, 'bg_cancer_rate.tif')
    if not os.path.isfile(background_cancer_rate_path):
        calc_background_cancer_rate_raster(background_cancer_rate_path)

    for scenario_key in ([_N_EXPORT_BASELINE_KEY] + _N_EXPORT_PATH_LIST):
        surface_noxn_path = os.path.join(
            noxn_dir, 'surface_noxn_{}.tif'.format(scenario_key))
        ground_noxn_path = os.path.join(
            noxn_dir, 'ground_noxn_{}.tif'.format(scenario_key))

        # set _NOXN_NODATA from noxn rasters
        surface_noxn_nodata = pygeoprocessing.get_raster_info(
            surface_noxn_path)['nodata'][0]
        ground_noxn_nodata = pygeoprocessing.get_raster_info(
            ground_noxn_path)['nodata'][0]
        assert numpy.isclose(surface_noxn_nodata, ground_noxn_nodata), (
            "nodata values differ between surface and ground noxn")
        global _NOXN_NODATA
        _NOXN_NODATA = surface_noxn_nodata

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


def calc_dir_change_masked_scenario_noxn(
        scenario_noxn_path, baseline_noxn_path,
        scenario_n_export_path, baseline_n_export_path,
        mask_path, masked_noxn_path):
    """Generate noxn raster masked by direction of change of n export.

    Substitute baseline values into the scenario raster in areas
    where modeled N export goes up, but nitrate concentration goes down; and in
    areas where modeled N export goes down, but nitrate concentration goes up.

    Parameters:
        scenario_noxn_path (string): path to raster containing predicted
            nitrate according to the scenario
        baseline_noxn_path (string) path to raster containing predicted
            nitrate according to baseline
        scenario_n_export_path (string): path to raster containing modeled n
            export for the scenario
        baseline_n_export_path (string): path to raster containing modeled n
            export for baseline
        mask_path (string): location where binary mask should be saved. In this
            raster, pixels with value=1 indicate that the baseline value should
            be substituted
        masked_noxn_path (string): path to location to save raster containing
            noxn for the scenario masked according to direction of change in
            n export.

    Side effects:
        creates a raster at `mask_path`
        creates a raster at `masked_noxn_path`

    Returns:
        None

    """
    def calc_mask_op(
            scen_noxn, baseline_noxn, scen_nexport, baseline_nexport):
        """Generate the mask where direction of change in noxn is `wrong`."""
        valid_mask = (
            (~numpy.isclose(scen_noxn, noxn_nodata)) &
            (~numpy.isclose(baseline_noxn, noxn_nodata)) &
            (~numpy.isclose(scen_nexport, nexport_s_nodata)) &
            (~numpy.isclose(baseline_nexport, nexport_b_nodata)))
        noxn_diff = scen_noxn - baseline_noxn
        n_export_diff = scen_nexport - baseline_nexport
        # N export goes up, but nitrate concentration goes down
        dir1_mask = (
            (n_export_diff > 0) &
            (noxn_diff < 0) &
            valid_mask)
        # N export goes down, but nitrate concentration goes up
        dir2_mask = (
            (n_export_diff < 0) &
            (noxn_diff > 0) &
            valid_mask)
        mask_ar = numpy.zeros(scen_noxn.shape, dtype=numpy.byte)
        mask_ar[dir1_mask] = 1
        mask_ar[dir2_mask] = 1
        return mask_ar

    def apply_mask_op(
            scen_noxn, baseline_noxn, mask):
        """Apply the mask. Substitute baseline into scenario where mask==1."""
        filter_mask = (
            (~numpy.isclose(scen_noxn, noxn_nodata)) &
            (~numpy.isclose(baseline_noxn, noxn_nodata)) &
            (mask == 1))
        result = scen_noxn.copy()
        result[filter_mask] = baseline_noxn[filter_mask]
        return result

    noxn_nodata = pygeoprocessing.get_raster_info(
        scenario_noxn_path)['nodata'][0]
    assert pygeoprocessing.get_raster_info(
        baseline_noxn_path)['nodata'][0] == noxn_nodata, (
        "Input noxn rasters must share nodata value")
    nexport_s_nodata = pygeoprocessing.get_raster_info(
        scenario_n_export_path)['nodata'][0]
    nexport_b_nodata = pygeoprocessing.get_raster_info(
        baseline_n_export_path)['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            scenario_noxn_path, baseline_noxn_path, scenario_n_export_path,
            baseline_n_export_path]],
        calc_mask_op, mask_path, gdal.GDT_Byte, -1)
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            scenario_noxn_path, baseline_noxn_path, mask_path]],
        apply_mask_op, masked_noxn_path, gdal.GDT_Float32, noxn_nodata)

def calc_masked_scenario_noxn(
        scenario_noxn_path, scenario_se_path, baseline_noxn_path,
        baseline_se_path, sig_diff_scenario_path, sig_diff_filled_path):
    """Generate raster of noxn masked by significant difference from baseline.

    Use standard error of prediction for noxn in scenario and baseline to
    identify areas where the scenario is significantly different from baseline
    (i.e., where confidence intervals of the scenario and baseline do not
    overlap). Create two masked versions of the predicted noxn for the
    scenario: one where areas in the scenario raster where confidence intervals
    overlap are set to nodata, and one where areas in the scenario raster where
    confidence intervals overlap are filled with baseline noxn values.

    Parameters:
        scenario_noxn_path (string): path to raster containing predicted
            nitrate according to the scenario
        scenario_se_path (string): path to raster containing standard error
            estimate for the scenario
        baseline_noxn_path (string) path to raster containing predicted
            nitrate according to baseline
        baseline_se_path (string): path to raster containing standard error
            estimate for baseline
        sig_diff_scenario_path (string): path to location to save raster
            containing noxn for the scenario only in areas where the scenario
            is significantly different from baseline
        sig_diff_filled_path (string): path to location to save raster
            containing noxn for the scenario only in areas where the scenario
            is significantly different from baseline, other areas filled with
            baseline values

    Side effects:
        creates a raster at `sig_diff_scenario_path`
        creates a raster at `sig_diff_filled_path`

    Returns:
        None

    """
    def sig_diff_masked_op(
            scenario_mean, scenario_se, baseline_mean, baseline_se):
        """Restrict scenario mean to areas of sig difference from baseline."""
        valid_mask = (
            (~numpy.isclose(scenario_mean, noxn_nodata)) &
            (~numpy.isclose(scenario_se, se_nodata)) &
            (~numpy.isclose(baseline_mean, noxn_nodata)) &
            (~numpy.isclose(baseline_se, se_nodata)))

        scenario_lower_ci = numpy.empty(
            scenario_mean.shape, dtype=numpy.float32)
        scenario_lower_ci[:] = _TARGET_NODATA
        scenario_lower_ci[valid_mask] = (
            scenario_mean[valid_mask] - 1.96 * scenario_se[valid_mask])
        scenario_upper_ci = numpy.empty(
            scenario_mean.shape, dtype=numpy.float32)
        scenario_upper_ci[:] = _TARGET_NODATA
        scenario_upper_ci[valid_mask] = (
            scenario_mean[valid_mask] + 1.96 * scenario_se[valid_mask])

        baseline_lower_ci = numpy.empty(
            baseline_mean.shape, dtype=numpy.float32)
        baseline_lower_ci[:] = _TARGET_NODATA
        baseline_lower_ci[valid_mask] = (
            baseline_mean[valid_mask] - 1.96 * baseline_se[valid_mask])
        baseline_upper_ci = numpy.empty(
            baseline_mean.shape, dtype=numpy.float32)
        baseline_upper_ci[:] = _TARGET_NODATA
        baseline_upper_ci[valid_mask] = (
            baseline_mean[valid_mask] + 1.96 * baseline_se[valid_mask])

        result = numpy.empty(scenario_mean.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        sig_smaller_mask = (
            (scenario_upper_ci < baseline_lower_ci) &
            valid_mask)
        sig_larger_mask = (
            (scenario_lower_ci > baseline_upper_ci) &
            valid_mask)
        result[sig_smaller_mask] = scenario_mean[sig_smaller_mask]
        result[sig_larger_mask] = scenario_mean[sig_larger_mask]
        return result

    def sig_diff_filled_op(
            scenario_mean, scenario_se, baseline_mean, baseline_se):
        """Fill nonsignificantly different areas with baseline mean."""
        valid_mask = (
            (~numpy.isclose(scenario_mean, noxn_nodata)) &
            (~numpy.isclose(scenario_se, se_nodata)) &
            (~numpy.isclose(baseline_mean, noxn_nodata)) &
            (~numpy.isclose(baseline_se, se_nodata)))

        scenario_lower_ci = numpy.empty(
            scenario_mean.shape, dtype=numpy.float32)
        scenario_lower_ci[:] = _TARGET_NODATA
        scenario_lower_ci[valid_mask] = (
            scenario_mean[valid_mask] - 1.96 * scenario_se[valid_mask])
        scenario_upper_ci = numpy.empty(
            scenario_mean.shape, dtype=numpy.float32)
        scenario_upper_ci[:] = _TARGET_NODATA
        scenario_upper_ci[valid_mask] = (
            scenario_mean[valid_mask] + 1.96 * scenario_se[valid_mask])

        baseline_lower_ci = numpy.empty(
            baseline_mean.shape, dtype=numpy.float32)
        baseline_lower_ci[:] = _TARGET_NODATA
        baseline_lower_ci[valid_mask] = (
            baseline_mean[valid_mask] - 1.96 * baseline_se[valid_mask])
        baseline_upper_ci = numpy.empty(
            baseline_mean.shape, dtype=numpy.float32)
        baseline_upper_ci[:] = _TARGET_NODATA
        baseline_upper_ci[valid_mask] = (
            baseline_mean[valid_mask] + 1.96 * baseline_se[valid_mask])

        result = numpy.empty(scenario_mean.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        result[valid_mask] = baseline_mean[valid_mask]  # fill with baseline
        sig_smaller_mask = (
            (scenario_upper_ci < baseline_lower_ci) &
            valid_mask)
        sig_larger_mask = (
            (scenario_lower_ci > baseline_upper_ci) &
            valid_mask)
        result[sig_smaller_mask] = scenario_mean[sig_smaller_mask]
        result[sig_larger_mask] = scenario_mean[sig_larger_mask]
        return result

    noxn_nodata = pygeoprocessing.get_raster_info(
        scenario_noxn_path)['nodata'][0]
    assert pygeoprocessing.get_raster_info(
        baseline_noxn_path)['nodata'][0] == noxn_nodata, (
        "Input noxn rasters must share nodata value")
    se_nodata = pygeoprocessing.get_raster_info(
        scenario_se_path)['nodata'][0]
    assert pygeoprocessing.get_raster_info(
        baseline_se_path)['nodata'][0] == se_nodata, (
        "Input se rasters must share nodata value")

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            scenario_noxn_path, scenario_se_path,
            baseline_noxn_path, baseline_se_path]],
        sig_diff_masked_op, sig_diff_scenario_path, gdal.GDT_Float32,
        _TARGET_NODATA)
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in [
            scenario_noxn_path, scenario_se_path,
            baseline_noxn_path, baseline_se_path]],
        sig_diff_filled_op, sig_diff_filled_path, gdal.GDT_Float32,
        _TARGET_NODATA)


def predict_noxn_and_endpoints():
    """Predict noxn in surface and groundwater, and calculate endpoints.

    Train random forest models for surface and groundwater from observed
    noxn concentrations. Use the trained models to predict noxn concentrations
    globally. Translate predicted noxn concentrations to health and financial
    endpoints.

    Returns:
        None

    """
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
    calc_endpoints(output_dir, output_dir)


def predict_endpoints_only():
    """Calculate endpoints from noxn in surface and groundwater.

    Translate predicted noxn concentrations to health and financial endpoints.

    Returns:
        None

    """
    noxn_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/R_ranger_pred"
    output_dir = os.path.join(noxn_dir, 'endpoints')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    global _PROCESSING_DIR
    _PROCESSING_DIR = os.path.join(noxn_dir, 'intermediate')
    if not os.path.exists(_PROCESSING_DIR):
        os.makedirs(_PROCESSING_DIR)
    calc_endpoints(noxn_dir, output_dir)


def confidence_interval_wrapper():
    """Generate lower and upper 95% confidence intervals for predicted noxn."""
    noxn_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.2.20/subset_2000_2015/R_ranger_pred"
    output_dir = noxn_dir
    for scenario_key in _N_EXPORT_PATH_DICT:
        surface_noxn_path = os.path.join(
            noxn_dir, 'surface_noxn_{}.tif'.format(scenario_key))
        surface_noxn_se_path = os.path.join(
            noxn_dir, 'surface_noxn_se_{}.tif'.format(scenario_key))
        basename = 'surface_{}'.format(scenario_key)
        generate_confidence_intervals(
            surface_noxn_path, surface_noxn_se_path, output_dir, basename)


        ground_noxn_path = os.path.join(
            noxn_dir, 'ground_noxn_{}.tif'.format(scenario_key))
        ground_noxn_se_path = os.path.join(
            noxn_dir, 'ground_noxn_se_{}.tif'.format(scenario_key))
        basename = 'ground_{}'.format(scenario_key)
        generate_confidence_intervals(
            ground_noxn_path, ground_noxn_se_path, output_dir, basename)


def dir_change_mask_endpoints_workflow():
    """Calc endpoints from noxn masked according to direction of change.

    New filtering criteria. Apply the mask to surface and ground water
    separately. Substitute baseline values into the scenario raster in areas
    where modeled N export goes up, but nitrate concentration goes down; and in
    areas where modeled N export goes down, but nitrate concentration goes up.

    Side effects:
        creates the following rasters for each scenario:
            the binary mask that was applied according to direction of change,
                for surface water
            the binary mask that was applied according to direction of change,
                for ground water
            masked surface water
            masked ground water
            cancer cases

    Returns:
        None

    """
    # directory containing raw nitrate concentrations
    predicted_noxn_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/R_ranger_pred"
    # N export for each scenario
    nexport_dir = "F:/NCI_NDR/Data NDR/updated_3.27.21/resampled_by_Becky/renamed"
    n_export_pattern = "compressed_{}.tif"
    # directory for aligned rasters
    aligned_dir = "C:/NCI_NDR/aligned_export_noxn"
    if not os.path.exists(aligned_dir):
        os.makedirs(aligned_dir)
    # directory for masks
    mask_dir = "C:/NCI_NDR/filter_mask"
    if not os.path.exists(mask_dir):
        os.makedirs(mask_dir)
    # directory containing masked, filled surface and ground concentration
    masked_filled_dir = "C:/NCI_NDR/noxn_dir_change_masked"
    if not os.path.exists(masked_filled_dir):
        os.makedirs(masked_filled_dir)
    # directory containing masked endpoints
    endpoint_dir = "C:/NCI_NDR/endpoints_dir_change_masked"
    if not os.path.exists(endpoint_dir):
        os.makedirs(endpoint_dir)

    # align n export and noxn rasters
    template_raster = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/subset_2000_2015/intermediate/template.tif"
    target_pixel_size = pygeoprocessing.get_raster_info(
        template_raster)['pixel_size']
    target_bb = pygeoprocessing.get_raster_info(
        template_raster)['bounding_box']
    full_scenario_list = [_N_EXPORT_BASELINE_KEY] + _N_EXPORT_PATH_LIST
    input_path_list = (
        [template_raster] + [
        os.path.join(nexport_dir, n_export_pattern.format(s)) for s in
        full_scenario_list] +
        [os.path.join(predicted_noxn_dir, 'surface_noxn_{}.tif'.format(s)) for
        s in full_scenario_list] +
        [os.path.join(predicted_noxn_dir, 'ground_noxn_{}.tif'.format(s)) for
        s in full_scenario_list])
    aligned_path_list = (
        [template_raster] +
        [os.path.join(aligned_dir, n_export_pattern.format(s)) for s in
        full_scenario_list] +
        [os.path.join(aligned_dir, 'surface_noxn_{}.tif'.format(s)) for
        s in full_scenario_list] +
        [os.path.join(aligned_dir, 'ground_noxn_{}.tif'.format(s)) for
        s in full_scenario_list])
    if not all([os.path.isfile(f) for f in aligned_path_list]):
        pygeoprocessing.align_and_resize_raster_stack(
            input_path_list, aligned_path_list,
            ['near'] * len(input_path_list), target_pixel_size, target_bb,
            raster_align_index=0)

    baseline_n_export_path = os.path.join(
        aligned_dir, n_export_pattern.format(_N_EXPORT_BASELINE_KEY))

    for fraction in ['surface', 'ground']:
        baseline_noxn_path = os.path.join(
            aligned_dir, '{}_noxn_{}.tif'.format(
                fraction, _N_EXPORT_BASELINE_KEY))
        # copy baseline noxn to masked_filled_dir
        shutil.copyfile(
            baseline_noxn_path, os.path.join(
                masked_filled_dir, '{}_noxn_{}.tif'.format(
                    fraction, _N_EXPORT_BASELINE_KEY)))
        for scenario_key in _N_EXPORT_PATH_LIST:
            scenario_noxn_path = os.path.join(
                aligned_dir, '{}_noxn_{}.tif'.format(
                    fraction, scenario_key))
            scenario_n_export_path = os.path.join(
                aligned_dir, n_export_pattern.format(scenario_key))
            mask_path = os.path.join(
                mask_dir, 'dir_change_mask_{}_{}.tif'.format(
                    fraction, scenario_key))
            masked_noxn_path = os.path.join(
                masked_filled_dir, '{}_noxn_{}.tif'.format(
                    fraction, scenario_key))
            if not os.path.exists(masked_noxn_path):
                calc_dir_change_masked_scenario_noxn(
                    scenario_noxn_path, baseline_noxn_path,
                    scenario_n_export_path, baseline_n_export_path,
                    mask_path, masked_noxn_path)

    global _PROCESSING_DIR
    _PROCESSING_DIR = "C:/NCI_NDR/intermediate"  # os.path.join(predicted_noxn_dir, 'intermediate')
    calc_endpoints(masked_filled_dir, endpoint_dir)


def masked_endpoints_workflow():
    """Generate endpoints for scenarios according to sig diff from baseline.

    For each scenario in _N_EXPORT_PATH_LIST, calculate estimated cancer cases
    and water treatment costs. Restrict predicted nitrate in surface and ground
    water for each scenario to be equal to baseline estimates, except in cases
    where predicted nitrate is significantly different from baseline according
    to standard error of the estimate.

    Side effects:
        creates the following rasters for each scenario:
            nitrate in drinking water
            cancer cases
            abatement costs

    Returns:
        None

    """
    # directory containing mean and se predictions
    predicted_noxn_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/R_ranger_pred"
    # directory containing noxn from which to calc endpoints
    sig_diff_filled_dir = "C:/NCI_NDR/sig_diff_filled_noxn"
    if not os.path.exists(sig_diff_filled_dir):
        os.makedirs(sig_diff_filled_dir)
    # directory containing scenario noxn only in areas significantly different
    # from baseline
    sig_diff_masked_dir = "C:/NCI_NDR/sig_diff_masked_scenario_noxn"
    if not os.path.exists(sig_diff_masked_dir):
        os.makedirs(sig_diff_masked_dir)
    # directory containing endpoints
    endpoint_dir = "C:/NCI_NDR/endpoints"
    if not os.path.exists(endpoint_dir):
        os.makedirs(endpoint_dir)
    for fraction in ['surface', 'ground']:
        baseline_noxn_path = os.path.join(
            predicted_noxn_dir, '{}_noxn_{}.tif'.format(
                fraction, _N_EXPORT_BASELINE_KEY))
        baseline_se_path = os.path.join(
            predicted_noxn_dir, '{}_noxn_se_{}.tif'.format(
                fraction, _N_EXPORT_BASELINE_KEY))
        # copy baseline noxn to sig_diff_filled_dir
        shutil.copyfile(
            baseline_noxn_path, os.path.join(
                sig_diff_filled_dir, '{}_noxn_{}.tif'.format(
                    fraction, _N_EXPORT_BASELINE_KEY)))
        for scenario_key in _N_EXPORT_PATH_LIST:
            if scenario_key == _N_EXPORT_BASELINE_KEY:
                continue  # no need to calculate masked scenario for baseline
            scenario_noxn_path = os.path.join(
                predicted_noxn_dir, '{}_noxn_{}.tif'.format(
                    fraction, scenario_key))
            scenario_se_path = os.path.join(
                predicted_noxn_dir,
                '{}_noxn_se_{}.tif'.format(fraction, scenario_key))
            sig_diff_scenario_path = os.path.join(
                sig_diff_masked_dir, '{}_noxn_{}.tif'.format(
                    fraction, scenario_key))
            sig_diff_filled_path = os.path.join(
                sig_diff_filled_dir, '{}_noxn_{}.tif'.format(
                    fraction, scenario_key))
            if not os.path.exists(sig_diff_filled_path):
                calc_masked_scenario_noxn(
                    scenario_noxn_path, scenario_se_path, baseline_noxn_path,
                    baseline_se_path, sig_diff_scenario_path,
                    sig_diff_filled_path)
    global _PROCESSING_DIR
    _PROCESSING_DIR = "C:/NCI_NDR/intermediate"  # os.path.join(predicted_noxn_dir, 'intermediate')
    calc_endpoints(sig_diff_filled_dir, endpoint_dir)


def zonal_stat_data_frame(
        raster_path, aggregate_vector_path, fid_field):
    """Calculate zonal sum and mean and return a data frame.

    Parameters:
        raster_path (string): path to raster containing values that should be
            summarized via zonal stats
        aggregate_vector_path (string): a path to a polygon vector containing
            zones to summarize values within
        fid_field (string): field in aggregate_vector_path that identifies
            features

    Returns:
        a data frame containing zonal stats: sum and mean

    """
    fid_to_objectid = map_FID_to_field(aggregate_vector_path, fid_field)
    zonal_stat_dict = pygeoprocessing.zonal_statistics(
        (raster_path, 1), aggregate_vector_path, polygons_might_overlap=False)
    objectid_zonal_stats_dict = {
        objectid: zonal_stat_dict[fid] for (fid, objectid) in
        fid_to_objectid.items()
    }
    zonal_df = pandas.DataFrame(
        {
            fid_field: [
                key for key, value in sorted(
                    objectid_zonal_stats_dict.items())],
            'sum': [
                value['sum'] for key, value in
                sorted(objectid_zonal_stats_dict.items())],
            'count': [
                value['count'] for key, value in
                sorted(objectid_zonal_stats_dict.items())]
        })
    zonal_df['mean'] = zonal_df['sum'] / zonal_df['count']
    return zonal_df


def combined_zonal_stats():
    """Calculate zonal stats by country for inputs and endpoints."""
    # scenarios for which to contrast zonal stats
    scenario_list = [_N_EXPORT_BASELINE_KEY] + _N_EXPORT_PATH_LIST
    countries_shp_path = "F:/NCI_NDR/Data world borders/TM_WORLD_BORDERS-0.3.shp"
    fid_field = 'NAME'  # 'ISO3'

    filled_noxn_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/filter_by_direction_of_change/noxn_dir_change_masked"
    raw_noxn_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/R_ranger_pred"
    endpoints_masked_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/filter_by_direction_of_change/endpoints_dir_change_masked"
    mask_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/filter_by_direction_of_change/filter_mask"
    rescaled_endpoint_dir = "F:/NCI_NDR/Results_backup/Results_3.27.21/filter_by_direction_of_change/endpoints_dir_change_masked/masked_protected_areas"
    df_list = []

    # N export by country: mean
    # n_export_pattern = "F:/NCI_NDR/Data NDR/updated_3.27.21/resampled_by_Becky/renamed/compressed_<scenario>.tif"
    # for scenario_id in scenario_list:
    #     n_export_path = n_export_pattern.replace('<scenario>', scenario_id)
    #     zonal_df = zonal_stat_data_frame(
    #         n_export_path, countries_shp_path, fid_field)
    #     mean_df = zonal_df[[fid_field, 'mean']]
    #     mean_df.rename(
    #         index=str,
    #         columns={'mean': 'n_export_mean_{}'.format(scenario_id)},
    #         inplace=True)
    #     df_list.append(mean_df)

    # nitrate in groundwater by country, masked: mean
    masked_ground_pattern = os.path.join(
        rescaled_endpoint_dir, "ground_noxn_<scenario>.tif")
    for scenario_id in scenario_list:
        gr_noxn_path = masked_ground_pattern.replace('<scenario>', scenario_id)
        zonal_df = zonal_stat_data_frame(
            gr_noxn_path, countries_shp_path, fid_field)
        mean_df = zonal_df[[fid_field, 'mean']]
        mean_df.rename(
            index=str,
            columns={
                'mean': 'ground_noxn_masked_rescaled_mean_{}'.format(
                    scenario_id)},
            inplace=True)
        df_list.append(mean_df)

    # nitrate in surfacewater by country, masked: mean
    masked_surf_pattern = os.path.join(
        rescaled_endpoint_dir, "surface_noxn_<scenario>.tif")
    for scenario_id in scenario_list:
        surf_noxn_path = masked_surf_pattern.replace(
            '<scenario>', scenario_id)
        zonal_df = zonal_stat_data_frame(
            surf_noxn_path, countries_shp_path, fid_field)
        mean_df = zonal_df[[fid_field, 'mean']]
        mean_df.rename(
            index=str,
            columns={
                'mean': 'surf_noxn_masked_rescaled_mean_{}'.format(
                    scenario_id)},
            inplace=True)
        df_list.append(mean_df)

    # nitrate in groundwater by country, not masked: mean
    # ground_pattern = os.path.join(raw_noxn_dir, "ground_noxn_<scenario>.tif")
    # for scenario_id in scenario_list:
    #     gr_noxn_path = ground_pattern.replace('<scenario>', scenario_id)
    #     zonal_df = zonal_stat_data_frame(
    #         gr_noxn_path, countries_shp_path, fid_field)
    #     mean_df = zonal_df[[fid_field, 'mean']]
    #     mean_df.rename(
    #         index=str,
    #         columns={'mean': 'ground_noxn_unmasked_mean_{}'.format(
    #             scenario_id)},
    #         inplace=True)
    #     df_list.append(mean_df)

    # # nitrate in surfacewater by country, not masked: mean
    # surf_pattern = os.path.join(
    #     raw_noxn_dir, "surface_noxn_<scenario>.tif")
    # for scenario_id in scenario_list:
    #     surf_noxn_path = surf_pattern.replace(
    #         '<scenario>', scenario_id)
    #     zonal_df = zonal_stat_data_frame(
    #         surf_noxn_path, countries_shp_path, fid_field)
    #     mean_df = zonal_df[[fid_field, 'mean']]
    #     mean_df.rename(
    #         index=str,
    #         columns={'mean': 'surf_noxn_unmasked_mean_{}'.format(scenario_id)},
    #         inplace=True)
    #     df_list.append(mean_df)

    # groundwater mask: % of pixels in country that are masked (i.e., direction
    # of change in predicted noxn is "wrong")
    # ground_mask_pattern = os.path.join(
    #     mask_dir, 'dir_change_mask_ground_{}.tif')
    # for scenario_id in scenario_list:
    #     if scenario_id == _N_EXPORT_BASELINE_KEY:
    #         continue  # no mask calculated for baseline
    #     ground_mask_path = ground_mask_pattern.format(scenario_id)
    #     zonal_df = zonal_stat_data_frame(
    #         ground_mask_path, countries_shp_path, fid_field)
    #     mean_df = zonal_df[[fid_field, 'mean']]
    #     mean_df.rename(
    #         index=str,
    #         columns={'mean': 'dir_change_%_masked_ground_{}'.format(
    #             scenario_id)},
    #         inplace=True)
    #     df_list.append(mean_df)

    # # surface mask: % of pixels in country that are masked (i.e., direction
    # # of change in predicted noxn is "wrong")
    # surf_mask_pattern = os.path.join(
    #     mask_dir, 'dir_change_mask_surface_{}.tif')
    # for scenario_id in scenario_list:
    #     if scenario_id == _N_EXPORT_BASELINE_KEY:
    #         continue  # no mask calculated for baseline
    #     surf_mask_path = surf_mask_pattern.format(scenario_id)
    #     zonal_df = zonal_stat_data_frame(
    #         surf_mask_path, countries_shp_path, fid_field)
    #     mean_df = zonal_df[[fid_field, 'mean']]
    #     mean_df.rename(
    #         index=str,
    #         columns={'mean': 'dir_change_%_masked_surface_{}'.format(
    #             scenario_id)},
    #         inplace=True)
    #     df_list.append(mean_df)

    # # cancer cases, not masked: sum
    # cases_pattern = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/endpoints_not_masked/cancer_cases_<scenario>.tif"
    # for scenario_id in scenario_list:
    #     cases_path = cases_pattern.replace(
    #         '<scenario>', scenario_id)
    #     zonal_df = zonal_stat_data_frame(
    #         cases_path, countries_shp_path, fid_field)
    #     sum_df = zonal_df[[fid_field, 'sum']]
    #     sum_df.rename(
    #         index=str,
    #         columns={'sum': 'cancer_cases_unmasked_sum_{}'.format(
    #             scenario_id)},
    #         inplace=True)
    #     df_list.append(sum_df)

    # cancer cases, masked
    cases_pattern = os.path.join(
        rescaled_endpoint_dir, "cancer_cases_<scenario>.tif")
    for scenario_id in scenario_list:
        cases_path = cases_pattern.replace('<scenario>', scenario_id)
        zonal_df = zonal_stat_data_frame(
            cases_path, countries_shp_path, fid_field)
        mean_df = zonal_df[[fid_field, 'mean']]
        mean_df.rename(
            index=str,
            columns={
                'mean': 'noxn_drinking_water_masked_rescaled_mean_{}'.format(
                    scenario_id)},
            inplace=True)
        df_list.append(mean_df)

    # drinking water: mean
    drink_pattern = os.path.join(
        rescaled_endpoint_dir, 'noxn_in_drinking_water_<scenario>.tif')
    for scenario_id in scenario_list:
        drink_path = drink_pattern.replace('<scenario>', scenario_id)
        zonal_df = zonal_stat_data_frame(
            drink_path, countries_shp_path, fid_field)
        sum_df = zonal_df[[fid_field, 'sum']]
        sum_df.rename(
            index=str,
            columns={'sum': 'cancer_cases_masked_rescaled_sum_{}'.format(
                scenario_id)},
            inplace=True)
        df_list.append(sum_df)

    # merge data frames together
    combined_df_path = "F:/NCI_NDR/Results_backup/Results_3.27.21/filter_by_direction_of_change/zonal_statistics_rescaled_mosaicked_summary.csv"
    combined_df = df_list[0]
    df_i = 1
    while df_i < len(df_list):
        combined_df = combined_df.merge(
            df_list[df_i], on=fid_field, suffixes=(False, False),
            validate="one_to_one")
        df_i = df_i + 1
    transposed_df = combined_df.transpose()
    try:
        transposed_df.to_csv(combined_df_path)
    except PermissionError:
        import pdb; pdb.set_trace()


def check_order_of_scenarios(endpoint_dir):
    """Check that cancer case scenarios share same order in each country.

    Parameters:
        endpoint_dir (string): directory containing rasters of cancer cases

    Side effects:
        writes two tables to endpoint_dir

    Returns:
        None

    """
    countries_shp_path = "F:/NCI_NDR/Data world borders/TM_WORLD_BORDERS-0.3.shp"
    fid_field = 'ISO3'
    cases_pattern = os.path.join(
        endpoint_dir, "cancer_cases_<scenario>.tif")

    df_list = []
    for scenario_id in _N_EXPORT_PATH_LIST:
        cases_path = cases_pattern.replace(
            '<scenario>', scenario_id)
        zonal_df = zonal_stat_data_frame(
            cases_path, countries_shp_path, fid_field)
        sum_df = zonal_df[[fid_field, 'sum']]
        sum_df.set_index(fid_field, inplace=True)
        sum_df.rename(index=str, columns={'sum': scenario_id}, inplace=True)
        df_list.append(sum_df)
    combined_df = df_list[0]
    df_i = 1
    while df_i < len(df_list):
        combined_df = combined_df.merge(
            df_list[df_i], on=fid_field, suffixes=(False, False),
            validate="one_to_one")
        df_i = df_i + 1
    combined_df.to_csv(
        os.path.join(endpoint_dir, "cancer_cases_by_country.csv"))
    combined_df_t = combined_df.transpose()
    # drop countries with 0 cancer cases in any scenario
    combined_df_t.drop(
        [col for col, val in combined_df_t.sum().iteritems() if val == 0],
        axis=1, inplace=True)
    order_dict = {}
    for label, content in combined_df_t.items():
        order_dict[label] = content.sort_values().index
    order_df = pandas.DataFrame(order_dict)
    order_df.to_csv(
        os.path.join(endpoint_dir, "scenario_order_by_country.csv"))


def resize_endpoint_rasters(
        endpoint_dir, rescaled_endpoint_dir, div=True):
    """Rescale endpoints to match ESA resolution.

    Optionally scale pixel values by pixel area.

    Parameters:
        endpoint_dir (string): path to directory containing endpoint
            rasters at ~10 km resolution
        rescaled_endpoint_dir (string): path to directory where rescaled
            endpoint rasters matching the resolution of ESA landcover
            should be created
        div (bool): should endpoint rasters be scaled by pixel size? This value
            should be false if the pixel value of endpoint rasters is
            independent of pixel area

    Side effects:
        creates a copy of each endpoint raster in `endpoint_dir` which is
            aligned with ESA landcover

    """
    def divide_op(value_raster, divisor):
        """Divide values in value raster by divisor."""
        valid_mask = (value_raster != _TARGET_NODATA)
        result = numpy.empty(value_raster.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        result[valid_mask] = value_raster[valid_mask] / float(divisor)
        return result

    if not os.path.exists(rescaled_endpoint_dir):
        os.makedirs(rescaled_endpoint_dir)
    ESA_path = "F:/KBA_ES_archive/ESA_landcover_2015/product/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif"
    target_pixel_size = pygeoprocessing.get_raster_info(
        ESA_path)['pixel_size']

    input_path_list = [
        os.path.join(endpoint_dir, f) for f in os.listdir(endpoint_dir) if
        os.path.isfile(os.path.join(endpoint_dir, f))]
    input_pixel_size = pygeoprocessing.get_raster_info(
        input_path_list[0])['pixel_size']
    divisor = (input_pixel_size[0] / target_pixel_size[0]) ** 2

    # store divided rasters in temporary directory
    intermediate_dir = tempfile.mkdtemp()
    if div:
        paths_to_align = []
        for in_path in input_path_list:
            target_path = os.path.join(intermediate_dir, os.path.basename(in_path))
            print(
                "Divide coarse scale raster: {}".format(os.path.basename(in_path)))
            pygeoprocessing.raster_calculator(
                [(in_path, 1), (divisor, 'raw')], divide_op, target_path,
                gdal.GDT_Float32, _TARGET_NODATA)
            paths_to_align.append(target_path)
    else:
        paths_to_align = input_path_list

    # align rescaled rasters with ESA landcover
    source_input_path_list = [ESA_path] + paths_to_align
    source_bn_list = [os.path.basename(f) for f in source_input_path_list]
    aligned_path_list = [
        os.path.join(
            rescaled_endpoint_dir, f) for f in source_bn_list]
    print("Align and resize rescaled rasters")
    pygeoprocessing.align_and_resize_raster_stack(
        source_input_path_list, aligned_path_list,
        ['near'] * len(source_input_path_list),
        target_pixel_size, 'union', raster_align_index=0)

    # clean up
    shutil.rmtree(intermediate_dir)


def mosaic_protected_areas(rescaled_endpoint_dir, masked_protected_areas_dir):
    """Mosaic cancer cases in protected areas from restoration scenario.

    Here we are patching a missing link in the scenario creation workflow.
    Protected areas were not treated correctly in the scenarios that went into
    running NDR, so instead we will mosaic final results in to
    every scenario from the restoration scenario, inside areas defined by
    protected areas mask.

    Parameters:
        rescaled_endpoint_dir (string): path to directory containing rescaled
            endpoint rasters matching the resolution of ESA landcover
        masked_protected_areas_dir (string): path to directory where masked
            endpoint rasters should be created

    Side effects:
        creates a copy of each raster in
            `masked_protected_areas_dir` which has endpoint values from
            restoration scenario inside protected areas

    """
    def mosaic_op(target_ar, restoration_ar, *mask_list):
        """Mosaic values according to a stack of masks.

        Where any raster in `mask_list` is 1, mosaic the value from
        restoration_ar into target_ar.

        """
        mosaic_mask = numpy.any(
            numpy.isclose(numpy.array(mask_list), 1), axis=0)
        result = numpy.copy(target_ar)
        result[mosaic_mask] = restoration_ar[mosaic_mask]
        return result

    pa_dir = "F:/NCI_NDR/Data protected areas"
    # applied to all layers
    all_mask_path = os.path.join(
        pa_dir, 'aligned_to_esa', 'wdpa_iucn_cat_i-iv.tif')
    # applied to all expansion layers
    expansion_mask_path = os.path.join(
        pa_dir,
        'wdpa_iucn_cat_v_full_md5_590504ff12a06f2ca3f080121d655cc4.tif')
    # applied to all layers except grazing
    allbutgrazing_mask_path = os.path.join(
        pa_dir, 'aligned_to_esa', 'wdpa_iucn_cat_vi.tif')

    if not os.path.exists(masked_protected_areas_dir):
        os.makedirs(masked_protected_areas_dir)

    # Mosaic protected areas in to all the endpoints
    for endpoint in [
            'noxn_in_drinking_water', 'ground_noxn', 'surface_noxn']:  # , 'cancer_cases'
        print("mosaic protected areas: {}, all scenarios \n".format(endpoint))
        restoration_path = os.path.join(
            rescaled_endpoint_dir, '{}_restoration.tif'.format(endpoint))

        # all three masks
        three_mask_list = [
            'extensification_bmps_irrigated', 'extensification_bmps_rainfed',
            'extensification_current_practices',
            'extensification_intensified_irrigated',
            'extensification_intensified_rainfed']
        mask_list = [
            all_mask_path, expansion_mask_path, allbutgrazing_mask_path]
        for scenario_key in three_mask_list:
            target_path = os.path.join(
                rescaled_endpoint_dir,
                '{}_{}.tif'.format(endpoint, scenario_key))
            mosaic_path = os.path.join(
                masked_protected_areas_dir,
                '{}_{}.tif'.format(endpoint, scenario_key))
            if not os.path.exists(mosaic_path):
                pygeoprocessing.raster_calculator(
                    [(path, 1) for path in [target_path, restoration_path] +
                    mask_list], mosaic_op, mosaic_path, gdal.GDT_Float32,
                    _TARGET_NODATA)

        two_mask_list = [
            'baseline_currentpractices', 'fixedarea_bmps_irrigated',
            'fixedarea_bmps_rainfed', 'fixedarea_intensified_irrigated',
            'fixedarea_intensified_rainfed', 'sustainable_currentpractices']
        mask_list = [all_mask_path, allbutgrazing_mask_path]
        for scenario_key in two_mask_list:
            target_path = os.path.join(
                rescaled_endpoint_dir,
                '{}_{}.tif'.format(endpoint, scenario_key))
            mosaic_path = os.path.join(
                masked_protected_areas_dir,
                '{}_{}.tif'.format(endpoint, scenario_key))
            if not os.path.exists(mosaic_path):
                pygeoprocessing.raster_calculator(
                    [(path, 1) for path in [target_path, restoration_path] +
                    mask_list], mosaic_op, mosaic_path, gdal.GDT_Float32,
                    _TARGET_NODATA)

        scenario_key = 'grazing_expansion'
        mask_list = [all_mask_path, expansion_mask_path]
        target_path = os.path.join(
            rescaled_endpoint_dir,
            '{}_{}.tif'.format(endpoint, scenario_key))
        mosaic_path = os.path.join(
            masked_protected_areas_dir,
            '{}_{}.tif'.format(endpoint, scenario_key))
        if not os.path.exists(mosaic_path):
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in [target_path, restoration_path] +
                mask_list], mosaic_op, mosaic_path, gdal.GDT_Float32,
                _TARGET_NODATA)


        # copy restoration scenario to mosaic dir, to be shared with Peter
        shutil.copyfile(
            restoration_path, os.path.join(
                masked_protected_areas_dir,
                '{}_restoration.tif'.format(endpoint)))


def main():
    """Program entry point."""
    # predict_endpoints_only()
    # confidence_interval_wrapper()
    # ground_ci_wrapper()
    # prepare_covariates_May2020_NDR()
    # masked_endpoints_workflow()
    # endpoint_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/endpoints"
    # endpoint_dir = "C:/NCI_NDR/endpoints"
    # check_order_of_scenarios(endpoint_dir)
    # rescaled_endpoint_dir = os.path.join(
        # endpoint_dir, 'rescaled_ESA_resolution')
    # resize_endpoint_rasters(endpoint_dir, rescaled_endpoint_dir)
    # masked_protected_areas_dir = os.path.join(
        # endpoint_dir, 'masked_protected_areas')
    # mosaic_protected_areas(rescaled_endpoint_dir, masked_protected_areas_dir)
    # endpoint_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_5.15.20/endpoints_not_masked"
    # rescaled_endpoint_dir = os.path.join(
    #     endpoint_dir, 'rescaled_ESA_resolution')
    # resize_endpoint_rasters(endpoint_dir, rescaled_endpoint_dir)
    # check_order_of_scenarios(rescaled_endpoint_dir)
    # prepare_covariates_March2021_NDR()
    # calculate endpoints not masked
    global _PROCESSING_DIR
    _PROCESSING_DIR = "C:/NCI_NDR/intermediate"
    # noxn_dir = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/Results_3.27.21/R_ranger_pred"
    # endpoint_dir = "C:/NCI_NDR/endpoints_not_masked"
    # calc_endpoints(noxn_dir, endpoint_dir)
    endpoint_dir = "F:/NCI_NDR/Results_backup/Results_3.27.21/filter_by_direction_of_change/endpoints_dir_change_masked"
    # dir_change_mask_endpoints_workflow()
    rescaled_endpoint_dir = os.path.join(
        endpoint_dir, 'rescaled_ESA_resolution')
    # resize_endpoint_rasters(endpoint_dir, rescaled_endpoint_dir)
    # resize noxn in ground and surface water too
    filtered_noxn_dir = "F:/NCI_NDR/Results_backup/Results_3.27.21/filter_by_direction_of_change/noxn_dir_change_masked"
    # resize_endpoint_rasters(
        # filtered_noxn_dir, rescaled_endpoint_dir, div=False)
    masked_protected_areas_dir = os.path.join(
        endpoint_dir, 'masked_protected_areas')
    # mosaic_protected_areas(rescaled_endpoint_dir, masked_protected_areas_dir)
    combined_zonal_stats()

if __name__ == '__main__':
    __spec__ = None  # for running with pdb
    main()
