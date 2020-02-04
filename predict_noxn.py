"""Use random forests to predict nitrate concentrations globally."""
import numpy
import pandas

from sklearn.ensemble import RandomForestRegressor
from joblib import dump, load

# covariate datasets at native resolution
_INPUT_COVARIATE_PATH_DICT = {}

# aligned covariates
_ALIGNED_COVARIATE_PATH_DICT = {}

# noxn and covariate observations for groundwater
_NOXN_PREDICTOR_GR_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_gr.csv"

# noxn and covariate observations for surface water
_NOXN_PREDICTOR_SURF_DF_PATH = "C:/Users/ginge/Documents/Python/nci_ndr/noxn_predictor_df_surf.csv"

def prepare_covariates():
    """Align covariate rasters."""
    # TODO align and resample covariates in _INPUT_COVARIATE_PATH_DICT,
    # populating _ALIGNED_COVARIATE_PATH_DICT


def predict_noxn(
        n_export_path, rf_pickle_filename, predictor_names, output_path):
    """Generate global nitrate concentrations from a set of predictors.

    Parameters:
        n_export_path (string): path to raster containing N export for a single
            scenario
        rf_pickle_filename (string): path to file on disk containing trained
            random forests model. This file should have the extension 'joblib'
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
    # align n_export_path with covariates if necessary
    rf_model = load(rf_pickle_filename)

    # load covariates in a block that can be contained in memory
    # must be in the same order as predictor_names
    # make predictions for a single block
    # noxn_predicted = rf_model.predict(covar_arr)

    # save predictions for that block


def train_rf_model(
        noxn_predictor_df_path, max_features, rf_pickle_filename,
        factor_list=None):
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
        factor_list (list): optional list of strings indicating variables in
            `noxn_predictor_df_path` that represent factors

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

    # generate dummy variables for factors
    if factor_list:
        dummy_df_list = []
        for var in factor_list:
            dummy_df_list.append(
                pandas.get_dummies(combined_df[var], prefix=var))
            combined_df.drop([var], axis=1, inplace=True)
        dummy_df_list.append(combined_df)
        combined_df = pandas.concat(dummy_df_list, axis=1)

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


def main():
    """Program entry point."""
    pass


if __name__ == '__main__':
    main()
