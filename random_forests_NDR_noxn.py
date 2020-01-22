"""Train random forests model to predict noxn from NDR outputs."""
import os
import argparse

import numpy
import pandas
from sklearn.ensemble import RandomForestRegressor


def random_forests_NDR_noxn(noxn_predictor_df_path):
    """Train and test random forests model.

    Parameters:
        noxn_predictor_df_path (string): path to csv containing noxn
            obserations and predictor values for a series of points. This csv
            must also include a column of boolean values, named "train",
            indicating which rows should be included in a training data subset.
            The complement of this subset is used for testing.

    Returns:
        None

    """
    # read data that was filtered and subsetted in R
    combined_df = pandas.read_csv(noxn_predictor_df_path)
    train_mask = numpy.array(combined_df.train) == 1
    combined_df.drop('train', axis=1)

    # generate dummy variables for factors
    combined_df = pandas.get_dummies(combined_df)

    # Separate data frame into response and predictors
    noxn_arr = numpy.array(combined_df['noxn'])
    predictor_df = combined_df.drop('noxn', axis=1)
    predictor_arr = numpy.array(predictor_df)

    # Saving feature names for later use
    predictor_names = list(predictor_df.columns)

    # separate into training and test subsets
    noxn_train_arr = noxn_arr[train_mask]
    noxn_test_arr = noxn_arr[~train_mask]

    predictor_train_arr = predictor_arr[train_mask]
    predictor_test_arr = predictor_arr[~train_mask]

    # train random forests model
    # Instantiate model with 1000 decision trees
    rf_model = RandomForestRegressor(n_estimators=1000, random_state=42)

    # train the model on training data
    rf_model.fit(predictor_train_arr, noxn_train_arr)

    # compute R squared of predictions on test data set
    test_r_sq = rf_model.score(predictor_test_arr, noxn_test_arr)

    # if needed for comparison with caret results: get predicted noxn for test
    # subset from trained model
    noxn_pred = rf_model.predict(predictor_test_arr)

    # get variable importance
    var_importances = rf_model.feature_importances_

    # print outputs to console
    print(
        "R squared (test partition, predicted vs actual): {}".format(test_r_sq))
    print("Variable importances: {}".format(var_importances))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Random forests model.')
    parser.add_argument(
        'noxn_predictor_df_path', type=str, nargs='+',
        help='path to csv file containing noxn and predictor covariates')
    args = parser.parse_args()
    noxn_predictor_df_path = args.noxn_predictor_df_path[0]
    random_forests_NDR_noxn(noxn_predictor_df_path)
