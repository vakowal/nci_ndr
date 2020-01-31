"""Train random forests model to predict noxn from NDR outputs."""
import os
import argparse

import numpy
import pandas
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, KFold


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
    # drop rows containing missing data
    combined_df.dropna(inplace=True)
    combined_df.drop('train', axis=1, inplace=True)

    # generate dummy variables for factors
    lake_dummy = pandas.get_dummies(combined_df['lake'], prefix='lake')
    river_dummy = pandas.get_dummies(combined_df['river'], prefix='river')
    combined_df.drop(['lake', 'river'], axis=1, inplace=True)
    combined_df = pandas.concat([combined_df, lake_dummy, river_dummy], axis=1)

    # Separate data frame into response and predictors
    noxn_arr = numpy.array(combined_df['noxn'])
    predictor_df = combined_df.drop('noxn', axis=1)
    predictor_arr = numpy.array(predictor_df)

    # Saving predictor names for later use
    predictor_names = list(predictor_df.columns)

    # train random forests model
    # Instantiate model, using parameters chosen to match those used in R
    rf_model = RandomForestRegressor(
        random_state=42, n_estimators=500, criterion="mse", max_features=10,
        min_samples_leaf=5, oob_score=True, bootstrap=True)

    # evaluate the model working on all data
    rf_model.fit(predictor_arr, noxn_arr)
    training_rsq = rf_model.score(predictor_arr, noxn_arr)
    print('R^2 Training Score: {}'.format(training_rsq))
    print('OOB Score: {}'.format(rf_model.oob_score_))

    # 10-fold cross-validation
    rf_cv_scores = cross_val_score(
        rf_model, predictor_arr, noxn_arr,
        cv=KFold(10, shuffle=True, random_state=42))
    mean_cv_score = numpy.asarray(rf_cv_scores).mean()
    print('mean cross validation R^2: {}'.format(mean_cv_score))

    # get variable importance
    var_importances = rf_model.feature_importances_
    var_imp_df = pandas.DataFrame(
        {'variable': predictor_names, 'importance': var_importances})
    print(var_imp_df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Random forests model.')
    parser.add_argument(
        'noxn_predictor_df_path', type=str, nargs='+',
        help='path to csv file containing noxn and predictor covariates')
    args = parser.parse_args()
    noxn_predictor_df_path = args.noxn_predictor_df_path[0]
    random_forests_NDR_noxn(noxn_predictor_df_path)
