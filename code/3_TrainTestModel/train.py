#!/usr/bin/env python
# ==============================================================================

# Name: train.py
# Author: Braden Yang
# Created: 03/13/2024
# Description: train an ensemble of classifiers to predict binary CN to MCI or
#   AD progression

# ==============================================================================

# *** toggle interactive mode ***
INTERACTIVE = False
# *** toggle interactive mode ***

# ==========================
# ===== IMPORT MODULES =====
# ==========================

import argparse
import os
import sys

import joblib
import numpy as np
from sklearn.model_selection import StratifiedKFold, LeaveOneGroupOut, GridSearchCV
from sklearn.metrics import roc_auc_score, accuracy_score, balanced_accuracy_score, f1_score
from sklearn.pipeline import Pipeline

sys.path.append("/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/code")
from project_utils import *


# ============================
# ===== DEFINE FUNCTIONS =====
# ============================


def get_metrics(model, X, y_true):

    if hasattr(model, "best_estimator_"):
        model = model.best_estimator_

    # get prediction scores
    if isinstance(model, RandomSubsetEnsembleClassifier):  # handle cases of custom RandomSubsetEnsembleClassifier
        if hasattr(model.base_estimator, "predict_proba"):
            y_score = model.predict_proba(X)[:,1]
        else:
            y_score = model.decision_function(X)
    elif hasattr(model, "predict_proba"):
        y_score = model.predict_proba(X)[:,1]
    else:
        y_score = model.decision_function(X)

    # get binary predictions
    y_pred = model.predict(X)

    # get confusion matrix
    conf_mat = confusion_matrix(y_true=y_true, y_pred=y_pred)
    if conf_mat.shape[0] != 2 or conf_mat.shape[1] != 2:
        TN, FN, FP, TP = np.nan, np.nan, np.nan, np.nan
    else:
        TN = conf_mat[0,0]; FN = conf_mat[1,0]; FP = conf_mat[0,1]; TP = conf_mat[1,1]

    # return accuracy metrics
    return {
        "roc_auc": roc_auc_score(y_true, y_score),
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "f1": f1_score(y_true, y_pred),
        "sensitivity": (TP) / (TP + FN),
        "specificity": (TN) / (TN + FP),
        "positive_predictive_value": (TP) / (TP + FP),
        "negative_predictive_value": (TN) / (TN + FN),
    }

def chen_2023_haufe_weights(train_X, y_pred):
    # Code is taken from JianZhong Chen's 2023 NeuroImage paper
    # https://doi.org/10.1016/j.neuroimage.2023.120115
    # https://github.com/ThomasYeoLab/CBIG/blob/424e2d9633521b2521d46e4bf22efa74d894f29f/stable_projects/predict_phenotypes/ChenOoi2023_ICCW/regression/RF/CBIG_ICCW_RF_interpretation.py#L276-L283
    # There, it was used to estimate Haufe feature importances
    # for nonlinear RF models.

    # The code is used mostly as is, except we don't redo
    # demeaning of X_train (this is already done for the model training)

    # Took me some time to understand this code, but think I get it now.
    # Covariance formula is as follows:
    #   cov(X, Y) = E[(X-E[X]) * (Y-E[Y])]

    # X-E[X] is the variable demeaned, which is achieved in the first two
    # lines of code for the training data and targets

    # The outer expectation is then the *mean* of the product of the
    # demeaned X and Y.  The last two lines achieve this, by
    # taking the dot product of each feature with Y and then dividing
    # by the number of observations (train_X.shape[0]).  The dot product
    # sums all the products, and the division takes the average.

    # Basically, we get the covariance between each feature and the model
    # predictions.  In Haufe 2014 (https://doi.org/10.1016/j.neuroimage.2013.10.067),
    # this is described in the "Simplifying conditions" section.  We can
    # take advantage of this since we only have one target (K=1).

    demean_y_pred = y_pred - np.mean(y_pred)
    demean_X_train = train_X.T

    fivals_Haufe = np.dot(demean_X_train, demean_y_pred) / train_X.shape[0]
    fivals_Haufe = fivals_Haufe[:, np.newaxis]

    return fivals_Haufe

def compute_haufe_ensemble(ensemble_model, X_train, approx = False):
    haufe_list = []
    for m in ensemble_model.estimators_:
        X_transform = m[:-1].transform(X_train)
        if approx:
            y_pred = m.predict(X_train)
        else:
            y_pred = m.decision_function(X_train)
        haufe_list.append(chen_2023_haufe_weights(X_transform, y_pred.flatten()))
    return np.concat(haufe_list, axis = 1)


# ======================
# ===== PARSE ARGS =====
# ======================


parser = argparse.ArgumentParser(
    description="Evaluate several ML models to predict binary CN to MCI/AD conversion on PAC data using nested cross validation",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument("-w", "--wdir", help="path to working directory")
parser.add_argument("-o", "--odir", help="path to output directory")
parser.add_argument("--prefix", help="prefix of output .joblib file")
parser.add_argument("-f", "--features", type=str, default="avn", help="features to include (a=amyloid, v=volume, n=nonimaging)")
parser.add_argument("-m", "--model", default="svm", type=str, help="name of model to train")

parser.add_argument("-t", "--max_time_to_progression", default=3, type=float, help="time window length for classifying CN-to-MCI progressors")
parser.add_argument("-s", "--min_time_stable", default=3, type=float, help="minimum time of stability for classifying CN-to-MCI stables")
parser.add_argument("--progressor_can_be_stable", action="store_true", help="if True, progressors that don't convert within max_time_to_progression but stay stable for min_time_stable are labeled as stable")

parser.add_argument("--random_state", type=int, default=None, help="set random seed")
parser.add_argument("--feature_reduction", type=str, default="none", choices=["relief", "pca", "none"], help="feature reduction technique, either Relief, PCA, or none (default=none)")
parser.add_argument("--control_standardize", action="store_true", help="use standard scaling with control group")
parser.add_argument("--combat", action="store_true", help="apply GAM-ComBat to harmonize tracers (only useful for leave-one-site-out validation)")
parser.add_argument("--ensemble", action="store_true", help="for logistic regression and SVM, create an ensemble classifier")
parser.add_argument("--n_estimators", type=int, default=100, help="number of estimators for ensemble classifier")
parser.add_argument("--leave_out_group", type=str, default="site", choices=["site", "tracer"], help="group to leave out for training/testing split")
parser.add_argument("--grid_search_cv_method", type=str, default="stratified_kfold", choices=["stratified_kfold", "leave_one_group_out"], help="cross-validation method for grid search; either 'stratified_kfold' or 'leave_one_group_out' is allowed")
parser.add_argument("--centiloid", action="store_true", help="use Centiloid-harmonized regional amyloid PET features instead of SUVRs")
parser.add_argument("--single_tracer", type=str, default=None, choices=["FBP", "PIB"], help="if specified, only use data from this tracer (e.g. FBP or PIB)")
parser.add_argument("--alt_stable", action="store_true", help="if True, use alternative stable definition where stables are allowed to progress beyond max_tim_to_progression but must remain stable for min_time_stable")

if INTERACTIVE:
    args = parser.parse_args(args = [])
    # define preset arguments
    args.wdir = "/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction"
    args.odir = "/scratch/b.y.yang/3_TrainTestModel/test"
    args.prefix = "model"
    args.model = "svm"
    args.features = "avn"
    
    args.max_time_to_progression = 3
    args.min_time_stable = 5

    args.control_standardize = True
    args.ensemble = True
    args.combat = True
    args.centiloid = False

    args.n_estimators = 2

    args.leave_out_group = "site"
else:
    args = parser.parse_args()

if not os.path.exists(args.odir): os.makedirs(args.odir, exist_ok = True)  # make output directory

# check arguments
if not args.model in ["lr", "svm", "rf"]: raise ValueError("invalid --model input")

# print args
print_args(args)


# ============================
# ===== DEFINE VARIABLES =====
# ============================


# define parameters used to train models
args_dict = vars(args)
param_dict = {
    "cv_splits_gridsearch": 5,
    "relief_inner_n": 20,
    "relief_outer_n": 20,
    "pca_variance_keep": 0.9,
    # "random_state": 42,
}
param_dict = param_dict | args_dict

# set random seed
np.random.seed(args.random_state)


# =====================
# ===== LOAD DATA =====
# =====================


print("+++++ Loading data +++++")
if args.centiloid:
    features_df = pd.read_csv(os.path.join(args.wdir, "data/features/features_cross_sectional_CL.csv"))
    features_df = features_df[features_df["site"].isin(["A4", "ADNI", "MCSA", "OASIS"])]
else:
    features_df = pd.read_csv(os.path.join(args.wdir, "data/features/features_cross_sectional.csv"))


# ===========================
# ===== PREPROCESS DATA =====
# ===========================


# one hot encode sex
ohe_sex = OneHotEncoder(categories = [["F", "M"]], drop="if_binary", sparse_output = False)
features_df["sex_M"] = ohe_sex.fit_transform(features_df["sex"].values.reshape(-1, 1))

# if single tracer, filter tracer
if args.single_tracer is not None:
    features_df = features_df[features_df["tracer"] == args.single_tracer]


# =============================
# ===== GET TARGET LABELS =====
# =============================


# two options:
# - use cleanest possible stable/progressor labels
#   - stables start at amyloid-positive CN and remain CN for *all* follow-up scans
#   - progressors start at amyloid-positive CN and progress to MCI or AD
#   - progressors cannot ever be classified as stable
# - recycle progressors who remain stable for n years as n-year stable

if args.alt_stable:
    data_df = get_stable_progressor_alt(
        features_df,
        max_time_to_progression = args.max_time_to_progression,
        min_time_stable = args.min_time_stable,
    )
else:
    data_df = get_stable_progressor(
        features_df,
        max_time_to_progression = args.max_time_to_progression,
        min_time_stable = args.min_time_stable,
        progressor_can_be_stable = args.progressor_can_be_stable
    )


# ============================
# ===== TRAIN/TEST SPLIT =====
# ============================


if args.leave_out_group == "site":
    hold_out_groups = features_df["site"].unique().tolist()
elif args.leave_out_group == "tracer":
    hold_out_groups = features_df["tracer"].unique().tolist()

logo_splitter = leave_one_group_out_splitter(
    data_df,
    data_df["y"].values,
    args.leave_out_group,
    hold_out_groups
)


# ===========================
# ===== DEFINE PIPELINE =====
# ===========================


# construct pipeline
feature_selector = get_feature_selector(
    args.features,
    control_standardize = args.control_standardize,
    relief = True if args.feature_reduction == "relief" else False,
    relief_n_features = param_dict["relief_inner_n"],
    pca = True if args.feature_reduction == "pca" else False,
    pca_variance_keep = param_dict["pca_variance_keep"]
)
estimator = config.model_params[args.model]["estimator"]
if args.combat:
    combat = GAMComBat(
        features = data_df.filter(regex=r".*\.amyloid").columns.values,  # get all amyloid regional SUVRs to harmonize
        covariates = ["age", "sex_M", "apoe"],
        site_name = "tracer",
    )
    pipeline = Pipeline([
        ("combat", combat),
        ("feature_selector", feature_selector),
        ("estimator", estimator)
    ])
else:
    pipeline = Pipeline([
        ("feature_selector", feature_selector),
        ("estimator", estimator)
    ])

# define ensemble model
if args.ensemble:
    model = RandomSubsetEnsembleClassifier(
        base_estimator = pipeline,
        n_estimators = args.n_estimators,
        sample_with_replacement = False,
    )
else:
    model = pipeline

# define grid search
if "param_grid" in config.model_params[args.model]:
    if args.grid_search_cv_method == "leave_one_group_out":
        cv = LeaveOneGroupOut()
    else:
        cv = StratifiedKFold(
            n_splits = param_dict["cv_splits_gridsearch"],
            shuffle = True
        )
    if args.ensemble:
        param_grid = config.model_params[args.model]["param_grid"]["ensemble"]
    else:
        param_grid = config.model_params[args.model]["param_grid"]["nonensemble"]
    model = GridSearchCV(
        model,
        param_grid = param_grid,
        scoring = "f1",
        cv = cv,
        n_jobs=-1
    )


# =====================
# ===== FIT MODEL =====
# =====================


print("+++++ Fitting models +++++")

trained_models = {}

for hold_out_group, X_train, X_test, y_train, y_test in logo_splitter:

    print(f"hold-out = {hold_out_group}")
    
    # train model
    m = clone(model)
    if isinstance(m, GridSearchCV) and args.grid_search_cv_method == "leave_one_group_out":
        m.fit(X_train, y_train, groups = X_train[args.leave_out_group].values)
    else:
        m.fit(X_train, y_train)

    # compute ROC AUC on training and test sets
    train_metrics = get_metrics(m, X_train, y_train)
    test_metrics = get_metrics(m, X_test, y_test)

    # compute Haufe feature importance
    if hasattr(m, "best_estimator_"):
        m_haufe = m.best_estimator_
    else:
        m_haufe = m
    
    if args.ensemble:
        haufe = compute_haufe_ensemble(m_haufe, X_train)
        haufe_approx = compute_haufe_ensemble(m_haufe, X_train, approx = True)
    else:
        X_transform = m_haufe[:-1].transform(X_train)
        haufe = chen_2023_haufe_weights(X_transform, m_haufe.decision_function(X_train))
        haufe_approx = chen_2023_haufe_weights(X_transform, m_haufe.predict(X_train))

    d = {
        "model": m,
        "X_train": X_train,
        "y_train": y_train,
        "X_test": X_test,
        "y_test": y_test,
        "train_metrics": train_metrics,
        "test_metrics": test_metrics,
        "haufe": haufe,
        "haufe_approx": haufe_approx,
    }
    trained_models[hold_out_group] = d

# final output dict
out_dict = {
    "trained_models": trained_models,
    "params": param_dict,
    "args": {k: v for k,v in vars(args).items()},
    "data": data_df
}


# =========================
# ===== PICKLE MODELS =====
# =========================

print("+++++ Model training complete; saving CV results +++++")

joblib.dump(
    value = out_dict,
    filename = os.path.join(args.odir, args.prefix + ".joblib")
)
