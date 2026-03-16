from copy import deepcopy
import os

from imblearn.under_sampling import RandomUnderSampler
import nibabel as nib
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from plotnine import *
from sklearn.model_selection import (GridSearchCV, StratifiedGroupKFold,
                                     StratifiedKFold, StratifiedShuffleSplit,
                                     cross_val_predict, cross_val_score,
                                     cross_validate)
from torch import Tensor
from torch.optim.lr_scheduler import LRScheduler
from torch.utils.data import Subset

from .data import MultimodalDataset


class MultiwindowBase:

    """
    Base class for anything that needs to handle multiwindow
    target variables
    """

    def coalesce_multiwindow(self, target):

        # handle nan
        valid_mask = ~torch.isnan(target)
        target_no_nan = torch.where(valid_mask, target, torch.tensor(0.0))  # replace nan with 0 (won't influence any())

        # coalesce
        return target_no_nan.any(axis = 1).to(int)


class MultiwindowRepeatedStratifiedGroupKFold(MultiwindowBase):

    """
    Class for performing a (repeated) StratifiedGroupKFold
    on a pytorch Dataset where the target variable is a
    multiwindow tensor
    """

    def __init__(self, n_splits=5, n_repeats=10, random_state=None):

        self.n_splits = n_splits
        self.n_repeats = n_repeats
        self.random_state = random_state

    def split_dataset(self, dataset, target, groups = None, stratify_on_tracer = None):

        """
        Split a pytorch Dataset (specifically a MultimodalDataset)
        into train and test. The split can be optionally stratified based on
        tracer in addition to the target variable

        This function will first split the df attribute of the MultimodalDataset,
        then create new MultimodalDataset's. A new column transformer
        will be fit on the train MultimodalDataset, while the fitted column
        transfomer will be copied to the test MultimodalDataset
        """

        # coalesce multiwindow target into single targets (handle nan)
        target_coalesce = self.coalesce_multiwindow(target)

        # if tracer is provided, then stratify on tracer
        if not stratify_on_tracer is None:
            stratify_var = pd.factorize(pd.DataFrame([target_coalesce, stratify_on_tracer]).T.astype(str).agg('-'.join, axis=1))[0]  # need to concat target with tracer label, then factorize to multi-target label
        else:
            stratify_var = target_coalesce

        # begin repeated split loop
        if not self.random_state is None:
            np.random.seed(self.random_state)
        for i in range(self.n_repeats):

            # define (grouped) stratified K-fold splitter
            if not groups is None:
                splitter = StratifiedGroupKFold(
                    n_splits = self.n_splits,
                    shuffle = True
                )
            else:
                splitter = StratifiedKFold(
                    n_splits = self.n_splits,
                    shuffle = True
                )

            # yield split datasets
            for train_idx, test_idx in splitter.split(np.arange(len(dataset)), stratify_var, groups = groups):
                
                df_train, df_test = dataset.df.iloc[train_idx,:], dataset.df.iloc[test_idx,:]

                dataset_train = deepcopy(dataset)
                dataset_train.df, dataset_train.y = df_train, dataset.y[train_idx,:]
                dataset_train.fit()

                dataset_test = MultimodalDataset(
                    df_test,
                    df_test["y"].astype(int).values,
                    def_from_train = dataset_train,
                    device = dataset.device,
                    multiwindow = dataset.multiwindow,
                    control_standardize = dataset.control_standardize,
                    time_to_progression = df_test[".time_to_progression"].values,
                    time_stable = df_test[".time_stable"].values
                )

                yield dataset_train, dataset_test


class MultiwindowStratifiedSplitterForSkorch:

    """
    class with a custom method to perform a strattified shuffle split on data with a 
    multi-window target tensor. implemented to use with skorch models to perform a
    single train/validation split

    if multiwindow target == 1 at any column, then it is defined as a positive
    observation, used to stratify

    this is a solution to acheive a stratified split on both the grid search CV and
    the training of the neural net for early stopping. very hacky but whatever. also
    implemented as a class to allow for pickling of the skorch model
    """

    def __init__(self, test_size = 0.2, random_state = None, is_slicedataset = False):

        self.test_size = test_size
        self.random_state = random_state
        self.is_slicedataset = is_slicedataset

    def split(self, dataset, y, **fit_params):

        # handle SliceDataset or raw target tensor
        if self.is_slicedataset:
            # if indices is None, then SliceDataset has not been subsetted
            if y.indices is None:
                _, y_subset_tensor = y.dataset[:]
            # else use indices to slice the dataset to obtain the subset y's
            else:
                _, y_subset_tensor = y.dataset[y.indices]
        else:
            y_subset_tensor = y

        # handle nan
        valid_mask = ~torch.isnan(y_subset_tensor)
        y_subset_tensor_no_nan = torch.where(valid_mask, y_subset_tensor, torch.tensor(0.0))  # replace nan with 0 (won't influence any())

        # indicate rows where any column has positive label
        y_is_positive = y_subset_tensor_no_nan.any(axis = 1).to(int)
        
        # get stratified split
        strat_split = StratifiedShuffleSplit(n_splits=1, test_size=self.test_size, random_state=self.random_state)
        train_idx, test_idx = next(strat_split.split(np.zeros(len(dataset)), y_is_positive))

        return Subset(dataset, train_idx), Subset(dataset, test_idx)


class MultiwindowStratifiedShuffleSplitForGridSearchCV(StratifiedShuffleSplit):
    
    """
    subclass of StratifiedShuffleSplit to handle multiwindow target tensor
    
    implemented to use with GridSearchCV
    """

    def __init__(self, **args):

        super().__init__(**args)

    def split(self, X, y=None, groups=None):

        """
        Assumes that the input X and y are skorch SliceDataset objects

        y's dataset should have the attribute .y, which contains the
        ground truth multi-window labels for each row of data
        """
        
        y_is_positive = y.dataset.y.any(axis = 1).to(int)  # collapse y
        return super().split(X, y_is_positive, groups)


class MaskedMultiwindowBCEWithLogitsLoss(nn.Module):

    """
    Custom loss function that computes nan-masked binary cross-entropy loss
    on a multi-time window target tensor
    
    Skeleton of code generated by Perplexity
    """

    def __init__(self, weight=None, reduction='mean', pos_weight=None):

        super().__init__()
        self.weight = weight
        self.reduction = reduction
        self.pos_weight = pos_weight

    def forward(self, logits, targets):

        n_samples = targets.shape[0]

        # Create a mask for valid (non-NaN) targets
        valid_mask = ~torch.isnan(targets)

        # Filter out invalid logits and targets using the mask
        logits = logits[valid_mask]
        targets = targets[valid_mask]

        # Hangle weight
        if not self.weight is None:
            weight_broadcast = self.weight.unsqueeze(0).repeat(n_samples, 1)
            weight = weight_broadcast[valid_mask]
        else:
            weight = None

        # Handle positive weight
        if not self.pos_weight is None:
            pos_weight_broadcast = self.pos_weight.unsqueeze(0).repeat(n_samples, 1)
            pos_weight = pos_weight_broadcast[valid_mask]
        else:
            pos_weight = None

        # Compute BCEWithLogitsLoss on valid data
        return F.binary_cross_entropy_with_logits(
            logits,
            targets,
            weight=weight,
            reduction=self.reduction,
            pos_weight=pos_weight,
        )


def grouped_random_undersample(X, y, groups, random_state = None):

    """
    Custom splitter function which returns a random subsample of
    data that organizes into groups (e.g., subject). A single random
    instance per group is sampled. Additionally, the majority class
    is undersampled to match the count of the minority class (uses
    imblearn RandomUnderSampler)
    """

    if not random_state is None:
        np.random.seed(random_state)

    # Get unique groups
    unique_groups = np.unique(groups)
    
    # Sample one index per group
    sampled_indices = []
    for group in unique_groups:
        group_indices = np.where(groups == group)[0]  # Get indices of the current group
        sampled_index = np.random.choice(group_indices)  # Randomly select one index
        sampled_indices.append(sampled_index)
    sampled_indices = np.array(sampled_indices)

    # perform random undersampling
    rus = RandomUnderSampler()
    idx_rus, _ = rus.fit_resample(sampled_indices.reshape(-1, 1), y[sampled_indices])
    idx_rus = idx_rus.flatten()

    # return samples
    return X.iloc[idx_rus,:], y[idx_rus], idx_rus


def leave_one_group_out_splitter(X, y, group_col, hold_out_groups):

    for g in hold_out_groups:
        if not isinstance(g,list):
            g = [g]
        idx = X[group_col].isin(g)
        X_train, X_test = X[~idx], X[idx]
        y_train, y_test = y[~idx], y[idx]
        s = "_".join(g)

        yield s, X_train, X_test, y_train, y_test


class CustomReduceLROnPlateau(LRScheduler):

    """Custom learning rate scheduler. This scheduler will reduce
    the learning rate of an optimizer when a validation score fails
    to increase for `patience` number of epochs in a row. It will
    signal to stop training once the validation score fails to
    increase for another `patience` number of epochs.
    
    Skeleton of code generated by Perplexity
    """
    
    def __init__(self, optimizer, factor=0.2, repeats=1, patience=5, maximize_loss=False, verbose=False):
        self.optimizer = optimizer
        self.factor = factor
        if repeats < 0:
            raise ValueError("number of repeats must be a non-negative integer")
        else:
            self.repeats = repeats
        self.patience = patience
        self.maximize_loss = False
        self.verbose = verbose
        
        # internal attributes
        self.stop_training = False
        self._best = None
        self._num_bad_epochs = 0
        self._is_better = lambda a, b: a > b if maximize_loss else a < b
        self._last_lr = [group['lr'] for group in self.optimizer.param_groups]
        self._repeats_left = repeats

    def step(self, metrics):
        current = metrics
        if self._best is None:
            self._best = current
            return

        if self._is_better(current, self._best):
            self._best = current
            self._reduced_lr = False  # reset so that LR can be reduced further
            self._num_bad_epochs = 0
        else:
            self._num_bad_epochs += 1

        if self._num_bad_epochs > self.patience:

            if self._repeats_left <= 0:
                self.stop_training = True
                if self.verbose:
                    print('Stopping training due to lack of improvement.')
            else:
                self._reduce_lr(self._num_bad_epochs)
                self._num_bad_epochs = 0
                self._repeats_left -= 1

    def _reduce_lr(self, epoch):
        for i, param_group in enumerate(self.optimizer.param_groups):
            old_lr = float(param_group['lr'])
            new_lr = old_lr * self.factor
            param_group['lr'] = new_lr
            self._last_lr[i] = new_lr
            if self.verbose:
                print(f'No improvement for {self.patience} epochs; reducing learning rate of group {i} to {new_lr:.4e}.')

    def get_last_lr(self):
        return self._last_lr


class LossTracker:

    """Custom loss function tracker class

    Skeleton of code generated from Perplexity
    """

    def __init__(self, loss_list):
        self.losses = {key: [] for key in loss_list}
        self.losses_sum = {key: 0.0 for key in loss_list}
        self.total_loss = 0.0
        
        self.losses_history = []
        self.losses_sum_history = []

    def update(self, **kwargs):

        """
        Append a loss to a list of batch losses; run this after computing
        the loss for a single batch
        """
        
        for key, value in kwargs.items():
            if isinstance(value, Tensor):
                value = value.item()
            self.losses[key].append(value)
            self.losses_sum[key] += value
            self.total_loss += value

    def reset(self):

        """
        Append the cumulative sum of loss to the history and reset
        to zero; run this after completion of an epoch
        """
        
        self._append_history()
        self.losses = {key: [] for key in self.losses.keys()}
        self.losses_sum = {key: 0.0 for key in self.losses.keys()}
        self.total_loss = 0.0

    def _append_history(self):

        self.losses_history.append(self.losses)
        
        d = self.losses_sum
        d["total"] = self.total_loss
        self.losses_sum_history.append(d)

    def get_average_losses(self, num_batches):
        return {key: value / num_batches for key, value in self.losses_sum.items()}

    def get_average_total_loss(self, num_batches):
        return self.total_loss / num_batches

    def get_loss_history(self, long = False):
        df = pd.DataFrame(self.losses_sum_history).reset_index().rename(columns = {"index": "epoch"})
        if long:
            df = df.melt(id_vars = "epoch", var_name = "loss type", value_name = "loss value")
        return df

    def plot_loss(self):

        df = self.get_loss_history(long = True)
        g = ggplot(df, aes(x = "epoch", y = "loss value")) \
            + geom_line()
        return g

    def __str__(self):
        return " | ".join([f"{key}: {value:.4f}" for key, value in self.losses_sum.items()])


def nested_cv(
    model,
    X, y,
    param_grid,
    inner_cv,
    outer_cv,
    outer_groups = None,
    inner_scoring = "roc_auc",
    outer_scoring = "roc_auc",
    return_predict = False,
    return_cv = False
):

    """
    Perform nested cross-validation and report estimator scoring for each outer
    loop testing evaluation

    Parameters
    ----------
    model : Any
        sklearn estimator
    X, y : np.ndarray, shape (m,n) and (m,1)
        array of features and target labels
    param_grid : dict
        dictionary of parameters to grid search over in inner CV
    inner_cv : BaseCrossValidator
        sklearn cross valdiator instance for inner CV (hyperparameter tuning)
    outer_cv : BaseCrossValidator
        sklearn cross validator instance for outer CV (generalizability eval)
    outer_groups : np.ndarray, shape (m,1) (default = None)
        array of groups to perform grouped CV on outer CV such as LeaveOneGroupOut;
        ignored if outer_cv is not a grouped CV method
    inner_scoring : str (default = "roc_auc)
        scoring method for inner CV, which determines the optimal hyperparameters of
        each inner estimator
    outer_scoring : str or list of str (default = "roc_auc")
        scoring method(s) for outer CV, which computes the generalizability of the model
        to the outer folds; only used when return_predict = False
    return_predict : bool (default = False)
        whether to use `cross_val_predict` and return predictions rather than score
    return_cv : bool (default = False)
        whether to use `cross_validation` and return full CV results
    
    Returns
    -------
    np.ndarray
        array of outer CV testing scores (or predictions if return_predict = True)
    -- or --
    dict (return_cv = True)
        CV dictionary

    Links
    -----
    [1] [nested CV](https://scikit-learn.org/stable/auto_examples/model_selection/plot_nested_cross_validation_iris.html)
    [2] [nested CV](https://inria.github.io/scikit-learn-mooc/python_scripts/cross_validation_nested.html)
    [3] [access group in LeaveOneGroupOut](https://stackoverflow.com/questions/71759832/how-to-access-the-group-value-for-each-cv-iteration-in-sklearn-leaveonegroupout)
    """
    
    # inner CV estimator is grid search
    inner_est = GridSearchCV(
        estimator = model,
        param_grid = param_grid,
        cv = inner_cv,
        scoring = inner_scoring
    )
    
    # run nested CV
    if return_predict:
        return cross_val_predict(
            estimator = inner_est,
            X = X, y = y,
            cv = outer_cv,
            groups = outer_groups
        )

    if return_cv:
        return cross_validate(
            estimator = inner_est,
            X = X, y = y,
            cv = outer_cv,
            groups = outer_groups,
            scoring = outer_scoring,
            return_train_score = True,
            return_estimator = True,
            return_indices = True
        )

    return cross_val_score(
        estimator = inner_est,
        X = X, y = y,
        cv = outer_cv,
        groups = outer_groups,
        scoring = outer_scoring
    )

def get_relief_features(relief_est, colnames):

    """
    Extract feature names selected by ReliefF

    Parameters
    ----------
    relief_est : skrebate.ReliefF estimator
        fitted ReliefF estimator
    colnames : np.ndarray
        array of column names; should be the same length as the
        number of columns inputted into ReliefF during fitting

    Returns
    -------
    np.ndarray
        names of columns selected by ReliefF

    Notes
    -----
    - this function will select the same number of column names
      that was selected by ReliefF (`n_features_to_select` param)
    
    """

    idx = relief_est.top_features_[:relief_est.n_features_to_select]
    return colnames[idx]

def get_relief_features_from_pipeline(pipeline):

    """
    Wrapper function for `get_relief_features` to operate on a
    sklearn.Pipeline, where the first transformer is a
    ColumnTransformer and the second transformer is a ReliefF
    """

    return get_relief_features(pipeline["relief"], pipeline["col"].get_feature_names_out())

class NiftiChecker:

    """
    Class that contains helper functions to check a list of nifti images
    """

    def __init__(self, nifti_list = None):

        self.nifti_list = nifti_list

    def find_nifti_files(self, input_dir):
        """
        Find all NIfTI files in the given directory and its subdirectories.
        Set the list of NIfTI files as an attribute.
        """
        nifti_files = []
        for root, _, files in os.walk(self.input_dir):
            for file in files:
                if file.endswith(('.nii', '.nii.gz')):
                    nifti_files.append(os.path.join(root, file))
        
        self.nifti_list = nifti_files
        self.input_dir = input_dir

    def get_shape(self):
        """
        Get the shape of each NIfTI file.
        Returns a dictionary with file paths as keys and shape tuples as values.
        """
        return {file: self.get_shape_single(file) for file in self.nifti_list}

    def count_frames(self):
        """
        Count the number of frames (third dimension or time points) for each NIfTI file.
        Returns a dictionary with file paths as keys and frame counts as values.
        """
        return {file: self.count_frames_single(file) for file in self.nifti_list}

    @staticmethod
    def get_shape_single(file):
        try:
            img = nib.load(file)
            return img.shape
        except FileNotFoundError:
            return None

    @staticmethod
    def count_frames_single(file):
        try:
            img = nib.load(file)
            return img.shape[-1] if len(img.shape) > 3 else 1
        except FileNotFoundError:
            return None