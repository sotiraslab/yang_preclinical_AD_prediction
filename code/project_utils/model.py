
import numpy as np
import torch
import torch.nn as nn
from sklearn.base import BaseEstimator, ClassifierMixin, OneToOneFeatureMixin, TransformerMixin, clone
from sklearn.compose import ColumnTransformer, make_column_selector
from sklearn.decomposition import PCA
from sklearn.pipeline import FeatureUnion, Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.utils.validation import check_is_fitted
from skrebate import ReliefF


class FullyConnectedNetwork(nn.Module):

    """
    Implementation of a vanilla fully connected network
    """

    def __init__(
        self,
        input_size,
        hidden_sizes,
        output_size
    ):

        super().__init__()

        self.input_size = input_size
        self.hidden_sizes = hidden_sizes
        self.output_size = output_size
        
        layers = []
        prev_size = input_size
        
        for h in hidden_sizes:
            layers.append(nn.Linear(prev_size, h))
            # layers.append(nn.BatchNorm1d(h))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(0.5))
            prev_size = h
        
        layers.append(nn.Linear(prev_size, output_size))
        
        self.network = nn.Sequential(*layers)
    
    def forward(self, x):

        return self.network(x)


class MultimodalNetwork(nn.Module):

    """
    Implementation of a generalized multimodal network
    
    Sub-networks may be inputted, which transform the data
    of each modality separately. The outputs of all sub-
    networks, along with additional features, is then passed
    through a final fully connected network to derive the
    output

    The size of the final output may be specified. If
    `monotonic` is set to True, the output will be restricted
    to be monotonically increasing. This is useful when you
    want the model to predict hazard rates over time (see
    Yala et al. Science 2021)

    A bias term is by default included in the final output. The
    bias is added to the output to result in the final output.
    This term is not inputted into ReLU, and therefore is
    allowed to be negative. This is useful if you want to
    allow the output to be negative but still monotonically
    increasing (for example, if you want to predict logit
    hazard rates)
    """

    def __init__(
        self,
        sub_networks,
        hidden_sizes,
        additional_features = 0,
        output_size = 1,
        monotonic = True
    ):
        
        super().__init__()
        
        self.sub_networks = nn.ModuleList(sub_networks)
        self.additional_features = additional_features
        self.hidden_sizes = hidden_sizes
        self.output_size = output_size
        self.monotonic = monotonic
        
        # Calculate the total size of all sub-network outputs plus additional features
        total_output_size = sum(net.network[-1].out_features for net in sub_networks) + additional_features
        
        # Create the combined network
        layers = []
        prev_size = total_output_size
        
        for h in hidden_sizes:
            layers.append(nn.Linear(prev_size, h))
            # layers.append(nn.BatchNorm1d(h))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(0.5))
            prev_size = h
        
        # # final linear layer
        # self.output_layer = nn.Linear(prev_size, output_size)

        # final cumulative probability layer
        if monotonic:
            # self.monotonic_relu = MonotonicReLU()
            self.cumulative_prob_layer = Cumulative_Probability_Layer(prev_size, output_size)
        else:
            self.cumulative_prob_layer = Cumulative_Probability_Layer(prev_size, output_size, make_probs_indep=True)

        # # bias term
        # if output_bias:
        #     self.bias_layer = nn.Linear(prev_size, 1)
        
        self.combined_network = nn.Sequential(*layers)
    
    def forward(self, sub_inputs, additional_input=None):
        # Process inputs through sub-networks
        sub_outputs = [net(input) for net, input in zip(self.sub_networks, sub_inputs)]
        
        # Concatenate sub-network outputs and additional features (if provded)
        if not additional_input is None:
            sub_outputs = sub_outputs + [additional_input]
        combined_input = torch.cat(sub_outputs, dim=1)
        
        # Process through combined network
        preoutput = self.combined_network(combined_input)

        # Run through final cumulative probability layer
        output = self.cumulative_prob_layer(preoutput)

        return output

        # # handle bias term
        # if self.output_bias:
        #     output_bias = self.bias_layer(preoutput)
        #     # handle monotonicity
        #     if self.monotonic:
        #         output = self.monotonic_relu(output)
            
        #     pred = output + output_bias
        #     prob = torch.sigmoid(pred)
        #     return pred, prob, output, output_bias
        # else:
        #     # handle monotonicity
        #     if self.monotonic:
        #         output = self.monotonic_relu(output)
        #     return output, torch.sigmoid(output)


# pulled from Yala et al. Science 2021
# https://github.com/yala/Mirai/blob/master/onconet/models/cumulative_probability_layer.py
class Cumulative_Probability_Layer(nn.Module):
    def __init__(self, num_features, max_followup, make_probs_indep=False):
        super(Cumulative_Probability_Layer, self).__init__()
        self.make_probs_indep = make_probs_indep
        self.hazard_fc = nn.Linear(num_features,  max_followup)
        self.base_hazard_fc = nn.Linear(num_features, 1)
        self.relu = nn.ReLU(inplace=True)
        mask = torch.ones([max_followup, max_followup])
        mask = torch.tril(mask, diagonal=0)
        mask = torch.nn.Parameter(torch.t(mask), requires_grad=False)
        self.register_parameter('upper_triagular_mask', mask)

    def hazards(self, x):
        raw_hazard = self.hazard_fc(x)
        pos_hazard = self.relu(raw_hazard)
        return pos_hazard

    def forward(self, x):
        if self.make_probs_indep:
            return self.hazards(x)
#        hazards = self.hazard_fc(x)
        hazards = self.hazards(x)
        B, T = hazards.size() #hazards is (B, T)
        expanded_hazards = hazards.unsqueeze(-1).expand(B, T, T) #expanded_hazards is (B,T, T)
        masked_hazards = expanded_hazards * self.upper_triagular_mask # masked_hazards now (B,T, T)
        cum_prob = torch.sum(masked_hazards, dim=1) + self.base_hazard_fc(x)
        return cum_prob


class MonotonicReLU(nn.Module):
    def forward(self, x):
        relu = nn.ReLU()
        x_monotonic = torch.cumsum(relu(x), dim = -1)
        return x_monotonic


class AmyloidVolumeNonimgNetwork(MultimodalNetwork):

    """
    Multimodal network using the following features:
    - regional amyloid SUVRs from PET
    - regional volume from MRI
    - non-imaging features (e.g. age, sex, APOE)
    """

    def __init__(
        self,
        amyloid_input_size,
        amyloid_hidden_sizes,
        amyloid_output_size,
        volume_input_size,
        volume_hidden_sizes,
        volume_output_size,
        additional_features,
        final_hidden_sizes,
        final_output_size = 1,
        monotonic = False
    ):

        # define subnetworks
        amyloid_subnetwork = FullyConnectedNetwork(amyloid_input_size, amyloid_hidden_sizes, amyloid_output_size)
        volume_subnetwork = FullyConnectedNetwork(volume_input_size, volume_hidden_sizes, volume_output_size)

        # init multimodal network
        super().__init__(
            [amyloid_subnetwork, volume_subnetwork],
            hidden_sizes = final_hidden_sizes,
            additional_features = additional_features,
            output_size = final_output_size,
            monotonic = monotonic,
        )

    def forward(self, x):

        """assumes that data is passed in as a tuple, where the
        first n-1 elements are passed to the sub-networks, and
        the last element is passed as additional inputs
        
        if additional inputs is 0, then just pass data into
        forward
        """

        if self.additional_features == 0:
            return super().forward(x)
        else:
            return super().forward(x[:-1], x[-1])


class ImgNonimgNetwork(MultimodalNetwork):

    """
    Single image modality network using the following features:
    - either one of (but not both):
        - regional amyloid SUVRs from PET
        - regional volume from MRI
    - non-imaging features (e.g. age, sex, APOE)
    """

    def __init__(
        self,
        img_input_size,
        img_hidden_sizes,
        img_output_size,
        additional_features,
        final_hidden_sizes,
        final_output_size = 1,
        monotonic = False
    ):

        # define subnetworks
        img_subnetwork = FullyConnectedNetwork(img_input_size, img_hidden_sizes, img_output_size)

        # init multimodal network
        super().__init__(
            [img_subnetwork],
            hidden_sizes = final_hidden_sizes,
            additional_features = additional_features,
            output_size = final_output_size,
            monotonic = monotonic,
        )

    def forward(self, x):

        if self.additional_features == 0:
            return super().forward(x)
        else:
            return super().forward(x[:-1], x[-1])


class RandomSubsetEnsembleClassifier(BaseEstimator, ClassifierMixin):

    """
    Ensemble estimator which trains a base classifier
    `n_estimators` number of times on a random subset of
    observations and aggregates the predicted probabilities
    of each classifier by taking the mean

    Skeleton of code generated by Perplexity
    """

    def __init__(
        self,
        base_estimator,
        n_estimators = 100,
        subset_frac = 0.5,
        sample_with_replacement = False,
        # subsetter,
        # group_col = None,
        random_state = None
    ):
        
        self.base_estimator = base_estimator
        # self.subsetter = subsetter
        self.n_estimators = n_estimators
        self.subset_frac = subset_frac
        self.sample_with_replacement = sample_with_replacement
        self.random_state = random_state
        # self.group_col = group_col

    def fit(self, X, y):

        if not self.random_state is None:
            np.random.seed(self.random_state)

        self.classes_, y = np.unique(y, return_inverse=True)   # https://scikit-learn.org/stable/developers/develop.html#estimator-types
        # if not self.group_col is None:
        #     groups = X.loc[:,self.group_col].values
        # else:
        #     groups = None

        self.estimators_ = []
        self.split_idx_ = []
        # n_samples = X.shape[0]
        
        for _ in range(self.n_estimators):
            
            # obtain random sample (with or without replacement)
            X_subset, y_subset, idx_subset = self.random_subsetter(
                X, y,
                subset_frac = self.subset_frac,
                sample_with_replacement = self.sample_with_replacement
            )
            
            # Clone and fit estimator
            estimator = clone(self.base_estimator)
            estimator.fit(X_subset, y_subset)
            self.estimators_.append(estimator)
            self.split_idx_.append(idx_subset)

        self.is_fitted_ = True
        return self
    
    def predict_proba(self, X):
        check_is_fitted(self, ['estimators_'])
        if hasattr(self.estimators_[0], "predict_proba"):
            probas = [estimator.predict_proba(X) for estimator in self.estimators_]
            return np.mean(probas, axis=0)
        else:
            raise AttributeError("base estimator does not have `predict_proba` method")
    
    def decision_function(self, X):
        check_is_fitted(self, ['estimators_'])
        decisions = [estimator.decision_function(X) for estimator in self.estimators_]
        return np.mean(decisions, axis=0)
    
    def predict(self, X):
        if hasattr(self.estimators_[0], "predict_proba"):
            D = self.predict_proba(X)
            return self.classes_[np.argmax(D, axis=1)]
        else:
            D = self.decision_function(X)
            return self.classes_[np.where(D > 0, 1, 0)]
        
    def random_subsetter(self, X, y, subset_frac = 0.5, sample_with_replacement = False, random_state = None):

        if not random_state is None:
            np.random.seed(random_state)

        # get random sample indices
        n = X.shape[0]
        idx = np.arange(n)
        idx_subset = np.random.choice(idx, size = np.round(n*subset_frac).astype(int), replace = sample_with_replacement)

        # subset data
        return X.iloc[idx_subset,:], y[idx_subset], idx_subset


# class StandardScalerControlGroup(StandardScaler):

#     """
#     Custom StandardScaler transformer to fit only on control group instances (y == 0).
#     This class was created so that this standard scaling operation follows scikit-learn's
#     API, which enables one to use this operation with things like cross-validation and
#     grid search

#     NOTE: this is an old version which computed the mean of the control group only.
#     However in the original Linn et al. paper, the pooled mean was used. This class
#     is now deprecated

#     Links
#     -----
#     [1] [tutorial](https://www.andrewvillazon.com/custom-scikit-learn-transformers/)
#     [2] [scikit-learn API guide](https://scikit-learn.org/stable/developers/develop.html)
#     """
    
#     def __init__(self, **kwargs):
#         super().__init__(**kwargs)
    
#     def fit(self, X, y):

#         """
#         Fit method where the base StandardScaler transformer is fit
#         only on instances from the negative class (i.e. y == 0),
#         assumed to be the control group
#         """

#         # subset is defined by the negative class of y
#         subset_idx = (y == 0)

#         # subset data matrix X
#         X = X[subset_idx]

#         # fit
#         super().fit(X)

#         return self

class StandardScalerControlGroup(OneToOneFeatureMixin, BaseEstimator, TransformerMixin):

    """
    Custom StandardScaler transformer to fit only on control group instances (y == 0).
    This class was created so that this standard scaling operation follows scikit-learn's
    API, which enables one to use this operation with things like cross-validation and
    grid search

    Standardizing by the control group has been shown to improve separability between
    groups (see Linn et al.)

    Links
    -----
    [1] [tutorial](https://www.andrewvillazon.com/custom-scikit-learn-transformers/)
    [2] [scikit-learn API guide](https://scikit-learn.org/stable/developers/develop.html)

    References
    ----------
    [1] Linn et al. Control-group feature normalization for multivariate pattern
    analysis of structural MRI data using the support vector machine
    https://doi.org/10.1016/j.neuroimage.2016.02.044
    """
    
    def __init__(self, **kwargs):
        self.mean_ = None
        self.var_ = None
        self.scale_ = None
    
    def fit(self, X, y):

        """
        Fit method where the base StandardScaler transformer is fit
        only on instances from the negative class (i.e. y == 0),
        assumed to be the control group
        """
        
        # set n_features_in_ (required to check that scaler is fitted)
        self.n_features_in_ = X.shape[1]

        # subset is defined by the negative class of y
        subset_idx = (y == 0)

        # subset data matrix X
        X_control = X[subset_idx]

        # compute mean and standard deviation
        self.mean_ = np.mean(X, axis = 0)  # pooled data
        self.var_ = np.var(X_control, axis = 0)  # control group only
        self.scale_ = np.sqrt(self.var_)

        return self

    def transform(self, X, y=None):
        return (X - self.mean_) / self.scale_

class ReliefFPipeline(Pipeline):

    """
    Custom sklearn.pipeline.Pipeline class that contains ReliefF as
    the last transformer in the pipeline. This class implements the
    `get_feature_names_out`, where it pulls the features selected by
    the ReliefF algorithm from the output features of the next-to-last
    transformer

    Links
    -----
    [1] https://stackoverflow.com/questions/48005889/get-features-from-sklearn-feature-union
    [2] https://stackoverflow.com/questions/42479370/getting-feature-names-from-within-a-featureunion-pipeline
    """

    def get_feature_names_out(self, input_features=None):

        # first check if pipeline has "top_features_" attribute
        if hasattr(self[-1], "top_features_"):

            idx = self[-1].top_features_[:self[-1].n_features_to_select]

            if not input_features is None:
                return input_features[idx]

            cols = self[:-1][0].get_feature_names_out()
            
            return cols[idx]

        # else, just run the default method (i.e. didn't use ReliefF)
        else:

            return super().get_feature_names_out()

def get_feature_selector(
    feature_subset: str,
    relief: bool = False,
    control_standardize: bool = False,
    relief_n_features: int = 20,
    pca: bool = False,
    pca_variance_keep: float = 0.9
):
    """
    Return feature selector as first step in sklearn Pipeline

    Option to include PCA to reduce dimensionality of selected features (PET and MRI features only)
    """

    # which StandardScaler to use?
    if control_standardize:
        myStandardizer = StandardScalerControlGroup
    else:
        myStandardizer = StandardScaler

    # define non-imaging feature selectors
    nonimg_cat_features = ("nonimg_cat", "passthrough", ["sex_M", "apoe"])
    nonimg_con_features = ("nonimg_cont", myStandardizer(), ["age"])
    nonimg_features = ("nonimg", ColumnTransformer([nonimg_cat_features, nonimg_con_features], verbose_feature_names_out=False))

    # define PET and MRI feature selectors
    if relief:
        # use Relief algorithm to select features
        amyloid_features = (
            "amyloid",
            ReliefFPipeline([
                ("col", ColumnTransformer([("amyloid", myStandardizer(), make_column_selector(pattern = ".*\\.amyloid"))], verbose_feature_names_out = False)),
                ("relief", ReliefF(n_features_to_select = relief_n_features, n_neighbors = 100, discrete_threshold = 10))
            ])
        )

        volume_features = (
            "volume",
            ReliefFPipeline([
                ("col", ColumnTransformer([("volume", myStandardizer(), make_column_selector(pattern = ".*\\.vol"))], verbose_feature_names_out = False)),
                ("relief", ReliefF(n_features_to_select = relief_n_features, n_neighbors = 100, discrete_threshold = 10))
            ])
        )
    elif pca:
        # learn PCA to reduce dimensionality of data
        pca_pipeline = Pipeline([
            ("scaler", myStandardizer()),
            ("pca", PCA(n_components = pca_variance_keep))
        ])
        amyloid_features = (
            "amyloid",
            ColumnTransformer([
                ("amyloid", pca_pipeline, make_column_selector(pattern = ".*\\.amyloid"))
            ])
        )
        volume_features = (
            "volume",
            ColumnTransformer([
                ("volume",  pca_pipeline, make_column_selector(pattern = ".*\\.vol"))
            ])
        )
    else:
        amyloid_features = (
            "amyloid",
            ColumnTransformer([
                ("amyloid", myStandardizer(), make_column_selector(pattern = ".*\\.amyloid"))
            ], verbose_feature_names_out=False)
        )
        volume_features = (
            "volume",
            ColumnTransformer([
                ("volume",  myStandardizer(), make_column_selector(pattern = ".*\\.vol"))
            ], verbose_feature_names_out=False)
        )

    # construct column transformer
    feature_list = []
    if "n" in feature_subset: feature_list.append(nonimg_features)
    if "a" in feature_subset: feature_list.append(amyloid_features)
    if "v" in feature_subset: feature_list.append(volume_features)
    
    return FeatureUnion(feature_list)
