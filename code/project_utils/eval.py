
from dataclasses import dataclass

import numpy as np
import pandas as pd
from plotnine import *
from sklearn.base import clone
from sklearn.covariance import EmpiricalCovariance
from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve
from sklearn.pipeline import Pipeline

# ====================================
# ===== HAUFE FEATURE IMPORTANCE =====
# ====================================

class HaufeTransform():

    def haufe_correction(self, X_train, coef = None, model = None, covariance_estimator = EmpiricalCovariance):

        """
        Rotate linear discriminator weights using covariance matrix of
        the data; see [1]

        Parameters
        ----------
        X_train : numpy.array or pandas.Series, shape (m,n)
            array of training data with m observations and n features; if
            pandas.Series, preserves the indices
        coef : numpy.ndarray, shape (n,k)
            array containing true model weights of n features (rows); multiple
            coefficient vectors may also be specified (k columns)
        model : sklearn estimator
            sklearn estimator with `coef_` attribute; if specified, `coef_`
            attribute will be used in place of the "coef" argument
        covariance_estimator : sklearn class
            sklearn instance of covariance estimator; defaults to empirical
            (maximum likelihood) estimator

        Returns
        -------
        numpy.ndarray, shape (n,k)
            weights rotated by covariance of X_train

        Reference
        ---------
        [1] S. Haufe, "On the interpretation of weight
            vectors of linear models in multivariate neuroimaging,"
            NeuroImage, vol. 87, pp. 96-110, 2014.
        """

        # TODO: write unit tests for this function
        # TODO: add axis argument to determine which dimension defines the different coef vectors

        # estimate covariance
        cov = covariance_estimator().fit(X_train).covariance_

        # get weights
        if not model is None:
            if hasattr(model, "coef_"):
                coef = model.coef_
            else:
                coef = model[-1].coef_
        elif coef is None:
            raise ValueError("one of `coef` or `model` must be specified")

        # handle pandas
        columns = None
        index = None
        if isinstance(coef, pd.Series):
            index = coef.index
            val = coef.values
        elif isinstance(coef, pd.DataFrame):
            index = coef.index
            columns = coef.columns
            val = coef.values
        else:
            val = coef

        # if 1D array or 2D row vector, reshape to 2D column vector
        if val.ndim == 1 or val.shape[0] == 1:
            val = val.reshape(-1, 1)
        
        # rotate weights
        coef_haufe = cov @ val

        # reformat to pandas, or reshape
        if not columns is None:
            coef_haufe = pd.DataFrame(coef_haufe, index = index, columns = columns)
        elif not index is None:
            coef_haufe = pd.Series(coef_haufe.flatten(), index = index)
        else:
            coef_haufe = coef_haufe.reshape(coef.shape)

        return coef_haufe

    def get_haufe_from_pipeline(self, pipeline, X):

        """
        Return covariance-rotated linear discriminator weights directly
        from sklearn pipeline

        Parameters
        ----------
        pipeline : sklearn.pipeline.Pipeline
            scikit-learn Pipeline where the last step is the estimator
            with `coef_` attribute
        X : numpy.ndarray
            input data to compute covariance matrix; this will be run
            through the pipeline (excluding last step) to transform
            before computing the covariance matrix

        Returns
        -------
        pd.Series
            Series containing covariance-rotated weights; if feature
            names are available, they are set to the index
        """

        # transform X
        X_transform = pipeline[:-1].transform(X)

        # get Haufe coefficients
        coef = haufe_correction(
            X_transform,
            coef = pipeline[-1].coef_
        )

        # return Series
        if hasattr(pipeline[:-1], "get_feature_names_out"):
            feature_names = pipeline[:-1].get_feature_names_out()
            coef = pd.Series(coef.flatten(), index = feature_names)
        else:
            coef = pd.Series(coef.flatten())

        return coef

        coef = gs_est[-1].coef_

    def permutation_test_coef(self, model, X_train, y_train, feature_names, num_permute=1000, is_pipeline=False, random_state=None):

        """
        Perform permutation testing to determine feature importance

        Note: if the model is a sklearn pipeline, then the features
        are transformed first by applying all steps of the pipeline
        except for the last, then is fed into the permutation testing
        loop

        Parameters
        ----------
        model :
            scikit-learn estimator with `coef_` attribute
        X_train : numpy.array, shape (m,n)
            array of training data with m observations and n features
        y_train : numpy.array, shape (m,)
            vector of labels associated with training data
        feature_names : numpy.array or list-like, shape (n,)
            names of features associated with data in X_train
        num_permute : int (default = 1000)
            (optional) number of permutations to perform
        is_pipeline : bool (default = False)
            indicates whether argument `model` is a scikit-learn
            pipeline object; if so, then X_train is transformed by
            all steps of the pipeline except for the last
        random_state : int
            (optional) random state of numpy permutation function

        Returns
        -------
        pandas.DataFrame, shape (num_permute, n)
            dataframe containing random model weights, where rows are
            permutations and columns are weights associated with each
            input feature
        pandas.Series, shape (n,)
            series containing true model weights
        """

        # create numpy RandomState
        # https://stackoverflow.com/questions/22994423/difference-between-np-random-seed-and-np-random-randomstate
        r = np.random.RandomState(random_state)

        # if pipeline, transform X_train
        if is_pipeline:
            X_train = model[:-1].transform(X_train)
            estimator = model[-1]
        else:
            estimator = model

        # get true coefficients
        true_coef = estimator.coef_.flatten()
        n_features = true_coef.size

        # permute y-labels, fit random models, and get null distribution of model coef
        coef_array = np.zeros((n_features, num_permute))
        for i in range(num_permute):
            # randomize labels
            y_train_shuff = r.permutation(y_train)
            
            # initialize new random model from inputted model
            estimator_rand = clone(estimator)

            # fit model on random data
            estimator_rand.fit(X_train, y_train_shuff)

            # get coefficients and store in dictionary
            coef_array[:, i] = estimator_rand.coef_.flatten()

        # reformat to dataframe, series
        coef_df = pd.DataFrame(coef_array.T, columns = feature_names)
        true_coef = pd.Series(true_coef, index = feature_names)

        return coef_df, true_coef

    def permutation_test_pval(self, rand_coef, true_coef, alpha = 0.05, multiple_comparison = None):

        """
        Compute p-values of true model weights (either corrected or
        uncorrected) compared to null distribution

        Parameters
        ----------
        rand_coef : pandas.DataFrame
            dataframe containing random model weights, where rows are
            permutations and columns are weights associated with each
            input feature
        true_coef : pandas.Series
            series containing true model weights
        alpha : float (default = 0.05)
            alpha value to determine significance; note that since the test
            is 2-sided, this will be divided by 2
        multiple_comparison : str
            (optional) indicate whether to perform multiple comparison
            correction; select a string from the `method` argument of
            the function statsmodels.stats.multitest.multipletests:
            https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html

        Returns
        -------
        numpy.array
            p-values
        numpy.array
            whether a feature is statistically significant, depending
            on given alpha level
        """

        # get number of instances where permuted coef is less than / greater than true coef
        lt_count = rand_coef.lt(true_coef).sum(axis=0) # compute count in null dist less than actual coef
        gt_count = rand_coef.gt(true_coef).sum(axis=0) # compute count in null dist greater than actual coef

        # compute 2-sided p-value
        pvals = np.minimum(lt_count, gt_count) / rand_coef.shape[0] # for 2-sided test, p-value is the minimum of the two divided by number of iterations

        # correct for multiple comparisons
        if not multiple_comparison is None:
            reject, pvals, _, _ = multipletests(
                pvals = pvals,
                alpha = alpha / 2,
                method = multiple_comparison
            )
            reject = pd.Series(reject, index=true_coef.index)
            pvals  = pd.Series(pvals,  index=true_coef.index)
        else:
            reject = pvals < (alpha / 2)
        
        return pvals, reject

    def permutation_test_wrapper(self, pipeline, X_train, y_train, num_permute = 1000, random_state = None):
        
        # wrapper function to generate random coefficients,
        # do Haufe correction and compute p-values

        # get feature names and modify
        pipeline_feature_names = pipeline[:-1].get_feature_names_out()

        # compute random coefficients
        rand_coef, true_coef = permutation_test_coef(
            model = pipeline,
            X_train = X_train,
            y_train = y_train,
            feature_names = pipeline_feature_names,
            num_permute = num_permute,
            random_state = random_state,
            is_pipeline = True
        )

        # apply Haufe linear weight correction
        # NOTE: since covariance is computed on transformed X matrix,
        #   we must first transform it before feeding it to function
        X_transform = pipeline[:-1].transform(X_train)
        rand_coef_haufe = haufe_correction(X_transform, coef = rand_coef.T).T
        true_coef_haufe = haufe_correction(X_transform, coef = true_coef)

        # compute p-values
        pval, _ = permutation_test_pval(rand_coef, true_coef)
        pval_haufe, _ = permutation_test_pval(rand_coef_haufe, true_coef_haufe)

        return pval, pval_haufe, true_coef, true_coef_haufe


@dataclass
class BaseEvaluator:

    """
    Base class for evaluating the performance of trained pipelines
    and visualizing their results

    Uses trained sklearn Pipelines and data to 
    """

    pipeline: Pipeline = None
    X: pd.DataFrame = None
    y: np.ndarray | pd.DataFrame = None

    def _default_X(self, X = None):

        if X is None:
            return self.X
        else:
            return X

    def _default_y(self, y = None):

        if y is None:
            return self.y
        else:
            return y

    def get_feature_names(self):

        return self.pipeline[-1].get_feature_names_out()

    def get_coef(self):

        return self.pipeline[-1].coef_
    
    def set_Xy(self, X, y):

        self.X = X
        self.y = y

    def get_prediction(self, X = None):

        """requires that pipeline's estimator has `predict` method"""

        if self.pipeline is None:
            raise ValueError("A scikit-learn pipeline must be inputed to the evaluator!")

        return self.pipeline.predict(self._default_X(X))

    def get_prediction_class(self, X = None, y = None, string_class = False):

        """
        Compute the classification results (true positive, true negative,
        false positive, false negative) for a set of true and predicted labels.
    
        Args:
            y_true (list or numpy.ndarray): The true labels.
            y_pred (list or numpy.ndarray): The predicted labels.
            
        Returns:
            numpy.ndarray: A vector of length len(y_true) with values 0, 1, 2, or 3, indicating:
                0 = true negative
                1 = true positive
                2 = false positive
                3 = false negative
        
        Notes:
            requires that pipeline's estimator has `predict` method
        """

        y_true = self._default_y(y)
        y_pred = self.get_prediction(X)

        # Compute the classification results
        true_positive = (y_true == 1) & (y_pred == 1)
        true_negative = (y_true == 0) & (y_pred == 0)
        false_positive = (y_true == 0) & (y_pred == 1)
        false_negative = (y_true == 1) & (y_pred == 0)
        
        # Combine the results into a single vector
        if string_class:
            results = np.zeros_like(y_true, dtype="U2")
            results[true_positive] = "TP"
            results[true_negative] = "TN"
            results[false_positive] = "FP"
            results[false_negative] = "FN"
        else:
            results = np.zeros_like(y_true, dtype=int)
            results[true_positive] = 1
            results[true_negative] = 0
            results[false_positive] = 2
            results[false_negative] = 3

        return results
    
    def get_scores_binary(self, X=None, y_true=None, y_pred=None, conf_mat=None):

        """
        Compute accuracy metrics for binary classification problem.
        Metrics include: accuracy, sensitivity (recall), specificity,
        ppv (precision), npv, balanced_accuracy, f1

        User may input either a pair of y_true and y_pred vectors, or
        a precomputed confusion matrix

        Parameters
        ----------
        X : pandas.DataFrame
            dataframe of data to input into pipeline
        y_true : numpy.array of int or float
            array of the true target labels, where 0 is the negative
            class and 1 is the positive class
        y_pred : numpy.array of int or float
            array of the predicted target labels
        conf_mat : numpy.array, shape (2,2) (optional)
            precomputed confusion matrix; if none specified, the
            confusion matrix will be derived from y_true and y_pred;
            y_true and y_pred will be ignored if conf_mat is
            specified

        Returns
        -------
        dict
            dictionary containing accuracy metrics
        """

        # initialize dictionary to store all metrics
        metric_dict = {}

        # compute confusion matrix
        if conf_mat is None:
            if y_pred is None:
                y_pred = self.get_prediction(self._default_X(X))
            conf_mat = confusion_matrix(y_true=self._default_y(y_true), y_pred=y_pred)
        
        # get TP, TN, FP, FN
        TN = conf_mat[0,0]; FN = conf_mat[1,0]; FP = conf_mat[0,1]; TP = conf_mat[1,1]

        # accuracy
        metric_dict["accuracy"] = (TN + TP) / (TN + FN + FP + TP)
        # sensitivity
        metric_dict["sensitivity"] = (TP) / (TP + FN)
        # specificty
        metric_dict["specificity"] = (TN) / (TN + FP)
        # PPV
        metric_dict["ppv"] = (TP) / (TP + FP)
        # NPV
        metric_dict["npv"] = (TN) / (TN + FN)
        # precision (=ppv)
        metric_dict["precision"] = metric_dict["ppv"]
        # recall (=sensitivity)
        metric_dict["recall"] = metric_dict["sensitivity"]
        # balanced accuracy
        metric_dict["balanced_accuracy"] = (metric_dict["sensitivity"] + metric_dict["specificity"]) / 2
        # F1 score
        metric_dict["f1"] = 2 * (metric_dict["precision"]*metric_dict["recall"]) / (metric_dict["precision"]+metric_dict["recall"])

        return metric_dict


@dataclass
class LinearWeightEvaluator(BaseEvaluator, HaufeTransform):

    """
    Evaluator class for plotting a set of linear discriminator weights.
    Also provides Haufe covariance-rotation functionality
    """

    haufe: bool = False

    def __post_init__(self):

        # get coefficient
        try:
            self.coef = self.pipeline[-1].coef_.flatten()
        except:
            raise ValueError("inputed pipeline's estimator (last step) does not have `coef_` attribute!")
        
        # Haufe transform
        if self.haufe:

            # X is required for this
            if self.X is None:
                raise ValueError("X must be inputed for Haufe correction to work")

            self.coef = self.haufe_correction(
                X_train = self.pipeline[:-1].transform(self.X),
                coef = self.coef
            ) 
        
        # get feature names
        try:

            # get feature names of pipeline
            self.feature_names = self.pipeline[:-1].get_feature_names_out()

        # else assign arbitrary feature names
        except AttributeError:

            n = self.coef.size
            self.feature_names = np.array([f"feature_{str(i+1).zfill(int(np.floor(np.log10(n)) + 1))}" for i in range(n)])

        # get coef as DataFrame
        self.coef_df = pd.Series(self.coef, index = self.feature_names) \
            .rename("coef") \
            .reset_index() \
            .sort_values(by = "coef") \
            .reset_index()
        self.coef_df[["modality", "feature_name"]] = self.coef_df["index"] \
            .str.split("__", expand = True)
        self.coef_df["coef_abs"] = abs(self.coef_df["coef"])


    def plot_col(self, n_vis = None, sort_by_abs = False, ascending = False, flip_xy = False, df = None):

        # get data
        if not df is None:
            df = df.copy()
        else:
            df = self.coef_df.copy()

        # sort by absolute value
        y_sort = "coef_abs" if sort_by_abs else "coef"

        # ascending or descending
        df["y_sort"] = df[y_sort] if ascending else -df[y_sort]
        if ascending:
            x_aes = f"reorder(feature_name, {y_sort})"
        else:
            x_aes = f"reorder(feature_name, -{y_sort})"

        # number of features to visualize
        df = df.nsmallest(n_vis, "y_sort") if not n_vis is None else df

        # plot
        g = ggplot(
            df,
            aes(x = x_aes, y = "coef", fill = "modality")
        ) \
            + geom_col() \
            + theme(
                axis_text_x = element_text(rotation = 90)
            ) \
            + labs(x = "feature name", y = "weight")

        if flip_xy: g = g + coord_flip()
        
        return g

    def plot_col_pos_neg(self, n_vis = None, sort_by_abs = False, ascending = False, flip_xy = False):

        # split positive and negative
        df_pos = self.coef_df[self.coef_df["coef"] >= 0].copy()
        df_pos["direction"] = "positive"
        df_neg = self.coef_df[self.coef_df["coef"] < 0].copy()
        df_neg["direction"] = "negative"
        
        # n most important
        y_sort = "coef_abs" if sort_by_abs else "coef"
        if not n_vis is None:
            df_pos = df_pos.nlargest(n_vis, y_sort)
            df_neg = df_neg.nsmallest(n_vis, y_sort)

        # concatenate
        df = pd.concat([df_pos, df_neg], ignore_index = True)

        # plot
        g = self.plot_col(
            df = df,
            sort_by_abs = sort_by_abs,
            ascending = ascending,
            flip_xy = flip_xy
        )

        if flip_xy:
            g = g + facet_wrap("direction", ncol = 1, scales = "free_y")
        else:
            g = g + facet_wrap("direction", nrow = 1, scales = "free_x")

        return g


@dataclass
class ConfusionEvaluator(BaseEvaluator):

    """
    Evaluator class for confusion matrix.

    The input pipeline's estimator (last step) is required to
    have a `predict` method.
    """

    @staticmethod
    def confusion_matrix_to_long(mat, row_labels, col_labels):

        """
        Convert an n-by-n confusion matrix to a long DataFrame.
        
        Parameters
        ----------
        mat : list of lists
            input matrix
        row_labels : list
            labels for the rows
        col_labels : list
            labels for the columns.
        
        Returns
        -------
        pd.DataFrame
            long DataFrame with columns for row label, column label, and value.
        """

        data = []
        for i, row in enumerate(mat):
            for j, count in enumerate(row):
                data.append({
                    'true label': row_labels[i],
                    'predicted label': col_labels[j],
                    'count': count
                })
        
        return pd.DataFrame(data)
    
    def get_confusion_matrix(
        self,
        X = None,
        y_true = None,
        y_pred = None,
        long = False,
        long_labels = ("stable", "progressor")
    ):

        """
        compute confusion matrix

        Parameters
        ----------
        X :
            input data; should be in same format as inputted into self.estimator
        y : 
            ground truth labels
        long : bool
            if True, return confusion matrix formatted as long DataFrame
        long_labels : list of str
            labels for confusion matrix rows/columns
        
        Returns
        -------
        numpy.ndarray, or pd.DataFrame
            confusion matrix; if long = True, returns in long format
        """

        if y_pred is None:
            y_pred = self.get_prediction(self._default_X(X))

        cm = confusion_matrix(
            y_true = self._default_y(y_true),
            y_pred = y_pred
        )

        if long:
            return self.confusion_matrix_to_long(cm, long_labels, long_labels)
        
        return cm

    def plot_confusion_matrix_from_long(self, df):

        """
        plot confusion matrix from a long formatted df; useful
        if multiple matrices need to be plotted in facets
        """

        g = ggplot(df, aes(x = "predicted label", y = "true label", fill = "count")) \
            + geom_tile() \
            + geom_text(aes(label = "count"), color = "red") \
            + coord_trans(x='reverse') \
            + theme(axis_text_x = element_text(angle = 45, hjust=1))
        return g

    def plot_confusion_matrix(self, X, y):

        cm_long = self.get_confusion_matrix(
            self._default_X(X),
            self._default_y(y),
            long = True
        )
        return self.plot_confusion_matrix_from_long(cm_long)


@dataclass
class ROCEvaluator(BaseEvaluator):

    def get_prediction_prob(self, X=None):

        if self.pipeline is None:
            raise ValueError("A scikit-learn pipeline must be inputed to the evaluator!")

        # compute classifier score
        if hasattr(self.pipeline, "decision_function"):
            y_score = self.pipeline.decision_function(self._default_X(X))
        else:
            y_score = self.pipeline.predict_proba(self._default_X(X))[:,1]
        
        return y_score
    
    def compute_roc(self, X=None, y_true=None, y_score=None, return_df = False):

        # get score if none provided
        if y_score is None:
            y_score = self.get_prediction_prob(X=self._default_X(X))

        # compute ROC curve coordinates
        fpr, tpr, _ = roc_curve(y_true=self._default_y(y_true), y_score=y_score)

        # compute ROC AUC
        auc = roc_auc_score(y_true=self._default_y(y_true), y_score=y_score)

        if return_df:
            df = pd.DataFrame({"fpr": fpr, "tpr": tpr})
            s_auc = pd.Series({"auc": auc})
            return df, s_auc
        else:
            return fpr, tpr, auc

    def plot_roc(self, X=None, y_true=None, y_score=None):

        df, auc = self.compute_roc(X=X, y_true=y_true, y_score=y_score, return_df=True)
        return self.plot_roc_from_df(df, title = f'ROC curve (AUC={np.round(auc,3)})')

    def plot_roc_from_df(self, df, color=None, title='ROC curve', auc_df = None):

        # TODO: plot from a dictionary of dataframes, concatenate and use keys as new column
        # TODO: add facetting options
        # TODO: add custom ROC legend

        if color is None:
            roc_aes = aes(x='fpr', y='tpr')
        else:
            roc_aes = aes(x='fpr', y='tpr', color=color, group=color)

        g = (ggplot(df, roc_aes)
            + geom_line(size = 1.5, alpha = 0.66)
            + geom_abline(intercept=0, slope=1, linetype='dashed', color='gray')
            + labs(title=title,
                    x='false positive rate',
                    y='true positive rate')
            + xlim(0,1) + ylim(0,1)
            + coord_fixed())
        
        if not auc_df is None:
            auc_df["auc_round"] = "AUC=" + auc_df["auc"].astype(float).round(3).astype(str)
            g = (g + geom_text(
                         data = auc_df,
                         mapping = aes(label = "auc_round"),
                         x = 0.75, y = 0.05,
                     ))

        return g
