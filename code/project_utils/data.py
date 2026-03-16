
import numpy as np
import pandas as pd
import torch
from sklearn.compose import ColumnTransformer, make_column_selector
from sklearn.preprocessing import StandardScaler
from torch.utils.data import Dataset

from .model import StandardScalerControlGroup


class MultimodalDataset(Dataset):

    """
    Given a pandas DataFrame with PET, MRI and non-imaging features, construct a
    pytorch dataset

    Features
    - columns ending with ".amyloid" are grouped as the PET features
    - columns ending with ".volume" are grouped as the MRI features
    - columns belonging in the list of non-imaging features are grouped as the non-imaging features
    """
    
    def __init__(
        self,
        df,
        y,
        multiwindow = False,
        time_to_progression = None,
        time_stable = None,
        num_time_windows = 5,
        features_to_extract = "avn",
        nonimg_cat_features = None,
        nonimg_con_features = None,
        control_standardize = False,
        combat_model = None,
        def_from_train = None,
        dtype = torch.float32,
        device = torch.device("cpu")
    ):
        
        self.df = df
        self.y = torch.tensor(y, device = device)
        self.dtype = dtype
        self.device = device

        # multiwindow
        if multiwindow:
            self.time_to_progression = time_to_progression
            self.time_stable = time_stable
            self.multiwindow = multiwindow
            self.num_time_windows = num_time_windows

        # set all attributes from a training dataset
        if not def_from_train is None:

            self.features_to_extract = features_to_extract
            self.control_standardize = def_from_train.control_standardize
            self.combat_model = def_from_train.combat_model
            
            self.nonimg_cat_features = def_from_train.nonimg_cat_features
            self.nonimg_con_features = def_from_train.nonimg_con_features
            self.nonimg_transform = def_from_train.nonimg_transform
            self.amyloid_transform = def_from_train.amyloid_transform
            self.volume_transform = def_from_train.volume_transform

            # update features
            self.features = self.transform()
        
        # else define transformer
        else:

            self.features_to_extract = features_to_extract
            self.control_standardize = control_standardize
            self.combat_model = combat_model
            
            # which StandardScaler to use?
            if control_standardize:
                myStandardizer = StandardScalerControlGroup
            else:
                myStandardizer = StandardScaler

            # define non-imaging feature selectors
            if "n" in self.features_to_extract:
                self.nonimg_cat_features = nonimg_cat_features
                self.nonimg_con_features = nonimg_con_features

                nonimg_cat_features_selector = ("nonimg_cat", "passthrough", nonimg_cat_features)  # for now, assume that categorical features have already been one-hot-encoded
                nonimg_con_features_selector = ("nonimg_cont", myStandardizer(), nonimg_con_features)
                self.nonimg_transform = ColumnTransformer([nonimg_cat_features_selector, nonimg_con_features_selector], verbose_feature_names_out=False)
            else:
                self.nonimg_cat_features = None
                self.nonimg_con_features = None
                self.nonimg_transform = None

            # define amyloid PET feature selectors
            if "a" in self.features_to_extract:
                self.amyloid_transform = ColumnTransformer([
                    ("amyloid", myStandardizer(), make_column_selector(pattern = ".*\\.amyloid"))
                ], verbose_feature_names_out=False)
            else:
                self.amyloid_transform = None
            
            # define volume feature selectors
            if "v" in self.features_to_extract:
                self.volume_transform = ColumnTransformer([
                    ("volume",  myStandardizer(), make_column_selector(pattern = ".*\\.vol"))
                ], verbose_feature_names_out=False)
            else:
                self.volume_transform = None

            # fit to data
            self.fit()

        # multiwindow target variable
        if multiwindow:
            self.y_agg = self.y
            self.y = torch.stack([self.get_y_vector(p, t, s, num_time_windows, use_time_stable_for_progressors=True) for p,t,s in zip(self.y==1, time_to_progression, time_stable)])

    def fit(self, df = None, y = None):

        if df is None:
            df = self.df
        if y is None:
            y = self.y

        # fit combat
        if not self.combat_model is None:
            self.combat_model.fit(df)
            df = self.combat_model.transform(df)  # harmonize data first before fitting column transformers

        # fit column transformers
        for t in self.get_transform_list():
            if not t is None:
                t.fit(df, y.detach().cpu().numpy())

        # after fitting, store transformed data in class
        self.features = self.transform()

    def transform(self, df = None, y = None, dtype = None, device = None):

        if df is None:
            df = self.df
        if dtype is None:
            dtype = self.dtype
        if device is None:
            device = self.device

        # apply combat
        if not self.combat_model is None:
            df = self.combat_model.transform(df)
        
        # apply column transformers
        return [torch.tensor(t.transform(df).astype(np.float32), dtype=dtype, device=device) for t in self.get_transform_list()]

    def get_y_vector(
        self,
        is_progressor: bool,
        time_to_progression: float,
        time_stable: float,
        num_time_windows: int,
        use_time_stable_for_progressors: bool = False
    ) -> torch.Tensor:
        # progressor y vector
        if is_progressor:
            if time_to_progression > num_time_windows:
                return torch.zeros(num_time_windows, device = self.device, dtype = self.dtype)
            else:
                arr = np.arange(num_time_windows) + 1
                arr = np.where(
                    arr >= time_to_progression,
                    1,
                    np.where(arr <= time_stable, 0, np.nan) if use_time_stable_for_progressors else 0
                )
                return torch.tensor(arr, device = self.device, dtype = self.dtype)
        # stable y vector
        else:
            arr = np.arange(num_time_windows) + 1
            arr = np.where(arr <= time_stable, 0, np.nan)
            return torch.tensor(arr, device = self.device, dtype = self.dtype)
    
    def get_transform_list(self):
        return [t for t in [self.amyloid_transform, self.volume_transform, self.nonimg_transform] if not t is None]

    def __len__(self):
        return len(self.df)
    
    def __getitem__(self, idx):
        return [x[idx] for x in self.features], self.y[idx]


def get_stable_progressor(
    df,
    max_time_to_progression = 3,
    min_time_stable = 3,
    recover_progressor = False,
    progressor_can_be_stable = False,
    label_col = ".progressor",
    amyloid_positive_col = "amyloid_positive",
    valid_baseline_col = ".is_valid_baseline",
    time_to_progression_col = ".time_to_progression",
    time_stable_col = ".time_stable",
    target_name = "y",
):
    """
    Filter PET dataframe to get stable and progressors given a
    maximum time to progression and minimum time of stability

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe containing PET data
    max_time_to_progression : float (default=3)
        maximum time to progression to be considered progressor
    min_time_stable : float (default=3)
        minimum time of stability to be considered stable
    recover_progressor : bool (default=False)
        if True, include recovered progressors (subjects who
        reverted back to cognitively normal but then go back to
        MCI or dementia)
    progressor_can_be_stable : bool (default=False)
        if True, progressors that don't convert within
        max_time_to_progression but stay stable for
        min_time_stable are labeled as stable; if False, then
        only select the cleanest stable cohort possible (i.e.
        CN in entire follow-up history)
    label_col : str
        column name containing stable/progressor labels
    amyloid_positive_col : str
        column name containing bool column indicating amyloid
        positivity of PET scan
    valid_baseline_col : str
        column name containing bool column indicating whether
        scan counts as a valid baseline (i.e. closest matching
        clinical diagnosis is CN)
    time_to_progression_col : str
        column name containing times to progression
    time_stable_col : str
        column name containing time of stability
    target_name : str
        name of new column containing target variable, where
        0 = stable and 1 = progressor

    Returns
    -------
    pandas.DataFrame
        filtered dataframe of stable and progressor subjects
    """

    neg_class = ["stable"]
    if recover_progressor:
        pos_class = ["progressor", "progressor_recovered"]
    else:
        pos_class = ["progressor"]

    # first set of filtering
    df_filter = df[
        (df[valid_baseline_col] &    # scan starts at CN
        df[amyloid_positive_col])    # scan is amyloid-positive
    ]

    # get progressors
    df_progressor = df_filter[
        (df_filter[label_col].isin(pos_class) &   # positive class
        (df_filter[time_to_progression_col] <= max_time_to_progression))  # max time to progression
    ]

    # get stables
    df_stable = df_filter[
        (df_filter[label_col].isin(neg_class) &   # positive class
        (df_filter[time_stable_col] >= min_time_stable))  # min time of stability
    ]

    # get additional stables from eventual progressors
    if progressor_can_be_stable:
        df_stable_additional = df_filter[
            (df_filter[label_col].isin(pos_class) &   # positive class
            (df_filter[time_to_progression_col] > max_time_to_progression) &   # progression is past max time to progression
            (df_filter[time_to_progression_col] < 5) &   # FOR REVISION: ensure that stable cohorts at the 5-year models are the same regardless of this option
            (df_filter[time_stable_col] >= min_time_stable))  # min time of stability
        ]
        df_stable = pd.concat([df_stable, df_stable_additional])

    # create new column
    df_progressor[target_name] = 1
    df_stable[target_name] = 0

    return pd.concat([df_progressor, df_stable])


def get_stable_progressor_alt(
    df,
    max_time_to_progression = 3,
    min_time_stable = 3,
    amyloid_positive_col = "amyloid_positive",
    valid_baseline_col = ".is_valid_baseline",
    time_to_progression_col = ".time_to_progression",
    time_stable_col = ".time_stable",
    target_name = "y",
):

    """
    alternate function where stables are allowed to be labeled as progressors

    this function skips the negative/positive class check (i.e. stables have
    to be stable for the entire history of their CDR followups) and instead
    just filters using the stable-time and time-to-progression
    """

    # neg_class = ["stable"]
    # if recover_progressor:
    #     pos_class = ["progressor", "progressor_recovered"]
    # else:
    #     pos_class = ["progressor"]

    # first set of filtering
    df_filter = df[
        (df[valid_baseline_col] &    # scan starts at CN
        df[amyloid_positive_col])    # scan is amyloid-positive
    ]

    # get progressors
    df_progressor = df_filter[
        (df_filter[time_to_progression_col] <= max_time_to_progression)  # max time to progression
    ]

    # get stables
    df_stable = df_filter[
        (df_filter[time_stable_col] >= min_time_stable) &  # min time of stability
        ~df_filter.index.isin(df_progressor.index)
    ]

    # # get additional stables from eventual progressors
    # if progressor_can_be_stable:
    #     df_stable_additional = df_filter[
    #         (df_filter[label_col].isin(pos_class) &   # positive class
    #         (df_filter[time_to_progression_col] > max_time_to_progression) &   # progression is past max time to progression
    #         (df_filter[time_to_progression_col] < 5) &   # FOR REVISION: ensure that stable cohorts at the 5-year models are the same regardless of this option
    #         (df_filter[time_stable_col] >= min_time_stable))  # min time of stability
    #     ]
    #     df_stable = pd.concat([df_stable, df_stable_additional])

    # create new column
    df_progressor[target_name] = 1
    df_stable[target_name] = 0

    return pd.concat([df_progressor, df_stable])


# class DataLoader:

#     def __init__(
#         self,
#         data_dir: str,
#         progressor_time_window: str = "5",
#         stable_min_time: str = "5",
#         progressor_baseline: str = None,
#         filter_baseline: str = None
#     ):

#         self.data_dir = data_dir
#         self.progressor_time_window = progressor_time_window.strip()
#         self.stable_min_time = stable_min_time.strip()
#         self.progressor_baseline = progressor_baseline
#         self.filter_baseline = filter_baseline

#     @staticmethod
#     def _merge_progressor(data, subj_list, progressor_col, time_to_progression_col):

#         merge_df = subj_list.loc[:, ["subj", progressor_col, time_to_progression_col]]
#         data = data.merge(
#             right = merge_df,
#             how = "left",
#             on = "subj",
#             suffixes = (None, None)
#         )
#         data = data.rename(columns = {progressor_col: "progressor", time_to_progression_col: "time_to_progression"})

#         return data

#     @staticmethod
#     def _read_rds(rds_path):

#         # https://stackoverflow.com/questions/40996175/loading-a-rds-file-in-pandas

#         from pyreadr import read_r

#         df = read_r(rds_path)[None]

#         return df

# class PACDataLoader(DataLoader):

#     """
#     Data loader class for PAC tidied data

#     Attributes
#     ----------
#     data_dir : str
#         path to data directory
#     progressor_time_window : int, str (default = 5)
#         time window length for marking progressor; used to select
#         stable/progressor subj list
#     stable_min_time : int, str (default = 5)
#         minimum time of stability for marking stable; used to select
#         stable/progressor subj list
#     progressor_baseline : str (default = "pac")
#         which baseline was used to mark stable/progressor; used to
#         select stable/progressor subj list
#     filter_baseline : str (default = "pac")
#         which baseline is used to slice cross-sectional observation
#         of each subject; if "pac", the closest visit to PAC baseline
#         is used (that occurs after baseline); else, the subject's
#         earliest PET observation is used
        

#     Notes
#     -----
#     - `progressor_baseline` and `filter_baseline` are not the same
#       thing. `progressor_baseline` is used to determine the stable
#       and progressor subjects, while `filter_baseline` is used to
#       select the cross-sectional row for each subject. These can
#       technically be two different baselines, although it is
#       recommended that they be kept the same
#     """

#     def __init__(
#         self,
#         data_dir: str,
#         progressor_time_window: str = "5",
#         stable_min_time: str = "5",
#         progressor_baseline: str = "pac",
#         filter_baseline: str = "pac",
#         pet_analysis: str = "suvr"
#     ):

#         super().__init__(data_dir, progressor_time_window, stable_min_time, progressor_baseline, filter_baseline)
#         self.pet_analysis = pet_analysis

#         # load data
#         self.data = self._read_rds(os.path.join(self.data_dir, "features/features.RDS"))
#         self.subj_list = self._read_rds(os.path.join(self.data_dir, "subj_list/progressor.RDS"))
        
#     def get_Xy(
#         self,
#         remove_na = False,
#         split_Xy = False,
#         longitudinal = False,
#         time_to_progression = False
#     ):

#         """
#         Select cross-sectional features, merge stable/progressor target
#         labels, and return as X and y

#         Parameters
#         ----------
#         remove_na : bool (default = False)
#             whether to remove rows which contain NA in any column of X
#         split_Xy : bool (default = False)
#             whether to separate X and y
#         longitudinal : bool (default = False)
#             if True, return all longitudinal data without selecting baseline
#         time_to_progression : bool (default = False)
#             if True, return time to progression (continuous) instead of
#             stable/progressor (binary); for stables, time to progression
#             will be NA

#         Returns
#         -------
#         X : np.ndarray, shape (m, n)
#             array of features, where m is the number of cross-sectional
#             observations and n is the number of features + other relevant
#             columns
#         y : np.ndarray, shape (m, 1)
#             array of binary labels
#         """

#         # filter PET analysis
#         Xy = self.data.copy()
#         Xy = Xy[Xy["pet_analysis"] == self.pet_analysis]

#         # select baseline
#         if not longitudinal:

#             if self.filter_baseline == "pac":

#                 _df = Xy.loc[:, ("subj", "PAC_VisitNo")].copy()
#                 _df["idx"] = _df["PAC_VisitNo"]
#                 _df.loc[_df["idx"] < 0, "idx"] = np.nan  # ignore negative PAC visits

#                 min_idx = _df.groupby("subj")["idx"].idxmin()  # get row idx of smallest positive visit no
#                 min_idx = min_idx[~min_idx.isna()]  # ignore if no min was found (i.e. all negative PAC visits)

#             else:

#                 min_idx = _df.groupby("subj")["age"].idxmin()

#             Xy = Xy.loc[min_idx].copy()

#         # merge subject list
#         suffix = f"time_window.{str(self.progressor_time_window)}__stable_min_time.{str(self.stable_min_time)}__bl.{str(self.progressor_baseline)}"
#         Xy = self._merge_progressor(
#             Xy,
#             self.subj_list,
#             f"progressor__{suffix}",
#             f"time_to_progression__{suffix}"
#         )

#         # change dtype of progressor col
#         Xy["progressor"] = pd.to_numeric(Xy["progressor"])

#         # filter out NA target labels
#         Xy.dropna(axis = 0, subset = "progressor", inplace = True)

#         # remove unnecesary columns
#         drop_col = Xy.loc[0, "diagnosis":"mean cortical"].index.values
#         Xy.drop(columns = np.concatenate([["pet_analysis", "PAC_VisitNo"], drop_col]), inplace = True)

#         # remove NA in X (ignore target variables)
#         if remove_na:
#             col_drop = [c for c in Xy.columns if not c in ("progressor", "time_to_progression")]
#             Xy.dropna(axis = 0, inplace = True, subset = col_drop)

#         # split into X, y
#         if split_Xy:
#             y = Xy.pop("progressor")
#             t = Xy.pop("time_to_progression")

#             if time_to_progression:
#                 return Xy, t.values
#             else:
#                 return Xy, y.values

#         return Xy


# class ADNIOASISDataLoader(DataLoader):

#     """
#     Data loader class for ADNI and OASIS data.
    
#     Loads subject baseline data from ADNI and OASIS datasets.
#     This loader depends on data processed by R scripts in
#     `BradenADLongitudinalPrediction`

#     Parameters
#     ----------
#     features_path : str
#         path to `features.csv`
#     from_RDS : bool
#         if True, load data from RDS
#     """

#     def __init__(
#         self,
#         data_dir: str,
#         progressor_time_window: str = "5",
#         stable_min_time: str = "5",
#         progressor_baseline: str = "amyloid_positive",
#         filter_baseline: str = "amyloid_positive"
#     ):

#         super().__init__(data_dir, progressor_time_window, stable_min_time, progressor_baseline, filter_baseline)

#         # load data
#         self.data = self._read_rds(os.path.join(self.data_dir, "adni_oasis/adni_oasis_features.RDS"))
#         self.subj_list = self._read_rds(os.path.join(self.data_dir, "subj_list/adni_oasis_progressor.RDS"))

#     def get_Xy(
#         self,
#         remove_na = False,
#         split_Xy = False,
#         split_adni_oasis = False,
#         longitudinal = False,
#         time_to_progression = False
#     ):

#         """
#         Select cross-sectional features, merge stable/progressor target
#         labels, and return as X and y

#         Parameters
#         ----------
#         remove_na : bool (default = False)
#             whether to remove rows which contain NA in any column of X
#         split_Xy : bool (default = False)
#             whether to separate X and y
#         split_adni_oasis : bool (default = False)
#             whether to split ADNI and OASIS
#         longitudinal : bool (default = False)
#             if True, return all longitudinal data without selecting baseline
#         time_to_progression : bool (default = False)
#             if True, return time to progression (continuous) instead of
#             stable/progressor (binary); for stables, time to progression
#             will be NA

#         Returns
#         -------
#         X : np.ndarray, shape (m, n)
#             array of features, where m is the number of cross-sectional
#             observations and n is the number of features + other relevant
#             columns
#         y : np.ndarray, shape (m, 1)
#             array of binary labels
#         """

#         Xy = self.data.copy()

#         # get AV45+ (baseline) rows
#         if not longitudinal:

#             idx = Xy \
#                 .query("tracer == 'AV45' and amyloid_positive == True") \
#                 .groupby("subj") \
#                 ["age"].idxmin() \
#                 .dropna()  # remove NA age
#             Xy = Xy.loc[idx, :]

#         # remove unnecessary columns (subj_id, cdr, etc.)
#         col_drop = ["mmse", "cdr", "amyloid_positive", "left.cerebellum.cortex.pet", "right.cerebellum.cortex.pet", "icv.fs"]
#         Xy = Xy.drop(labels=col_drop, axis=1)

#         # merge subject list
#         suffix = f"time_window.{str(self.progressor_time_window)}__stable_min_time.{str(self.stable_min_time)}"
#         Xy = self._merge_progressor(
#             Xy,
#             self.subj_list,
#             f"progressor__{suffix}",
#             f"time_to_progression__{suffix}"
#         )

#         # change dtype of progressor col
#         Xy["progressor"] = pd.to_numeric(Xy["progressor"])

#         # filter out NA target labels
#         Xy.dropna(axis = 0, subset = "progressor", inplace = True)

#         # remove NA in X (ignore target variables)
#         if remove_na:
#             col_drop = [c for c in Xy.columns if not c in ("progressor", "time_to_progression")]
#             Xy.dropna(axis = 0, inplace = True, subset = col_drop)
        
#         # split datasets
#         if split_adni_oasis:
            
#             adni_df = Xy[Xy["study"] == "ADNI"]
#             oasis_df = Xy[Xy["study"] == "OASIS"]

#             # split into X, y
#             if split_Xy:
                
#                 adni_y = adni_df.pop("progressor")
#                 oasis_y = oasis_df.pop("progressor")
#                 adni_t = adni_df.pop("time_to_progression")
#                 oasis_t = oasis_df.pop("time_to_progression")
                
#                 if time_to_progression:
#                     return adni_df, oasis_df, adni_t.values, oasis_t.values
#                 else:
#                     return adni_df, oasis_df, adni_y.values, oasis_y.values
            
#             return adni_df, oasis_df

#         else:

#             # split into X, y
#             if split_Xy:
#                 y = Xy.pop("progressor")
#                 t = Xy.pop("time_to_progression")
#                 if time_to_progression:
#                     return Xy, t.values
#                 else:
#                     return Xy, y.values
        
#             return Xy
