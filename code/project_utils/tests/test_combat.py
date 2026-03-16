import unittest

import numpy as np
import pandas as pd

import neuroHarmonize
from ..combat import GAMComBat


# TODO: test for robustness of missing batch labels

class TestGAMComBat(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.df, cls.data, cls.covar, cls.sites = cls._create_dummy_data(random_state = 100)
        cls.df_test, cls.data_test, cls.covar_test, cls.sites_test = cls._create_dummy_data(random_state = 101)

    def _create_dummy_data(
        n_row=100,
        n_features=4,
        n_covar=2,
        n_site=3,
        random_state = None
    ):
        if not random_state is None:
            np.random.seed(random_state)
        
        data = np.random.rand(n_row, n_features)
        covar = np.random.rand(n_row, n_covar)
        sites = np.random.choice([f"site{i}" for i in range(1, n_site+1)], size=n_row)
        df = pd.DataFrame(np.concat([data, covar], axis=1), columns=[f"feature{i}" for i in range(1, n_features+1)] + [f"covar{i}" for i in range(1, n_covar+1)])
        df["SITE"] = sites

        return df, data, covar, sites

    def test_combat_linear(self):

        # manually train ComBat
        model, _ = neuroHarmonize.harmonizationLearn(
            self.data,
            self.df.loc[:,("SITE", "covar1", "covar2")],
            seed = 42
        )
        result1 = neuroHarmonize.harmonizationApply(
            self.data_test,
            self.df_test.loc[:,("SITE", "covar1", "covar2")],
            model
        )

        # use GAMComBat class to train
        gam_combat = GAMComBat(
            features=["feature1", "feature2", "feature3", "feature4"],
            covariates=["covar1", "covar2"],
            site_name = "SITE",
            random_state = 42
        )
        gam_combat.fit(self.df)
        result2 = gam_combat.transform(self.df_test).loc[:, ["feature1", "feature2", "feature3", "feature4"]].values

        # check for equality
        self.assertTrue(np.allclose(result1, result2))

    def test_combat_gam(self):

        # TODO: for some reason, results are not deterministic even when seed is set; open issue on GitHub

        # manually train ComBat
        model, _ = neuroHarmonize.harmonizationLearn(
            self.data,
            self.df.loc[:,("SITE", "covar1", "covar2")],
            smooth_terms = ["covar1"],
            smooth_term_bounds = (None, None),
            seed = 42
        )
        result1 = neuroHarmonize.harmonizationApply(
            self.data_test,
            self.df_test.loc[:,("SITE", "covar1", "covar2")],
            model
        )

        # use GAMComBat class to train
        gam_combat = GAMComBat(
            features=["feature1", "feature2", "feature3", "feature4"],
            covariates=["covar1", "covar2"],
            site_name = "SITE",
            smooth_terms = ["covar1"],
            random_state = 42
        )
        gam_combat.fit(self.df)
        result2 = gam_combat.transform(self.df_test).loc[:, ["feature1", "feature2", "feature3", "feature4"]].values

        # check for equality
        self.assertTrue(np.allclose(result1, result2))

    def test_orig_data(self):

        gam_combat = GAMComBat(
            features=["feature1", "feature2"],
            covariates=["covar1", "covar2"],
            site_name = "SITE",
            random_state = 42
        )
        gam_combat.fit(self.df)
        df_test_harm = gam_combat.transform(self.df_test)

        # check that column order is preserved
        self.assertListEqual(df_test_harm.columns.tolist(), self.df_test.columns.tolist())

        # check that the original features are not modified
        self.assertTrue(np.allclose(df_test_harm.loc[:,["feature3", "feature4"]].values, self.df_test.loc[:,["feature3", "feature4"]].values))
