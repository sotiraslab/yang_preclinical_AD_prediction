import unittest

import numpy as np

from ..model import StandardScalerControlGroup


class TestStandardScalerControlGroup(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        
        cls.X = np.array([[0.35, 0.19, 0.26],
                          [0.92, 0.82, 0.33],
                          [0.62, 0.14, 0.49],
                          [0.62, 0.7 , 0.01],
                          [0.85, 0.26, 0.48],
                          [0.07, 0.95, 0.67],
                          [0.2 , 0.75, 0.52],
                          [0.13, 0.55, 0.45],
                          [0.54, 0.28, 0.62],
                          [0.  , 0.03, 0.53]])
        cls.y = np.array([1,1,0,0,1,0,0,0,0,0])

    def test_standardization(self):

        m = np.mean(self.X, axis=0)
        s = np.sqrt(np.var(self.X[self.y==0], axis=0))
        X_std1 = (self.X - m) / s

        sscg = StandardScalerControlGroup()
        X_std2 = sscg.fit_transform(self.X, self.y)

        self.assertTrue(np.isclose(X_std1, X_std2).all())
