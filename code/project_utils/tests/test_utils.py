import unittest

import torch

from ..utils import MaskedMultiwindowBCEWithLogitsLoss


class TestMaskedMultiwindowBCEWithLogitsLoss(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        
        cls.y_true = torch.Tensor([[0,0,0,0,1],[0,0,torch.nan,1,1]])
        cls.y_pred = torch.Tensor([[0.1,0.2,0.3,0.4,0.9],[0.2,0.4,0.6,0.7,0.9]])
        cls.pos_weight = torch.Tensor([1,2,3,4,5])
        cls.weight = torch.Tensor([2,4,6,8,10])

    def test_outputs_equal(self):

        mmbce = MaskedMultiwindowBCEWithLogitsLoss()
        loss1 = mmbce(self.y_pred, self.y_true)

        pos_log_sigmoid = torch.log(torch.sigmoid(self.y_pred))
        neg_log_sigmoid = torch.log(1 - torch.sigmoid(self.y_pred))
        loss2 = -torch.concatenate([pos_log_sigmoid[self.y_true == 1], neg_log_sigmoid[self.y_true == 0]]).mean()

        self.assertTrue(torch.isclose(loss1, loss2))
    
    def test_outputs_equal_pos_weight(self):

        mmbce = MaskedMultiwindowBCEWithLogitsLoss(pos_weight = self.pos_weight)
        loss1 = mmbce(self.y_pred, self.y_true)

        pos_log_sigmoid = torch.log(torch.sigmoid(self.y_pred)) * self.pos_weight
        neg_log_sigmoid = torch.log(1 - torch.sigmoid(self.y_pred))
        loss2 = -torch.concatenate([pos_log_sigmoid[self.y_true == 1], neg_log_sigmoid[self.y_true == 0]]).mean()

        self.assertTrue(torch.isclose(loss1, loss2))

    def test_outputs_equal_weight(self):

        mmbce = MaskedMultiwindowBCEWithLogitsLoss(weight = self.weight)
        loss1 = mmbce(self.y_pred, self.y_true)

        pos_log_sigmoid = torch.log(torch.sigmoid(self.y_pred)) * self.weight
        neg_log_sigmoid = torch.log(1 - torch.sigmoid(self.y_pred)) * self.weight
        loss2 = -torch.concatenate([pos_log_sigmoid[self.y_true == 1], neg_log_sigmoid[self.y_true == 0]]).mean()

        self.assertTrue(torch.isclose(loss1, loss2))
    
    def test_outputs_equal_weight_pos_weight(self):

        mmbce = MaskedMultiwindowBCEWithLogitsLoss(weight = self.weight, pos_weight = self.pos_weight)
        loss1 = mmbce(self.y_pred, self.y_true)

        pos_log_sigmoid = torch.log(torch.sigmoid(self.y_pred)) * self.weight * self.pos_weight
        neg_log_sigmoid = torch.log(1 - torch.sigmoid(self.y_pred)) * self.weight
        loss2 = -torch.concatenate([pos_log_sigmoid[self.y_true == 1], neg_log_sigmoid[self.y_true == 0]]).mean()

        self.assertTrue(torch.isclose(loss1, loss2))