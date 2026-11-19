import unittest
import numpy as np
from sklearn.linear_model import Ridge
from pdv3.image_calibration import donor_weights,constrain_composition,fit_grid,outer_fold


class ImageCalibrationTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(82)
        self.x=rng.normal(size=(32,5));self.y=rng.dirichlet([2,3,1],size=32)
        self.groups=np.repeat(['a','b','c','d'],[5,7,9,11])

    def test_donors_have_equal_total_training_weight(self):
        w=donor_weights(self.groups)
        np.testing.assert_allclose([w[self.groups==g].sum() for g in np.unique(self.groups)],8)

    def test_raw_predictions_match_independent_sklearn_solver(self):
        model=fit_grid(self.x,self.y,self.groups,[13.])[13.]
        z=(self.x-model.mean)/model.scale
        other=Ridge(alpha=13.,fit_intercept=True,solver='cholesky').fit(z,self.y,sample_weight=donor_weights(self.groups))
        _,raw=model.predict(self.x)
        np.testing.assert_allclose(raw,other.predict(z),atol=1e-12)

    def test_holdout_targets_cannot_change_own_predictions_or_alpha(self):
        first=outer_fold(self.x,self.y,self.groups,'a',alphas=[1.,10.])
        changed=self.y.copy();changed[self.groups=='a']=np.roll(changed[self.groups=='a'],1,axis=1)
        second=outer_fold(self.x,changed,self.groups,'a',alphas=[1.,10.])
        np.testing.assert_array_equal(first['prediction'],second['prediction'])
        self.assertEqual(first['selected_alpha'],second['selected_alpha'])
        self.assertNotIn('a',first['training_donors'])
        for row in first['inner']:
            self.assertNotIn('a',row['training_donors'])
            self.assertNotIn(row['inner_test_donor'],row['training_donors'])

    def test_constant_features_return_donor_balanced_baseline(self):
        model=fit_grid(np.ones_like(self.x),self.y,self.groups,[1.])[1.]
        prediction,_=model.predict(np.ones((2,5)))
        expected=np.mean([self.y[self.groups==g].mean(0) for g in np.unique(self.groups)],axis=0)
        np.testing.assert_allclose(prediction,np.tile(expected,(2,1)),atol=1e-12)

    def test_nonnegative_constraint_and_zero_fallback(self):
        p=constrain_composition(np.array([[-1.,2.,2.],[-2.,-1.,-3.]]),np.array([.2,.3,.5]))
        np.testing.assert_allclose(p,[[0,.5,.5],[.2,.3,.5]])


if __name__=='__main__':unittest.main()
