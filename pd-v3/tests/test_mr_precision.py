import unittest

import numpy as np
from scipy.stats import chi2, ncx2, norm

from pdv3.mr_precision import conditional_power, detectable_effect, covariance_delta_se, bounded_fieller_interval


class PrecisionTests(unittest.TestCase):
    def test_normal_power_against_noncentral_chisquare(self):
        for alpha in [.05, .05/11655183, 1e-12]:
            for effect in [-3., -.2, 0., .2, 3.]:
                actual = conditional_power(effect, .13, alpha=alpha)
                expected = ncx2.sf(chi2.isf(alpha, 1), 1, (effect/.13)**2)
                self.assertAlmostEqual(actual, expected, places=12)
            self.assertAlmostEqual(conditional_power(0, 1., alpha=alpha), alpha, delta=alpha*1e-12)

    def test_detectable_effect_inversion_and_units(self):
        se = np.array([.0001, .03, 1., 1000.])
        for alpha in [.05, .05/11655183]:
            for target in [.8, .9]:
                mde = detectable_effect(se, target_power=target, alpha=alpha)
                np.testing.assert_allclose(conditional_power(mde, se, alpha=alpha), target, atol=1e-12)
                np.testing.assert_allclose(detectable_effect(se*100, target_power=target, alpha=alpha), mde*100)

    def test_covariance_zero_matches_existing_delta(self):
        bx, sx, by, sy = -.2, .01, .08, .006
        expected = np.sqrt(sy**2/bx**2 + by**2*sx**2/bx**4)
        self.assertAlmostEqual(covariance_delta_se(bx, sx, by, sy), expected, places=14)
        self.assertEqual(covariance_delta_se(.2, .01, .1, .005, rho=1), 0.)

    def test_fieller_independent_polynomial_roots(self):
        rng = np.random.default_rng(991)
        for _ in range(100):
            bx = rng.choice([-1, 1])*rng.uniform(.1, 1)
            sx = abs(bx)/rng.uniform(7, 50)
            by, sy = rng.normal(0, .1), rng.uniform(.005, .04)
            for rho in [-1, -.5, 0, .5, 1]:
                lower, upper = bounded_fieller_interval(bx, sx, by, sy, rho=rho)
                q = norm.isf(.025)**2
                roots = np.sort(np.roots([bx*bx-q*sx*sx, -2*bx*by+2*q*rho*sx*sy, by*by-q*sy*sy]))
                np.testing.assert_allclose([lower, upper], roots, atol=1e-12, rtol=1e-10)
                for theta in [lower, upper]:
                    self.assertAlmostEqual((by-theta*bx)**2,
                        q*(sy*sy+theta*theta*sx*sx-2*theta*rho*sx*sy), places=12)
                self.assertEqual(bool(lower <= 0 <= upper), bool(abs(by/sy) <= norm.isf(.025)))

    def test_fieller_allele_reversal_and_scale(self):
        base = bounded_fieller_interval(.2, .01, -.08, .006, rho=.5)
        np.testing.assert_allclose(base, bounded_fieller_interval(-.2, .01, .08, .006, rho=.5))
        np.testing.assert_allclose(np.array(base)/2, bounded_fieller_interval(.4, .02, -.08, .006, rho=.5))
        np.testing.assert_allclose([-base[1], -base[0]], bounded_fieller_interval(.2, .01, .08, .006, rho=-.5))

    def test_refuse_invalid_or_unbounded_inputs(self):
        for args in [(0., .01, .1, .01), (.1, -.01, .1, .01), (np.nan, .01, .1, .01), (.01, .02, .1, .01)]:
            with self.assertRaises(ValueError):
                bounded_fieller_interval(*args)
        with self.assertRaises(ValueError):
            covariance_delta_se(.2, .01, .1, .01, rho=1.01)
        for alpha in [0, 1, np.nan]:
            with self.assertRaises(ValueError):
                conditional_power(.1, .1, alpha=alpha)
        with self.assertRaises(ValueError):
            detectable_effect(.1, target_power=.01)


if __name__ == '__main__':
    unittest.main()
