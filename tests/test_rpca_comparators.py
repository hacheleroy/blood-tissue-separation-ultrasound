"""Small complex-valued regression checks for the three comparison arms."""
import unittest
import numpy as np

from btbench.rpca_comparators import stable_rpca, cauchy_rpca_core, batch_irls_rpca


class TestRpcaComparators(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(47)
        tissue = np.outer(rng.normal(size=12) + 1j*rng.normal(size=12),
                          np.exp(1j*np.linspace(0, 1, 16)))
        blood = np.zeros((12, 16), complex)
        blood[2:5] = 0.4*(rng.normal(size=(3, 16)) + 1j*rng.normal(size=(3, 16)))
        self.iq = (tissue + blood).T.reshape(16, 3, 4)

    def test_global_complex_phase_and_scale_equivariance(self):
        methods = [(stable_rpca, {"tau": 0.1, "sparse_lambda": 0.03}),
                   (cauchy_rpca_core, {"sparse_lambda": 0.1, "gamma": 0.1}),
                   (batch_irls_rpca, {"low_rank_lambda": 0.1, "sparse_mu": 0.01,
                                      "rank": 3, "epsilon": 0.02})]
        scalar = 2.7*np.exp(1j*0.62)
        for method, params in methods:
            with self.subTest(method=method.__name__):
                y, _ = method(self.iq, max_iter=12, **params)
                rotated, _ = method(scalar*self.iq, max_iter=12, **params)
                self.assertEqual(y.shape, self.iq.shape)
                self.assertTrue(np.isfinite(y).all())
                self.assertLess(np.linalg.norm(rotated-scalar*y)/max(np.linalg.norm(scalar*y), 1e-30), 1e-8)

    def test_irls_one_step_complex_normal_equations(self):
        iq = self.iq
        x = iq.reshape(len(iq), -1).T
        scale = np.sqrt(np.mean(abs(x)**2))
        x = x/scale
        rank, lam, mu, epsilon = 3, 0.1, 0.01, 0.02
        left, singular, vh = np.linalg.svd(x, full_matrices=False)
        u = left[:, :rank]*np.sqrt(singular[:rank])
        v = vh[:rank].conj().T*np.sqrt(singular[:rank])
        s = (x-u@v.conj().T)/(1+mu/epsilon)
        d = 1/np.sqrt(np.sum(abs(u)**2, axis=0)+np.sum(abs(v)**2, axis=0)+epsilon**2)
        gram_u = u.conj().T@u+lam*np.diag(d)
        v_new = (x-s).conj().T@u@np.linalg.inv(gram_u)
        gram_v = v_new.conj().T@v_new+lam*np.diag(d)
        u_new = (x-s)@v_new@np.linalg.inv(gram_v)
        self.assertLess(np.linalg.norm(v_new@gram_u-(x-s).conj().T@u), 1e-10)
        self.assertLess(np.linalg.norm(u_new@gram_v-(x-s)@v_new), 1e-10)
        output, diag = batch_irls_rpca(iq, lam, mu, rank, epsilon, max_iter=1)
        np.testing.assert_allclose(output.reshape(len(iq), -1).T/scale, s, atol=1e-10)
        self.assertEqual(diag["iterations"], 1)


if __name__ == "__main__":
    unittest.main()
