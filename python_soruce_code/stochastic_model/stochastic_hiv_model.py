"""
    Here we implement numerical methods to approximate the solution of our
    stochastic fBM HIV-AIDS model---which is a stochastic extension of the
    HIV-AIDS model reported in [1].

    Using the same notation as in [1], our model reads:
    \begin{equation}
        \begin{aligned}
            \frac{dS}{dt} &=
                \mu K - c \beta (I + b J) S - \mu S,
                \\
            \frac{dI}{dt} &=
                c \beta (I + b J) S - (\mu + k_1)  I + \alpha J,
                \\
            \frac{dJ}{dt} &=
                k_1 I - (\mu + k_2 + \alpha) J,
                \\
            \frac{dA}{dt} &=
                k_2 J - (\mu + mean_mu_d) A.

        \end{aligned}
    \end{equation}
    [1]
        Stability analysis of an HIV/AIDS epidemic model with treatment.
        Liming Cai, Xuezhi Li, Mini Ghoshc, Baozhu Guod.
        Journal of Computational and Applied Mathematics.
"""
import numpy as np
from fbm import FBM


class StochasticHIVAIDSMODEL(object):

    def __init__(self, t_0=0.0, t_f=300.0,
                 b=0.3, c=3.0, mean_mu_d=1.5,
                 alpha=0.01, beta=0.0005,
                 mu=0.02, k_1=0.01, k_2=0.02,
                 gamma_d=1.0, sigma_d=1.0, hurst=0.5,
                 mu_d_zero=3.0
                 ):
        # Parameters for the test example
        self.t_0 = t_0
        self.t_f = t_f
        #
        self.b = b
        self.c = c
        self.mean_mu_d = mean_mu_d
        self.alpha = alpha
        self.beta = beta
        self.mu = mu
        #
        self.k_1 = k_1
        self.k_2 = k_2
        #  Initial conditions
        self.s_zero = 150.0
        self.i_zero = 30.0
        self.j_zero = 20.0
        self.a_zero = 0.0
        self.n_total = self.s_zero + self.i_zero + self.j_zero + self.a_zero
        num = c * beta * self.n_total * (mu + k_2 + alpha + b * k_1)
        den = (mu + k_1) * (mu + k_2) * (mu * alpha)
        self.r_zero_d = num / den
        #
        self.gamma_d = gamma_d
        self.mu_d_zero = mu_d_zero
        self.sigma_d = sigma_d
        self.hurst = hurst

    def set_parameters(self, b, c, mean_mu_d, alpha, beta, mu, k_1, k_2):
        self.b = b
        self.c = c
        self.mean_mu_d = mean_mu_d
        self.alpha = alpha
        self.beta = beta
        self.mu = mu
        self.k_1 = k_1
        self.k_2 = k_2
        self.r_zero_d = 0.0

    def compute_det_r_zero(self):
        b = self.b
        c = self.c
        alpha = self.alpha
        beta = self.beta
        mu = self.mu
        k_1 = self.k_1
        k_2 = self.k_2
        n_total = self.n_total

        num = c * beta * n_total * (mu + k_2 + alpha + b * k_1)
        den = (mu + k_1) * (mu + k_2) * (mu * alpha)
        r_zero = num / den
        self.r_zero_d = r_zero
        return r_zero

    #
    @staticmethod
    def f_rhs(x, t, b, c, mu_d, alpha, beta, mu, k_1, k_2, n_total):
        s = x[0]
        i = x[1]
        j = x[2]
        a = x[3]
        f_s = mu * n_total - c * beta * (i + b * j) * s - mu * s
        f_i = c * beta * (i + b * j) * s - (mu + k_1) * i + alpha * j
        f_j = k_1 * i - (mu + k_2 + alpha) * j
        f_a = k_2 * j - (mu + mu_d) * a
        del t
        rhs = np.array([f_s, f_i, f_j, f_a])
        return rhs

    def ou_drift(self, x):
        gamma_d = self.gamma_d
        mean_mu_d = self.mean_mu_d
        f_ou = gamma_d * (mean_mu_d - x)
        return f_ou


class StochasticHIVAIDSMODELNumerics(StochasticHIVAIDSMODEL):

    def __init__(self, eps=0.0001, n_max=2 ** 12, dynamic_dim=4):
        super(StochasticHIVAIDSMODELNumerics, self).__init__()
        self.n_max = n_max
        self.eps = eps
        self.dynamic_dim = dynamic_dim
        self.t = np.linspace(self.t_0, self.t_f, n_max)
        self.x = np.zeros([n_max, dynamic_dim])
        self.mu_d = np.zeros([n_max, 1])

    def fbm(self, n, hurst, length=1, method="daviesharte"):
        """One off sample of fBm."""
        f = FBM(n, hurst, length, method)
        return f.fbm()

    def mean_ou_solution(self):
        n_max = self.n_max
        x_j = self.mu_d_zero
        self.mu_d[0] = x_j
        mu_d = self.mu_d
        gamma_d = self.gamma_d
        mean_mu_d = self.mean_mu_d
        h = self.t[1] - self.t[0]
        sigma_d = self.sigma_d
        gamma_d_mean_mu_d = gamma_d * mean_mu_d
        hurst = self.hurst
        fbm_generator = self.fbm(n=n_max, hurst=0.5, length=1,
                                 method="daviesharte")
        for j in np.arange(n_max - 1):
            x_j = mu_d[j]
            x_jp1 = gamma_d * mean_mu_d + np.exp(-gamma_d * h) \
                    * (x_j - gamma_d_mean_mu_d)
            mu_d[j + 1] = x_jp1 + \
                          sigma_d * (fbm_generator[j + 1] - fbm_generator[j])
        self.mu_d = mu_d
        return mu_d
