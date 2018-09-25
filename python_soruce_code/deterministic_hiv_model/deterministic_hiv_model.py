"""
    Here we reproduce the numerical solution of the HIV-AIDS model
    reported in [1].

    Using the same notation as in [1], the model reads:
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
                k_2 J - (\mu + d) A.
        \end{aligned}
    \end{equation}
    [1]
        Stability analysis of an HIV/AIDS epidemic model with treatment.
        Liming Cai, Xuezhi Li, Mini Ghoshc, Baozhu Guod.
        Journal of Computational and Applied Mathematics.
"""
import numpy as np
from scipy import integrate


class DeterministicHIVAIDSMODEL(object):

    def __init__(self, t_0=0.0, t_f=300.0,
                 b=0.3, c=3.0, d=0.0,
                 alpha=0.01, beta=0.0005,
                 mu=0.02, k_1=0.01, k_2=0.02
                 ):
        # Parameters for the test example
        self.t_0 = t_0
        self.t_f = t_f
        #
        self.b = b
        self.c = c
        self.d = d
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
        self.r_zero = num / den

    def set_parameters(self, b, c, d, alpha, beta, mu, k_1, k_2):
        self.b = b
        self.c = c
        self.d = d
        self.alpha = alpha
        self.beta = beta
        self.mu = mu
        self.k_1 = k_1
        self.k_2 = k_2
        self.r_zero = 0.0

    def compute_r_zero(self):
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
        self.r_zero = r_zero
        return r_zero

    #
    @staticmethod
    def f_rhs(x, t, b, c, d, alpha, beta, mu, k_1, k_2, n_total):
        s = x[0]
        i = x[1]
        j = x[2]
        a = x[3]
        f_s = mu * n_total - c * beta * (i + b * j) * s - mu * s
        f_i = c * beta * (i + b * j) * s - (mu + k_1) * i + alpha * j
        f_j = k_1 * i - (mu + k_2 + alpha) * j
        f_a = k_2 * j - (mu + d) * a
        del t
        rhs = np.array([f_s, f_i, f_j, f_a])
        return rhs


class DeterministicHIVAIDSMODELNumerics(DeterministicHIVAIDSMODEL):

    def __init__(self, eps=0.0001, n_max=1000, dynamic_dim=4):
        super(DeterministicHIVAIDSMODELNumerics, self).__init__()
        self.n_max = n_max
        self.eps = eps
        self.dynamic_dim = dynamic_dim
        self.t = np.linspace(self.t_0, self.t_f, n_max)
        self.x = np.zeros([n_max, dynamic_dim])

    def lsoda_solution(self):
        # lsoda solution from scipy library
        t = self.t
        y_0 = [self.s_zero, self.i_zero, self.j_zero, self.a_zero]
        b = self.b
        c = self.c
        d = self.d
        alpha = self.alpha
        beta = self.beta
        mu = self.mu
        k_1 = self.k_1
        k_2 = self.k_2
        n_total = self.n_total
        #
        y = integrate.odeint(self.f_rhs, y_0, t,
                             args=(b, c, d, alpha, beta,
                                   mu, k_1, k_2, n_total))
        self.x = y
        return y
