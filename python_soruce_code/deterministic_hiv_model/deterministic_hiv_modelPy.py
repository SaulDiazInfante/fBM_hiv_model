import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from deterministic_hiv_model import DeterministicHIVAIDSMODELNumerics

#
#
hiv_aids_solver = DeterministicHIVAIDSMODELNumerics()
x = hiv_aids_solver.lsoda_solution()
t = hiv_aids_solver.t
s = x[:, 0]
i = x[:, 1]
j = x[:, 2]
a = x[:, 3]
#
#
fig = plt.figure()
gs = GridSpec(3, 3, figure=fig)
#
#
ax0 = fig.add_subplot(gs[0, 2])
plt.plot(t, s)
ax0.set_xlabel(r'$t$')
ax0.set_ylabel(r'Suceptibles $S(t)$')
#
#
ax1 = fig.add_subplot(gs[:, 0: 2])
ax1.set_ylabel(r'Infected populations $I(t), J(t)$')
plt.plot(t, i, label='Asymptomatic $I(t)$')
plt.plot(t, j, label='Symptomatic $J(t)$')
plt.legend(loc=0)
ax2 = fig.add_subplot(gs[2, 2])
plt.plot(t, a)
ax2.set_ylabel(r'AIDS $A(t)$')
plt.tight_layout()
plt.show()
