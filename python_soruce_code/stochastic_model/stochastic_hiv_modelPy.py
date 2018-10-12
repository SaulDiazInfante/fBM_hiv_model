import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from stochastic_hiv_model import StochasticHIVAIDSMODELNumerics

#
#
hiv_aids_solver = StochasticHIVAIDSMODELNumerics()
ou = hiv_aids_solver.mean_ou_solution()
plt.plot(ou)
plt.show()
