import pbeis
import pybamm
import numpy as np
import matplotlib.pyplot as plt
import scipy

# Set up model
model = pybamm.lithium_ion.DFN()
parameter_values = pybamm.ParameterValues("Chen2020")

# Compute impedances
frequencies = np.logspace(-4, 2, 30)
eis_sim = pbeis.EISSimulation(model, parameter_values=parameter_values)
impedances_freq = eis_sim.solve(frequencies, 'direct')

plt.ion()
ax = pbeis.nyquist_plot(impedances_freq, alpha=0.7, )
plt.suptitle(f"{model.name}")
scipy.io.savemat("pybamm_chen_impedances.mat", {"impedances": impedances_freq})
