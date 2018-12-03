# Example -- Optimally Sample an Embedded star Cluster.
# An example Python 3 code which use galIMF.py to sample the stellar masses in one star cluster with optimal sampling.

# Made by: Zhiqiang & Tereza
# -----------------------------------------------------------------------

import galIMF  # Main part of the GalIMF code for generating and sampling Galaxy-wide stellar Initial Mass Function.
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scipy.integrate import quad


# figure output settings:
fig0 = plt.figure(figsize=(4, 3))  # size for one column plot
gs1 = GridSpec(1, 1)
ax0 = plt.subplot(gs1[0])

# input parameters:
StarClusterMass = float(input("\n    ============================\n"
                              "    === example_star_cluster ===\n"
                              "    ============================\n\n"
                              "    This code generate the stellar masses of one star cluster with given mass"
                              " applying optimal sampling.\n\n"
                              "    Please type in the cluster mass in solar mass unit then hit return:"))
alpha3_model = 1
Fe_over_H = float(input("    Please type in the initial [Fe/H] then hit return:"))

print("\n    - Sampling the star cluster with {} solar mass and [Fe/H] = {} -".format(StarClusterMass, Fe_over_H))

# setup alpha_3 value:
alpha3_change = galIMF.Function_alpha_3_change(1, StarClusterMass, 0)

# apply galIMF to optimally sample stars from IMF:
galIMF.function_sample_from_IMF(StarClusterMass, 1, 0.08, 1.3, 0.5, 2.3, 1, alpha3_change, 150)

print("\n    - Sampling completed -")
# followings are all sampled results:

# most massive stellar mass in the cluster:
print("    The most massive star in this star cluster has {} solar mass".format(round(galIMF.list_M_str_i[0])))

# All of the sampled stellar masses in solar mass unit are (from massive to less massive):
list_stars = np.array(galIMF.list_M_str_i)

# NOTE! Multiple stars can be represented by a same stellar mass if they have similar masses,
# The number of stars represented by the stellar masses above are:
list_orgen = galIMF.list_n_str_i
if list_orgen[-1] == 0:
    del list_orgen[-1]
n_stars = np.array(list_orgen)

# formating a figure output to compare the optimally sampled result (label: OS) with canonical IMF (label: IMF):

# bining the sampled star number:
bins = np.logspace(np.log10(0.08), np.log10(150), 20, base=10)
vals0 = np.zeros(len(bins))

for i, b in enumerate(bins):
    if i == len(bins)-1:
        break
    else:
        star_array = (list_stars[np.logical_and(list_stars >= b, list_stars < bins[i+1])])
        n_array = (n_stars[np.logical_and(list_stars >= b, list_stars < bins[i+1])])
        len_array = 0
        for j, n in enumerate(n_array):
            len_array = len_array+n
        vals0[i] = len_array/(bins[i+1]-bins[i])

ax0.step(np.log10(bins), np.log10(vals0+1.e-3), color='blue', where='post', zorder=1, lw=1.5, label="OS")

# constracting the canonical IMF:
N = 100
can_imf = np.zeros(N)
masses = np.logspace(np.log10(0.08), np.log10(150), N, base=10)

for i, m in enumerate(masses):
    if m <= 0.5:
        can_imf[i] = m ** (-1.3)
    else:
        can_imf[i] = 0.5*m ** (-2.3)


def imf(mass, k, alpha):
    return k*mass*mass**(-alpha)

Norm = quad(imf, 0.08, 0.5, args=(1, 1.3))[0] + quad(imf, 0.5, 120, args=(0.5, 2.3))[0]
can_imf = np.array(can_imf)*StarClusterMass/Norm
ax0.plot(np.log10(masses), np.log10(can_imf), color='black', lw=1.5, label='IMF', zorder=0)

# plot settings:
ax0.set_ylabel(r'$\log_{\rm 10}(\xi, [\#_{\star}/M_\odot])$')
ax0.set_xlabel(r'$\log_{\rm 10}(m, [M_\odot])$')
plt.legend()
plt.tight_layout()

# save the plot:
plt.savefig('cluster_optimal_sample.pdf', dpi=300)

# end of the example:
print("    The sampling results are plotted in file: 'cluster_optimal_sample.pdf'\n\n"
      "    ============================\n")

# show the plot
plt.show()
