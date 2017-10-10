# Example -- Optimally Sample an Embedded star Cluster.
# An example Python 3 code which use galIMF.py to sample the stellar masses in one star cluster with optimal sampling.

# Made by: Zhiqiang & Tereza

import galIMF  # Main part of the GalIMF code for generating and sampling Galaxy-wide stellar Initial Mass Function.
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scipy.integrate import quad

# -----------------------------------------------------------------------
# figure output settings:
fig0 = plt.figure(figsize=(4, 3))  # size for one column plot
gs1 = GridSpec(1, 1)
ax0 = plt.subplot(gs1[0])

# -----------------------------------------------------------------------


StarClusterMass = 1.e4
alpha3_model = 1
Fe_over_H = 0

print("\n - Sampling One star cluster with: -")
print("mass = {} solar mass;".format(StarClusterMass))
print("alpha3_model = {} (see Function_alpha_3_change in the file 'galIMF.py' for details);".format(alpha3_model))
print("[Fe/H] = {}.\n".format(Fe_over_H))


alpha3_change = galIMF.Function_alpha_3_change(1, StarClusterMass, 0)  # read in alpha_3 value

# Run galIMF and sample stars from IMF:
galIMF.function_sample_from_IMF(StarClusterMass, 1, 0.08, 1.3, 0.5, 2.3, 1, alpha3_change, 150)

# Sampling results:
print(" - The most massive stellar mass in solar mass unit in this star cluster is: -")
print(galIMF.list_M_str_i[0])

# The stellar masses in solar mass unit are (from massive to less massive):
list_stars = np.array(galIMF.list_M_str_i)

# NOTE! Multiple stars can be represented by a same stellar mass if they have similar masses,
# The number of stars represented by the stellar masses above are:
n_stars = np.array(galIMF.list_n_str_i)

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

ax0.set_ylabel(r'$\log_{\rm 10}(\xi, [\#_{\star}/M_\odot])$')
ax0.set_xlabel(r'$\log_{\rm 10}(m, [M_\odot])$')

plt.legend()
plt.tight_layout()
plt.savefig('cluster_optimal_sample.pdf', dpi=300)

file = open('stellar_masses_in_a_star_cluster.txt', 'w')
file.write("# Optimally sampled stellar masses in a star cluster with:\n")
file.write("# mass = {} solar mass\n".format(StarClusterMass))
file.write("# alpha3_model = {}\n".format(alpha3_model))
file.write("# [Fe/H] = {}\n".format(Fe_over_H))
file.write("# The stellar masses are:\n")
for item in galIMF.list_M_str_i:
    file.write("%s\n" % item)
file.close()

print("\n - Sampling completed -")
print(" - The results are saved in file: "
      "'stellar_masses_in_a_star_cluster.txt' and 'cluster_optimal_sample.pdf' -")
