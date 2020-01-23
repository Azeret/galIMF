# Python3 code, last update Wed 20 Dec 2018

# Example file for sampling the stellar masses of every star in the star cluster.

# Made by: Yan Zhiqiang & Tereza Jerabkova

# The outputs of this example are:

#  - a comparison plot of generated variable IMF and canonical IMF ('star_cluster_IMF_plot.pdf');
#  - a .txt file containing the stellar masses ('Stellar_masses_for_a_star_cluster.txt').

# --------------------------------------------------------------------------------------------------------------------------------
# Import modules and libraries
# --------------------------------------------------------------------------------------------------------------------------------

import galimf  # Main part of the GalIMF code for generating and sampling Galaxy-wide stellar Initial Mass Function.
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import csv  # csv and izip/zip are used to create output files
try:
    from itertools import izip as zip
except ImportError:  # will be python 3.x series
    pass

# -----------------------------------------------------------------------

# figure output settings:
fig0 = plt.figure(figsize=(4, 3))  # size for one column plot
gs1 = GridSpec(1, 1)
ax0 = plt.subplot(gs1[0])

# input parameters:
StarClusterMass = float(input("\n    ================================\n"
                              "    === example_star_cluster_IMF ===\n"
                              "    ================================\n\n"
                              "    This code generate the stellar masses of one star-cluster given the total "
                              "star-cluster mass applying optimal sampling.\n\n"
                              "    Please type in the cluster mass in solar mass unit then hit return:"))
M_over_H = float(input("    Please type in the initial metallicity, [M/H], then hit return:"))

print("\n    - Sampling the star cluster with {} solar mass and [M/H] = {} -".format(StarClusterMass, M_over_H))

# setup alpha values:
alpha_2 = 2.3
alpha_1 = 1.3
alpha3_model = 'R14'  # 2
alpha2_model = 'R14'  # 1
alpha1_model = 0  # 1
alpha3_change = galimf.function_alpha_3_change(alpha3_model, StarClusterMass, M_over_H)
alpha2_change = galimf.function_alpha_2_change(alpha_2, alpha2_model, M_over_H)
alpha1_change = galimf.function_alpha_1_change(alpha_1, alpha1_model, M_over_H)

# apply galIMF to optimally sample stars from IMF:
galimf.function_sample_from_imf(StarClusterMass, 1, 0.08, alpha1_change, 0.5, alpha2_change, 1, alpha3_change, 150)

# apply galIMF to draw IMF analytically:
galimf.function_draw_xi_str(0.08, StarClusterMass, 1, 0.08, alpha1_change, 0.5, alpha2_change, 1, alpha3_change, 150)
List_M_str_for_xi_str = galimf.x_IMF
List_xi_str = galimf.y_IMF

print("\n    - Sampling completed -\n")
# followings are all sampled results:

# most massive stellar mass in the cluster:
print("    The most massive star in this star cluster has {} solar mass".format(round(galimf.list_M_str_i[0], 2)))

# All of the sampled stellar masses in solar mass unit are (from massive to less massive):
list_stellar_masses = np.array(galimf.list_M_str_i)

# NOTE! Multiple stars can be represented by a same stellar mass if they have similar masses,
# The number of stars represented by the stellar masses above are:
list_stellar_numbers = galimf.list_n_str_i
if list_stellar_numbers[-1] == 0:
    del list_stellar_numbers[-1]
n_stars = np.array(list_stellar_numbers)

# save the sampled stellar mass in a txt file:
with open('Stellar_masses_for_a_star_cluster.txt', 'w') as file:
    writer = csv.writer(file, delimiter=' ')
    file.write(
        "# Output file of the generated stellar masses for a star cluster with given mass and metallicity.\n"
        "# The columns are:\n# Mass in solar mass unit; "
        "Number of stars in this star cluster have mass close to this value\n\n")
    writer.writerows(
        zip(list_stellar_masses, list_stellar_numbers))
print("\n    Stellar masses of every star in the star cluster is saved in the file: "
      "Stellar_masses_for_a_star_cluster.txt\n")

# formatting a figure output to compare the optimally sampled result (label: OS) with canonical IMF (label: IMF):

# binning the sampled star number:
bins = np.logspace(np.log10(0.08), np.log10(150), 20, base=10)
vals0 = np.zeros(len(bins))

for i, b in enumerate(bins):
    if i == len(bins)-1:
        break
    else:
        star_array = (list_stellar_masses[np.logical_and(list_stellar_masses >= b, list_stellar_masses < bins[i+1])])
        n_array = (n_stars[np.logical_and(list_stellar_masses >= b, list_stellar_masses < bins[i+1])])
        len_array = 0
        for j, n in enumerate(n_array):
            len_array = len_array+n
        vals0[i] = len_array/(bins[i+1]-bins[i])

ax0.step(np.log10(bins), np.log10(vals0+1.e-3), color='blue', where='post', zorder=1, lw=1.5, label="Optimally sampled stellar masses")

# constructing the canonical IMF:
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


Norm = quad(imf, 0.08, 0.5, args=(1, 1.3))[0] + quad(imf, 0.5, 150, args=(0.5, 2.3))[0]
can_imf = np.array(can_imf)*StarClusterMass/Norm
ax0.plot(np.log10(masses), np.log10(can_imf), color='black', lw=1.5, label='Canonical IMF', zorder=0, ls='dotted')

# plot analytical IMF:
ax0.plot(np.log10(List_M_str_for_xi_str), np.log10(List_xi_str), color='red', label='analytical IMF', zorder=0, ls='dashed')


# plot settings:
ax0.set_ylabel(r'$\log_{\rm 10}(\xi, [\#_{\star}/M_\odot])$')
ax0.set_xlabel(r'$\log_{\rm 10}(m, [M_\odot])$')
plt.legend()
plt.tight_layout()

# save the plot:
plt.savefig('star_cluster_IMF_plot.pdf', dpi=300)

# end of the example:
print("    The plot is saved in the file: star_cluster_IMF_plot.pdf\n\n"
      "    ============================\n")

# show the plot
plt.show()



### Plot the mmax--Mecl relation:
#
# alpha_2 = 2.3
# alpha_1 = 1.3
# alpha3_model = 2
# alpha2_model = 1
# alpha1_model = 1
# M_over_H_list = [-3, -2, -1, 0, 1]
# for j in range(5):
#     M_over_H = M_over_H_list[j]
#     # The lowest possible star cluster mass has a limit when M_ecl=m_max, depending on the assumed IMF.
#     if M_over_H < 0.1:
#         lower_cluster_mass_limit = 0.15 - M_over_H / 11
#     else:
#         lower_cluster_mass_limit = 0.15 - M_over_H / 28
#     alpha2_change = galimf.function_alpha_2_change(alpha_2, alpha2_model, M_over_H)
#     alpha1_change = galimf.function_alpha_1_change(alpha_1, alpha1_model, M_over_H)
#     StarClusterMass_list = np.arange(lower_cluster_mass_limit, 10, 0.1).tolist()
#     for i in range(49):
#         StarClusterMass_list += [10 ** ((i + 10) / 10)]
#
#     M_max_list = []
#
#     for i in range(len(StarClusterMass_list)):
#         StarClusterMass = StarClusterMass_list[i]
#         alpha3_change = galimf.function_alpha_3_change(alpha3_model, StarClusterMass, M_over_H)
#         galimf.function_sample_from_imf(StarClusterMass, 1, 0.08, alpha1_change, 0.5, alpha2_change, 1, alpha3_change,
#                                         150)
#         galimf.function_draw_xi_str(0.08, StarClusterMass, 1, 0.08, alpha1_change, 0.5, alpha2_change, 1, alpha3_change,
#                                     150)
#         List_M_str_for_xi_str = galimf.x_IMF
#         List_xi_str = galimf.y_IMF
#         M_max_list.append(galimf.list_M_str_i[0])
#     plt.loglog(StarClusterMass_list, M_max_list, label="[Z]={}".format(M_over_H))
#
# plt.loglog([5, 5], [0.1, 10], lw=0.5, label=r'cluster mass = 5 [M$_\odot$]')
# plt.loglog([0.1, 10], [1, 1], ls='dotted', c='0.5')
#
# plt.xlabel(r"Star cluster mass [M$_\odot$]")
# plt.ylabel(r"Most massive star mass [M$_\odot$]")
# plt.legend()
# plt.tight_layout()
# plt.show()
