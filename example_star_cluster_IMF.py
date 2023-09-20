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
                              "    This code generates the stellar masses of one star-cluster given the total "
                              "star-cluster mass applying optimal sampling.\n\n"
                              "    Please type in the cluster mass in the solar mass unit then hit return:"))
M_over_H = float(input("\n    The code assumes an empirical relation between the IMF slopes for low-mass stars and metallicity.\n"
                       "    Canonical IMF is recovered with solar metallicity, i.e., [M/H]=0.\n"
                       "    Please type in the initial metallicity of the cluster, [M/H], then hit return to sample stellar masses:"))

age = float(input("\n    The code calculate the mass of the most massive star at a given age according to PARSEC v1.2 stellar evolution model.\n"
                  "    Currently, the age resolution is 10 Myr.\n"
                  "    Please type in the age of the star cluster in [yr], then hit return to sample stellar masses (input 0 if age is not a concern):"))


# Calculate lifetime according to (mass, metallicity)
Z_list_value = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.008, 0.01, 0.02, 0.03, 0.04]
Z_list_index = np.argmin(np.abs(np.array(Z_list_value) - 0.02*10**M_over_H))
stellar_Z_extrapolated = Z_list_value[Z_list_index]
if stellar_Z_extrapolated == 0.0001:
    data_AGB = np.loadtxt('Mass_lifetime_relation/PARSEC/Mass_lifetime_relation_Z_0.0001.txt')
elif stellar_Z_extrapolated == 0.0002:
    data_AGB = np.loadtxt('Mass_lifetime_relation/PARSEC/Mass_lifetime_relation_Z_0.0002.txt')
elif stellar_Z_extrapolated == 0.0005:
    data_AGB = np.loadtxt('Mass_lifetime_relation/PARSEC/Mass_lifetime_relation_Z_0.0005.txt')
elif stellar_Z_extrapolated == 0.001:
    data_AGB = np.loadtxt('Mass_lifetime_relation/PARSEC/Mass_lifetime_relation_Z_0.001.txt')
elif stellar_Z_extrapolated == 0.002:
    data_AGB = np.loadtxt('Mass_lifetime_relation/PARSEC/Mass_lifetime_relation_Z_0.002.txt')
elif stellar_Z_extrapolated == 0.004:
    data_AGB = np.loadtxt('Mass_lifetime_relation/PARSEC/Mass_lifetime_relation_Z_0.004.txt')
elif stellar_Z_extrapolated == 0.008:
    data_AGB = np.loadtxt('Mass_lifetime_relation/PARSEC/Mass_lifetime_relation_Z_0.008.txt')
elif stellar_Z_extrapolated == 0.01:
    data_AGB = np.loadtxt('Mass_lifetime_relation/PARSEC/Mass_lifetime_relation_Z_0.01.txt')
elif stellar_Z_extrapolated == 0.02:
    data_AGB = np.loadtxt('Mass_lifetime_relation/PARSEC/Mass_lifetime_relation_Z_0.02.txt')
elif stellar_Z_extrapolated == 0.03:
    data_AGB = np.loadtxt('Mass_lifetime_relation/PARSEC/Mass_lifetime_relation_Z_0.03.txt')
elif stellar_Z_extrapolated == 0.04:
    data_AGB = np.loadtxt('Mass_lifetime_relation/PARSEC/Mass_lifetime_relation_Z_0.04.txt')

def function_mass_boundary(this_time, data_AGB):
    logAge_mass_boundary = np.round(data_AGB[:, 0], 5)
    logAge_value = np.log10(this_time)
    logAge_list_value = np.round(sorted(set(logAge_mass_boundary)), 5)
    logAge_list_index = np.argmin(np.abs(np.array(logAge_list_value) - np.round(logAge_value, 5)))
    logAge_value_extrapolated = logAge_list_value[logAge_list_index]
    index = np.where((logAge_mass_boundary == np.round(logAge_value_extrapolated, 5)))
    index = index[0]
    AGB_mass_boundary = 10**data_AGB[index, 2]
    star_mass_boundary = 10**data_AGB[index, 3]
    return AGB_mass_boundary, star_mass_boundary

(AGB_mass_boundary, star_mass_boundary) = function_mass_boundary(age, data_AGB)
print("    The most massive star with the given age (with a time resolution of 10 Myr) and metallicity "
      "can have an initial mass of {} solar mass, according to PARSCE.".format(star_mass_boundary))


# setup alpha values:
alpha_2 = 2.3
alpha_1 = 1.3
alpha3_model = 2
alpha2_model = 'Z'  # or 1 for our publications before 2020
alpha1_model = 'Z'  # or 1 for our publications before 2020
alpha3_change = galimf.function_alpha_3_change(alpha3_model, StarClusterMass, M_over_H)
alpha2_change = galimf.function_alpha_2_change(alpha_2, alpha2_model, M_over_H)
alpha1_change = galimf.function_alpha_1_change(alpha_1, alpha1_model, M_over_H)


print("\n    - Sampling the star cluster with {} solar mass and [M/H] = {} -".format(StarClusterMass, M_over_H))

# apply galIMF to optimally sample stars from IMF:
galimf.function_sample_from_imf(StarClusterMass, 1, 0.08, alpha1_change, 0.5, alpha2_change, 1, alpha3_change, 150)

# apply galIMF to draw IMF analytically:
galimf.function_draw_xi_str(0.08, StarClusterMass, 1, 0.08, alpha1_change, 0.5, alpha2_change, 1, alpha3_change, 150)
List_M_str_for_xi_str = galimf.x_IMF
List_xi_str = galimf.y_IMF

print("\n    - Sampling completed -\n")
# followings are all sampled results:

# most massive stellar mass in the cluster:
print("    The most massive star formed in this star cluster has {} solar mass.".format(round(galimf.list_M_str_i[0], 2)))

length_list = len(galimf.list_M_str_i)
i__ = length_list
while i__ > 0 and galimf.list_M_str_i[0] > star_mass_boundary:
    del galimf.list_M_str_i[0]
    del galimf.list_n_str_i[0]
    (i__) = (i__-1)

print("    The most massive star still alive at {} Myr has {} solar mass.".format(age/1e6, round(galimf.list_M_str_i[0], 2)))

# All of the sampled stellar masses in the solar mass unit are (from massive to less massive):
list_stellar_masses = np.array(galimf.list_M_str_i)

# # The bolometric luminosity is estimated according to Yan et al. 2019, 2022:
# L_bol_tot = 0
# for mass in list_stellar_masses:
#     log_mass = math.log(mass, 10)
#     if  log_mass < -0.37571790416:  # < log0.421
#         log_L_bol = 2.3 * log_mass -0.63827216398
#     elif log_mass < 0.29225607135:
#         log_L_bol = 4 * log_mass
#     elif log_mass < 1.74358815016:
#         log_L_bol = 3.5 * log_mass + 0.14612803567
#     else:
#         log_L_bol = log_mass + 4.50514997832
#     L_bol_tot += 10**log_L_bol
# print("    The total (ZAMS) bolometric luminosity of all the optimally-sampled stars is estimated to be: {} L_sun.".format(round(L_bol_tot, 2)))

# NOTE! Multiple stars can be represented by the same stellar mass if they have similar masses,
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
print("\n    Stellar masses of every star still alive at {} Myr in the star cluster are saved in the file: Stellar_masses_for_a_star_cluster.txt".format(age/1e6))

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

# Constructing the canonical IMF:
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
