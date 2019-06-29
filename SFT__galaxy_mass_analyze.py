import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import griddata
import numpy as np
import math
import itertools
import random

canonical_IMF_assumption = "Kroupa"

# plot observations ###

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
fig = plt.figure(1, figsize=(4, 3.5))

plt.xlabel(r'log$_{10}$(M$_{*}$ [M$_\odot$])')
plt.ylabel(r'[Z/X] or [Mg/Fe]')


data_Arrigoni2010_dynamical_mass = [9.24, 9.18, 9.25, 9, 9.36, 9.22, 10.3, 8.32, 10.61, 11.85, 11.89, 10.69, 10.48, 11.16, 10.85, 10.27, 11.41, 10.02, 10.19, 10.5, 9.45, 10.59, 9.34, 10.26, 10.85, 9.21, 10.12, 11.91, 10.89, 9.74, 10.72, 11.45, 9.88, 11.84, 11.09, 11.13, 10.26, 10.39, 10.56, 9.9, 10.5, 9.7, 9.34, 11.03, 11.07, 10.78, 10.73, 9.71, 10.08, 10.33, 9.56, 10.09, 11.4, 10.72, 10.68, 10.72, 11.21, 10.63, 11.38, 10.39, 11.25, 11, 10.68, 11.54, 10.2, 11.29, 11.58, 11.36, 11.19]
data_Arrigoni2010_dynamical_mass_error = [0.43, 0.39, 0.48, 0.45, 0.46, 0.5, 0.33, 0.28, 0.15, 0.35, 0.34, 0.3, 0.34, 0.27, 0.28, 0.42, 0.22, 0.31, 0.4, 0.31, 0.45, 0.31, 0.4, 0.32, 0.28, 0.4, 0.39, 0.21, 0.32, 0.45, 0.29, 0.39, 0.67, 0.36, 0.32, 0.38, 0.94, 0.35, 0.32, 0.28, 0.28, 0.44, 0.49, 0.35, 0.33, 0.28, 0.35, 0.42, 0.41, 0.45, 0.56, 0.38, 0.28, 0.23, 0.29, 0.33, 0.26, 0.31, 0.23, 0.36, 0.41, 0.3, 0.32, 0.38, 0.31, 0.35, 0.3, 0.37, 0.39]
data_Arrigoni2010_Z_H = [-0.431, -0.732, -0.412, -0.266, -0.904, 0.357, 0.049, -0.398, 0.281, 0.316, 0.359, 0.318, 0.16, 0.448, 0.21, -0.084, 0.26, -0.367, -0.079, -0.117, -0.313, 0.114, 0.068, -0.02, 0.157, -0.012, 0.104, 0.346, 0.187, -0.506, -0.062, 0.294, 0.155, 0.409, 0.335, 0.404, 0.226, 0.243, 0.384, 0.038, 0.159, 0.332, 0.029, 0.247, 0.285, 0.224, 0.361, 0.311, 0.153, 0.418, 0.074, 0.173, 0.284, 0.027, 0.185, 0.403, 0.1, 0.36, 0.296, 0.215, 0.266, 0.439, 0.18, 0.179, -0.128, 0.203, 0.32, 0.336, 0.266]
data_Arrigoni2010_Z_H_error_p = [0.083, 0.159, 0.098, 0.053, 0.22, 0.083, 0.068, 0.038, 0.068, 0.038, 0.068, 0.023, 0.053, 0.114, 0.038, 0.053, 0.053, 0.053, 0.053, 0.053, 0.068, 0.068, 0.053, 0.053, 0.098, 0.114, 0.053, 0.068, 0.068, 0.098, 0.023, 0.038, 0.083, 0.053, 0.023, 0.053, 0.068, 0.053, 0.053, 0.038, 0.038, 0.053, 0.083, 0.068, 0.068, 0.053, 0.038, 0.083, 0.098, 0.083, 0.068, 0.053, 0.038, 0.053, 0.053, 0.038, 0.053, 0.023, 0.038, 0.114, 0.038, 0.038, 0.053, 0.038, 0.053, 0.038, 0.023, 0.023, 0.038]
data_Arrigoni2010_Z_H_error_m = [0.098, 0.144, 0.098, 0.114, 0.189, 0.098, 0.083, 0.038, 0.038, 0.023, 0.053, 0.023, 0.053, 0.083, 0.023, 0.053, 0.038, 0.053, 0.053, 0.038, 0.083, 0.114, 0.068, 0.053, 0.068, 0.098, 0.038, 0.038, 0.038, 0.189, 0.023, 0.023, 0.053, 0.023, 0.023, 0.023, 0.053, 0.053, 0.053, 0.038, 0.023, 0.038, 0.083, 0.053, 0.068, 0.038, 0.023, 0.098, 0.083, 0.038, 0.083, 0.038, 0.023, 0.038, 0.023, 0.098, 0.038, 0.008, 0.023, 0.114, 0.023, 0.023, 0.023, 0.038, 0.038, 0.008, 0.023, 0.008, 0.008]
# data_Arrigoni2010_Mg_Fe = [-0.005, -0.044, -0.076, 0.045, -0.164, 0.085, 0.162, -0.026, 0.153, 0.239, 0.179, 0.172, 0.093, 0.249, 0.137, 0.121, 0.109, 0.122, 0.173, 0.14, 0.111, 0.149, -0.045, 0.154, 0.127, -0.041, 0.092, 0.212, 0.128, 0.173, 0.104, 0.182, 0.134, 0.224, 0.114, 0.212, 0.096, 0.166, 0.095, 0.131, 0.168, 0.045, -0.011, 0.06, 0.098, 0.136, 0.186, 0.087, 0.088, 0.253, 0.109, 0.04, 0.175, 0.09, 0.162, 0.188, 0.197, 0.135, 0.215, 0.089, 0.217, 0.102, 0.11, 0.213, 0.031, 0.144, 0.173, 0.218, 0.151]
# data_Arrigoni2010_Mg_Fe_error_p = [0.066, 0.126, 0.076, 0.045, 0.187, 0.045, 0.025, 0.025, 0.025, 0.015, 0.035, 0.015, 0.015, 0.056, 0.015, 0.035, 0.025, 0.035, 0.035, 0.025, 0.035, 0.025, 0.035, 0.025, 0.035, 0.056, 0.025, 0.025, 0.025, 0.045, 0.015, 0.015, 0.025, 0.025, 0.015, 0.015, 0.025, 0.025, 0.025, 0.025, 0.015, 0.025, 0.045, 0.035, 0.035, 0.015, 0.015, 0.045, 0.056, 0.025, 0.045, 0.015, 0.015, 0.025, 0.015, 0.005, 0.015, 0.015, 0.015, 0.056, 0.015, 0.015, 0.015, 0.015, 0.025, 0.015, 0.005, 0.005, 0.015]
# data_Arrigoni2010_Mg_Fe_error_m = [0.045, 0.096, 0.056, 0.035, 0.5, 0.035, 0.015, 0.005, 0.005, 0.005, 0.015, 0.005, 0.015, 0.025, 0.005, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.025, 0.015, 0.015, 0.035, 0.005, 0.015, 0.005, 0.025, 0.005, 0.005, 0.015, 0.015, 0.005, 0.005, 0.015, 0.005, 0.025, 0.005, 0.005, 0.015, 0.035, 0.015, 0.025, 0.005, 0.005, 0.035, 0.035, 0.015, 0.035, 0.015, 0.005, 0.005, 0.005, 0.025, 0.005, 0.005, 0.005, 0.045, 0.005, 0.015, 0.005, 0.005, 0.015, 0.005, 0.005, 0.005, 0.005]


plt.errorbar(data_Arrigoni2010_dynamical_mass, data_Arrigoni2010_Z_H, xerr=[data_Arrigoni2010_dynamical_mass_error, data_Arrigoni2010_dynamical_mass_error],
             yerr=[data_Arrigoni2010_Z_H_error_m, data_Arrigoni2010_Z_H_error_p], capsize=1, elinewidth=0.3, capthick=0.5, fmt='none', c='0.5')
# plt.errorbar(data_Arrigoni2010_dynamical_mass, data_Arrigoni2010_Mg_Fe, xerr=[data_Arrigoni2010_dynamical_mass_error, data_Arrigoni2010_dynamical_mass_error],
#              yerr=[data_Arrigoni2010_Mg_Fe_error_m, data_Arrigoni2010_Mg_Fe_error_p], capsize=1, elinewidth=0.3, capthick=0.5, fmt='none', c='0.5')

plt.scatter(data_Arrigoni2010_dynamical_mass, data_Arrigoni2010_Z_H, marker='x', c='royalblue', label='Arrigoni10 [Z/X]', zorder=6)
# plt.scatter(data_Arrigoni2010_dynamical_mass, data_Arrigoni2010_Mg_Fe, marker='x', c='chocolate', label='Arrigoni10 [Mg/Fe]', zorder=6)


mass_generate_from_Arrigoni2010 = []
Z_H_generate_from_Arrigoni2010 = []
MgFe_generate_from_Arrigoni2010 = []
j = 0
while j < 1000:
    i = 0
    length = len(data_Arrigoni2010_dynamical_mass)
    while i < length:
        mass = data_Arrigoni2010_dynamical_mass[i]
        sigma_mass = data_Arrigoni2010_dynamical_mass_error[i]
        random_mass_1 = random.normalvariate(mass, sigma_mass)
        random_mass_2 = random.normalvariate(mass, sigma_mass)
        Z_H = data_Arrigoni2010_Z_H[i]
        sigma_Z_H_p = data_Arrigoni2010_Z_H_error_p[i]
        sigma_Z_H_m = data_Arrigoni2010_Z_H_error_m[i]
        random_error_Z_H_p = abs(random.normalvariate(0, sigma_Z_H_p))
        random_error_Z_H_m = abs(random.normalvariate(0, sigma_Z_H_m))
        random_Z_H_1 = Z_H + random_error_Z_H_p
        random_Z_H_2 = Z_H - random_error_Z_H_m
        mass_generate_from_Arrigoni2010.append(random_mass_1)
        mass_generate_from_Arrigoni2010.append(random_mass_2)
        Z_H_generate_from_Arrigoni2010.append(random_Z_H_1)
        Z_H_generate_from_Arrigoni2010.append(random_Z_H_2)

        # Mg_Fe = data_Arrigoni2010_Mg_Fe[i]
        # sigma_MgFe_p = data_Arrigoni2010_Mg_Fe_error_p[i]
        # sigma_MgFe_m = data_Arrigoni2010_Mg_Fe_error_m[i]
        # random_error_MgFe_p = abs(random.normalvariate(0, sigma_MgFe_p))
        # random_error_MgFe_m = abs(random.normalvariate(0, sigma_MgFe_m))
        # random_MgFe_1 = Mg_Fe + random_error_MgFe_p
        # random_MgFe_2 = Mg_Fe - random_error_MgFe_m
        # MgFe_generate_from_Arrigoni2010.append(random_MgFe_1)
        # MgFe_generate_from_Arrigoni2010.append(random_MgFe_2)
        (i) = (i+1)
    (j) = (j+1)
# plt.scatter(mass_generate_from_Arrigoni2010, Z_H_generate_from_Arrigoni2010, s=0.1)
# plt.scatter(mass_generate_from_Arrigoni2010, MgFe_generate_from_Arrigoni2010, s=0.1)

mass__list = []
mean_data_Z_H__list = []
up_data_Z_H__list = []
down_data_Z_H__list = []
mean_data_MgFe__list = []
up_data_MgFe__list = []
down_data_MgFe__list = []
# std_data_Z_H__list = []
# std_data_MgFe__list = []
mass__ = 8
while mass__ < 12.501:
    mass__upperlimit = mass__ + 0.5
    mass__lowerlimit = mass__ - 0.5
    data_mass__ = []
    data_Z_H__ = []
    data_MgFe__ = []
    i = 0
    length = len(mass_generate_from_Arrigoni2010)
    while i < length:
        if mass_generate_from_Arrigoni2010[i] < mass__upperlimit and mass_generate_from_Arrigoni2010[i] > mass__lowerlimit:
            data_Z_H__.append(Z_H_generate_from_Arrigoni2010[i])
            # data_MgFe__.append(MgFe_generate_from_Arrigoni2010[i])
        (i) = (i + 1)
    mean_data_Z_H__ = np.mean(data_Z_H__)
    # mean_data_MgFe__ = np.mean(data_MgFe__)
    std_data_Z_H__ = np.std(data_Z_H__)
    # std_data_MgFe__ = np.std(data_MgFe__)
    mass__list.append(mass__)
    mean_data_Z_H__list.append(mean_data_Z_H__)
    up_data_Z_H__list.append(mean_data_Z_H__ + std_data_Z_H__)
    down_data_Z_H__list.append(mean_data_Z_H__ - std_data_Z_H__)
    # mean_data_MgFe__list.append(mean_data_MgFe__)
    # up_data_MgFe__list.append(mean_data_MgFe__+std_data_MgFe__)
    # down_data_MgFe__list.append(mean_data_MgFe__-std_data_MgFe__)
    # std_data_Z_H__list.append(std_data_Z_H__)
    # std_data_MgFe__list.append(std_data_MgFe__)
    (mass__) = (mass__+0.01)

plt.plot(mass__list, mean_data_Z_H__list, c='k', label='Arrigoni10 mean', zorder=7)
plt.plot(mass__list, up_data_Z_H__list, c='k', ls='dotted', label='Arrigoni10 std', zorder=5)
plt.plot(mass__list, down_data_Z_H__list, c='k', ls='dotted', zorder=5)
plt.fill_between(mass__list, down_data_Z_H__list, up_data_Z_H__list, alpha=0.3, facecolor='royalblue', linewidth=0)

# plt.plot(mass__list, mean_data_MgFe__list, c='k', zorder=7)
# plt.plot(mass__list, up_data_MgFe__list, c='k', ls='dotted', zorder=5)
# plt.plot(mass__list, down_data_MgFe__list, c='k', ls='dotted', zorder=5)
# plt.fill_between(mass__list, down_data_MgFe__list, up_data_MgFe__list, alpha=0.3, facecolor='chocolate', linewidth=0)


mass_thomas = np.arange(6, 15, 1)
alpha_thomas = -0.459 + 0.062 * mass_thomas
alpha_thomas_up = -0.459 + 0.062 * mass_thomas + 0.1
alpha_thomas_low = -0.459 + 0.062 * mass_thomas - 0.1
plt.plot(mass_thomas, alpha_thomas, c='k', ls='dashed', zorder=7, label='Thomas05')
plt.plot(mass_thomas, alpha_thomas_up, c='k', ls='dotted', zorder=5)
plt.plot(mass_thomas, alpha_thomas_low, c='k', ls='dotted', zorder=5)
plt.fill_between(mass_thomas, alpha_thomas_low, alpha_thomas_up, alpha=0.3, facecolor='chocolate', linewidth=0)


############## fit discrete relation with a function ##############

mass_obs_array = np.array(mass__list)
metal_obs_array_high = np.array(up_data_Z_H__list)
metal_obs_array = np.array(mean_data_Z_H__list)
metal_obs_array_low = np.array(down_data_Z_H__list)
fun_obs_metal_mass_high = interpolate.interp1d(mass_obs_array, metal_obs_array_high, kind='cubic')
fun_obs_metal_mass = interpolate.interp1d(mass_obs_array, metal_obs_array, kind='cubic')
fun_obs_metal_mass_low = interpolate.interp1d(mass_obs_array, metal_obs_array_low, kind='cubic')

alpha_obs_array_high = np.array(up_data_MgFe__list)
alpha_obs_array = np.array(mean_data_MgFe__list)
alpha_obs_array_low = np.array(down_data_MgFe__list)
# fun_obs_alpha_mass_high = interpolate.interp1d(mass_obs_array, alpha_obs_array_high, kind='cubic')
# fun_obs_alpha_mass = interpolate.interp1d(mass_obs_array, alpha_obs_array, kind='cubic')
# fun_obs_alpha_mass_low = interpolate.interp1d(mass_obs_array, alpha_obs_array_low, kind='cubic')


# mass_obs       = [9, 10, 10.5, 11, 12]
# metal_obs_high = [-0.1, 0.1, 0.2, 0.25, 0.3]
# metal_obs      = [-0.8, -0.3, -0.02, 0.1, 0.15]
# metal_obs_low  = [-2.5, -1, -0.375, -0.1, 0]


# mass_obs       = [9, 10, 10, 10.3, 10.3, 10.5, 10.5, 10.8, 10.8, 11, 11, 12]
# metal_obs_high = [-0.35, -0.35, 0.1, 0.1, 0.2, 0.2, 0.4, 0.4, 0.41, 0.41, 0.3, 0.3]
# metal_obs      = []
# metal_obs_low  = [-1, -1, -0.75, -0.75, -0.5, -0.5, -0.2, -0.2, -0.08, -0.08, 0 ,0]
# for i in range(len(mass_obs)):
#     metal_obs.append((metal_obs_high[i]+metal_obs_low[i])/2)

# plt.plot(mass_obs, metal_obs_high, c='k', ls='dashed')
# plt.plot(mass_obs, metal_obs_low, c='k', ls='dashed')
# plt.plot(mass_obs, metal_obs, c='k')
# # plt.plot(mass_obs, metal_obs, 'k', color='royalblue', label='Luminosity-weighted stellar [Z] obs')
# # plt.fill_between(mass_obs, metal_obs_low, metal_obs_high, alpha=0.5, facecolor='royalblue', linewidth=0)


# mass_obs_array = np.array([9, 10, 10.5, 11, 12])
# metal_obs_array_high = np.array([-0.1, 0.1, 0.2, 0.25, 0.3])
# metal_obs_array = np.array([-0.8, -0.3, -0.02, 0.1, 0.15])
# metal_obs_array_low = np.array([-2.5, -1, -0.375, -0.1, 0])
# # mass_obs       = [0, 8, 9, 9.5, 10, 10.15, 10.3, 10.4, 10.5, 10.65, 10.8, 10.9, 11, 11.5, 12, 13, 20]
# # metal_obs_high = [-0.35, -0.35, -0.35, -0.35, (0.1-0.35)/2, 0.1, (0.3)/2, 0.2, (0.6)/2, 0.4, (0.81)/2, 0.41, (0.71)/2, 0.3, 0.3, 0.3, 0.3]
# # metal_obs      = []
# # metal_obs_low  = [-1, -1, -1, -1, (-1.75)/2, -0.75, -1.25/2, -0.5,  -0.7/2, -0.2, -0.28/2, -0.08, -0.08/2, 0, 0, 0, 0]
# # for i in range(len(mass_obs)):
# #     metal_obs.append((metal_obs_high[i]+metal_obs_low[i])/2)
# mass_obs_array = np.array(mass_obs)
# metal_obs_array_high = np.array(metal_obs_high)
# metal_obs_array = np.array(metal_obs)
# metal_obs_array_low = np.array(metal_obs_low)

# mass_new = np.arange(9, 12.01, 0.01)
# metal_new_high = fun_obs_metal_mass_high(mass_new)
# metal_new = fun_obs_metal_mass(mass_new)
# metal_new_low = fun_obs_metal_mass_low(mass_new)
# plt.plot(mass_new, metal_new_high, c='k', ls='dashed')
# plt.plot(mass_new, metal_new, c='k')
# plt.plot(mass_new, metal_new_low, c='k', ls='dashed')

# plt.plot(mass_new, metal_new, 'k', color='royalblue', label='Luminosity-weighted stellar [Z] obs')
# plt.fill_between(mass_new, metal_new_low, metal_new_high, alpha=0.5, facecolor='royalblue', linewidth=0)


# mass_obs  =     [ 9,    10,    11,   12]
# alpha_obs_high= [ 0.08, 0.205, 0.27, 0.32]
# alpha_obs =     [-0.01, 0.125, 0.2,  0.21]
# alpha_obs_low = [-0.1,  0.045, 0.13, 0.10]
# # mass_obs  =     [ 9, 9.5, 9.5, 10.5, 10.5, 11.5, 11.5, 12]
# # alpha_obs_high= [ 0.08, 0.08, 0.205, 0.205, 0.27, 0.27, 0.32, 0.32]
# # alpha_obs =     [-0.01, -0.01, 0.125, 0.125, 0.2, 0.2, 0.21, 0.21]
# # alpha_obs_low = [-0.1, -0.1, 0.045, 0.045, 0.13, 0.13, 0.10, 0.10]
# # plt.plot(mass_obs, alpha_obs_high, c='k', ls='dashed')
# # plt.plot(mass_obs, alpha_obs_low, c='k', ls='dashed')
# # plt.plot(mass_obs, alpha_obs, c='k')
# # # plt.plot(mass_obs, alpha_obs, 'k', color='chocolate', label='Luminosity-weighted stellar [Mg/Fe] obs')
# # # plt.fill_between(mass_obs, alpha_obs_low, alpha_obs_high, alpha=0.5, facecolor='chocolate', linewidth=4, antialiased=True)


# mass_obs_array = np.array([9,   10,    11,   12])
# alpha_obs_array_high = np.array([0.08,  0.205,  0.27, 0.32])
# alpha_obs_array = np.array([-0.01, 0.125, 0.2, 0.21])
# alpha_obs_array_low = np.array([-0.1,  0.045,   0.13, 0.10])
# fun_obs_alpha_mass_high = interpolate.interp1d(mass_obs_array, alpha_obs_array_high, kind='cubic')
# fun_obs_alpha_mass = interpolate.interp1d(mass_obs_array, alpha_obs_array, kind='cubic')
# fun_obs_alpha_mass_low = interpolate.interp1d(mass_obs_array, alpha_obs_array_low, kind='cubic')

# mass_new = np.arange(9, 12.01, 0.01)
# alpha_new_high = fun_obs_alpha_mass_high(mass_new)
# alpha_new = fun_obs_alpha_mass(mass_new)
# alpha_new_low = fun_obs_alpha_mass_low(mass_new)
# # plt.plot(mass_new, alpha_new, 'k', color='chocolate', label='Luminosity-weighted stellar [Mg/Fe] obs')
# # plt.fill_between(mass_new, alpha_new_low, alpha_new_high, alpha=0.5, facecolor='chocolate', linewidth=4, antialiased=True)
# plt.plot(mass_new, alpha_new_high, c='k', ls='dashed')
# plt.plot(mass_new, alpha_new, c='k')
# plt.plot(mass_new, alpha_new_low, c='k', ls='dashed')

########################
######## Fig 1 #########

plt.xlim(9, 12)
plt.legend()
plt.tight_layout()
# plt.savefig('ZX_MgFe_obs.pdf', dpi=250)




# import simulation data points ###

file = open('simulation_results_from_galaxy_evol/Metal_mass_relation.txt', 'r')
data = file.readlines()
file.close()

file = open('simulation_results_from_galaxy_evol/Metal_mass_relation_igimf.txt', 'r')
data_igimf = file.readlines()
file.close()

# plot for STF = 0.3

line = []
line_igimf = []
IMF = []
IMF_igimf = []
Log_SFR = []
Log_SFR_igimf = []
SFEN = []
SFEN_igimf = []
STF = []
STF_igimf = []
alive_stellar_mass = []
alive_stellar_mass_igimf = []
Dynamical_mass = []
Dynamical_mass_igimf = []
Mass_weighted_stellar_Mg_over_Fe = []
Mass_weighted_stellar_Mg_over_Fe_igimf = []
Mass_weighted_stellar_Z_over_X = []
Mass_weighted_stellar_Z_over_X_igimf = []
luminosity_weighted_stellar_Mg_over_Fe = []
luminosity_weighted_stellar_Mg_over_Fe_igimf = []
luminosity_weighted_stellar_Z_over_X = []
luminosity_weighted_stellar_Z_over_X_igimf = []
gas_Mg_over_Fe = []
gas_Mg_over_Fe_igimf = []
gas_Z_over_X = []
gas_Z_over_X_igimf = []
middle_Mg_over_Fe = []
middle_Mg_over_Fe_igimf = []
middle_Z_over_X = []
middle_Z_over_X_igimf = []
error_Mg_over_Fe = []
error_Mg_over_Fe_igimf = []
error_Z_over_X = []
error_Z_over_X_igimf = []
points = np.array([[0, 0, 0]])  # [Dynamical_mass, Mass_weighted_stellar_Mg_over_Fe, Mass_weighted_stellar_Z_over_X]
points_igimf = np.array([[0, 0, 0]])  # [Dynamical_mass, Mass_weighted_stellar_Mg_over_Fe, Mass_weighted_stellar_Z_over_X]
metal_values = np.array([[0]])  #
metal_values_igimf = np.array([[0]])  #
alpha_values = np.array([[0]])  #
alpha_values_igimf = np.array([[0]])  #


length_data = len(data)
i = 1
while i < length_data:
    line_i = [float(x) for x in data[i].split()]
    if (canonical_IMF_assumption == "Kroupa" and line_i[0] == 0) \
            or (canonical_IMF_assumption == "Salpeter" and line_i[0] == 2):
        if True: #line_i[2] == 100:
            line.append(line_i)
            IMF.append(line_i[0])
            Log_SFR.append(line_i[1])
            SFEN.append(line_i[2])
            STF.append(line_i[3])
            alive_stellar_mass.append(line_i[4])
            Dynamical_mass.append(line_i[5])
            Mass_weighted_stellar_Mg_over_Fe.append(line_i[6])
            Mass_weighted_stellar_Z_over_X.append(line_i[7])
            points = np.append(points, [[line_i[2], line_i[3], line_i[5]]], axis=0)
            metal_values = np.append(metal_values, [[line_i[7]]], axis=0)
            alpha_values = np.append(alpha_values, [[line_i[6]]], axis=0)
            # metal_values = np.append(metal_values, [[(line_i[7] + line_i[9]) / 2]], axis=0)
            # alpha_values = np.append(alpha_values, [[(line_i[6] + line_i[8]) / 2]], axis=0)
            if len(line_i) > 8:
                gas_Mg_over_Fe.append(line_i[8])
                gas_Z_over_X.append(line_i[9])
                middle_Mg_over_Fe.append((line_i[6] + line_i[8]) / 2)
                middle_Z_over_X.append((line_i[7] + line_i[9]) / 2)
                error_Mg_over_Fe.append(abs((line_i[6] - line_i[8]) / 2))
                error_Z_over_X.append(abs((line_i[7] - line_i[9]) / 2))
            else:
                gas_Mg_over_Fe.append(None)
                gas_Z_over_X.append(None)
            if len(line_i) > 10:
                luminosity_weighted_stellar_Mg_over_Fe.append(line_i[10])
                luminosity_weighted_stellar_Z_over_X.append(line_i[11])
            else:
                luminosity_weighted_stellar_Mg_over_Fe.append(None)
                luminosity_weighted_stellar_Z_over_X.append(None)
    (i) = (i+1)

length_data_igimf = len(data_igimf)
i = 1
while i < length_data_igimf:
    line_i = [float(x) for x in data_igimf[i].split()]
    if True:  # line_i[2] == 100:
        line_igimf.append(line_i)
        IMF_igimf.append(line_i[0])
        Log_SFR_igimf.append(line_i[1])
        SFEN_igimf.append(line_i[2])
        STF_igimf.append(line_i[3])
        alive_stellar_mass_igimf.append(line_i[4])
        Dynamical_mass_igimf.append(line_i[5])
        Mass_weighted_stellar_Mg_over_Fe_igimf.append(line_i[6])
        Mass_weighted_stellar_Z_over_X_igimf.append(line_i[7])
        points_igimf = np.append(points_igimf, [[line_i[2], line_i[3], line_i[5]]], axis=0)
        metal_values_igimf = np.append(metal_values_igimf, [[line_i[7]]], axis=0)
        alpha_values_igimf = np.append(alpha_values_igimf, [[line_i[6]]], axis=0)
        # metal_values_igimf = np.append(metal_values_igimf, [[(line_i[7] + line_i[9]) / 2]], axis=0)
        # alpha_values_igimf = np.append(alpha_values_igimf, [[(line_i[6] + line_i[8]) / 2]], axis=0)
        gas_Mg_over_Fe_igimf.append(line_i[8])
        gas_Z_over_X_igimf.append(line_i[9])
        middle_Mg_over_Fe_igimf.append((line_i[6] + line_i[8]) / 2)
        middle_Z_over_X_igimf.append((line_i[7] + line_i[9]) / 2)
        error_Mg_over_Fe_igimf.append(abs((line_i[6] - line_i[8]) / 2))
        error_Z_over_X_igimf.append(abs((line_i[7] - line_i[9]) / 2))
        if len(line_i) > 10:
            luminosity_weighted_stellar_Mg_over_Fe_igimf.append(line_i[10])
            luminosity_weighted_stellar_Z_over_X_igimf.append(line_i[11])
        else:
            luminosity_weighted_stellar_Mg_over_Fe_igimf.append(None)
            luminosity_weighted_stellar_Z_over_X_igimf.append(None)
            # points_igimf.append([line_i[5], line_i[6], line_i[7]])
    (i) = (i+1)
points = np.delete(points, 0, axis=0)
points_igimf = np.delete(points_igimf, 0, axis=0)
metal_values = np.delete(metal_values, 0, axis=0)
metal_values_igimf = np.delete(metal_values_igimf, 0, axis=0)
alpha_values = np.delete(alpha_values, 0, axis=0)
alpha_values_igimf = np.delete(alpha_values_igimf, 0, axis=0)


#################### setup interpretation ####################

grid_SFEN, grid_STF, grid_Dynamical_mass = np.mgrid[5:800:100j, 0:1:100j, 9:12:100j]
# the range of grid_Dynamical_mass is limited by the observation data
number_of_SFEN = 100
number_of_STF = 100
number_of_Dynamical_mass = 100

# grid_alpha0 = griddata(points, alpha_values, (grid_SFEN, grid_STF, grid_Dynamical_mass), method='nearest')
grid_alpha1 = griddata(points, alpha_values, (grid_SFEN, grid_STF, grid_Dynamical_mass), method='linear')

# grid_metal0 = griddata(points, metal_values, (grid_SFEN, grid_STF, grid_Dynamical_mass), method='nearest')
grid_metal1 = griddata(points, metal_values, (grid_SFEN, grid_STF, grid_Dynamical_mass), method='linear')

# grid_alpha_igimf0 = griddata(points_igimf, alpha_values_igimf, (grid_SFEN, grid_STF, grid_Dynamical_mass), method='nearest')
grid_alpha_igimf1 = griddata(points_igimf, alpha_values_igimf, (grid_SFEN, grid_STF, grid_Dynamical_mass), method='linear')

# grid_metal_igimf0 = griddata(points_igimf, metal_values_igimf, (grid_SFEN, grid_STF, grid_Dynamical_mass), method='nearest')
grid_metal_igimf1 = griddata(points_igimf, metal_values_igimf, (grid_SFEN, grid_STF, grid_Dynamical_mass), method='linear')

######### define functions that calculate likelihood #########

def calculate_likelihood(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha):
    if STF_j < len(grid_SFEN[SFEN_i]):
        SFEN_interpolate = grid_SFEN[SFEN_i][STF_j][Dynamical_mass_k]
        STF_interpolate = grid_STF[SFEN_i][STF_j][Dynamical_mass_k]
        Dynamical_mass_interpolate = grid_Dynamical_mass[SFEN_i][STF_j][Dynamical_mass_k]

        alpha_likelihood = calculate_likelihood_alpha(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_alpha)
        metal_likelihood = calculate_likelihood_metal(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal)

        # total_likelihood = metal_likelihood
        # total_likelihood = alpha_likelihood
        total_likelihood = alpha_likelihood * metal_likelihood
        return SFEN_interpolate, STF_interpolate, Dynamical_mass_interpolate, total_likelihood
    else:
        return 0, 0, 0, 0


def calculate_likelihood_alpha(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_alpha):
    Dynamical_mass_interpolate__ = grid_Dynamical_mass[SFEN_i][STF_j][Dynamical_mass_k]

    if IMF == "igimf":
        alpha_interpolate__ = grid_alpha_igimf1[SFEN_i][STF_j][Dynamical_mass_k]
    if IMF == "Kroupa":
        alpha_interpolate__ = grid_alpha1[SFEN_i][STF_j][Dynamical_mass_k]

    # alpha_interpolate_obs__ = fun_obs_alpha_mass(Dynamical_mass_interpolate__) + error_alpha
    # std_alpha_obs = abs(
    #     fun_obs_alpha_mass_high(Dynamical_mass_interpolate__) - fun_obs_alpha_mass_low(Dynamical_mass_interpolate__))/2

    alpha_interpolate_obs__ = -0.459+0.062*Dynamical_mass_interpolate__ + error_alpha
    std_alpha_obs = 0.1

    mismatch_alpha = abs(alpha_interpolate__ - alpha_interpolate_obs__)
    x_alpha = mismatch_alpha / std_alpha_obs
    alpha_likelihood = math.erfc(x_alpha/2**0.5)

    return alpha_likelihood



def calculate_likelihood_metal(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal):
    Dynamical_mass_interpolate = grid_Dynamical_mass[SFEN_i][STF_j][Dynamical_mass_k]

    if IMF == "igimf":
        metal_interpolate = grid_metal_igimf1[SFEN_i][STF_j][Dynamical_mass_k]
    if IMF == "Kroupa":
        metal_interpolate = grid_metal1[SFEN_i][STF_j][Dynamical_mass_k]

    metal_interpolate_obs = fun_obs_metal_mass(Dynamical_mass_interpolate) + error_metal

    std_metal_obs = abs(
        fun_obs_metal_mass_high(Dynamical_mass_interpolate) - fun_obs_metal_mass_low(Dynamical_mass_interpolate))/2
    # std_metal_obs = std_alpha_obs
    mismatch_metal = abs(metal_interpolate - metal_interpolate_obs)
    x_metal = mismatch_metal / std_metal_obs
    metal_likelihood = math.erfc(x_metal/2**0.5)
    return metal_likelihood


######### SFT--galaxy-mass relations given by Thomas05 and Recchi09 #########

def f1(t):
    return math.exp(3.67-0.37*t)
def f2(t):
    return math.exp(2.38-0.24*t)
t = [8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12]
y1 = [f1(8.5), f1(9), f1(9.5), f1(10), f1(10.5), f1(11), f1(11.5), f1(12)]
y2 = [f2(8.5), f2(9), f2(9.5), f2(10), f2(10.5), f2(11), f2(11.5), f2(12)]


######### assume typical systematic error in observation to be 0.3 dex (see Conroy et al. 2014) #########
######### and setup modifications on the observational abundances #########

error_metal = 0
error_alpha = 0
standard_deviation_systematic_metal = 0.3
standard_deviation_systematic_alpha = 0.3


IMF = "Kroupa"
error_metal_list = [-0.2, 0, 0.2]
error_alpha_list = [-0.2, 0, 0.2]


# IMF = "igimf"
# error_metal_list = [-0.2, 0, 0.2]
# error_alpha_list = [-0.2, 0, 0.2]


# IMF = "Kroupa"
# error_metal_list = [0]
# error_alpha_list = [0]


# IMF = "igimf"
# error_metal_list = [0.2]
# error_alpha_list = [0.3]






### plot data ###

# if IMF == "igimf":
#     plt.scatter(Dynamical_mass_igimf, Mass_weighted_stellar_Z_over_X_igimf, s=3, label='[Z/X] igimf')
#     plt.scatter(Dynamical_mass_igimf, Mass_weighted_stellar_Mg_over_Fe_igimf, s=3, label='[Mg/Fe] igimf')
#     # plt.errorbar(Dynamical_mass_igimf, middle_Z_over_X_igimf, yerr=error_Z_over_X_igimf, capsize=4, linestyle="None", fmt='o', label='[Z/X] igimf')
#     # plt.errorbar(Dynamical_mass_igimf, middle_Mg_over_Fe_igimf, yerr=error_Mg_over_Fe_igimf, capsize=4, linestyle="None", fmt='o', label='[Mg/Fe] igimf')
# elif IMF == "Kroupa":
#     plt.scatter(Dynamical_mass, Mass_weighted_stellar_Z_over_X, s=3, label='[Z/X] Kroupa-IMF')
#     plt.scatter(Dynamical_mass, Mass_weighted_stellar_Mg_over_Fe, s=3, label='[Mg/Fe] Kroupa-IMF')
#     # plt.errorbar(Dynamical_mass, middle_Z_over_X, yerr=error_Z_over_X, capsize=4, linestyle="None", fmt='o', label='[Z/X] Kroupa-IMF')
#     # plt.errorbar(Dynamical_mass, middle_Mg_over_Fe, yerr=error_Mg_over_Fe, capsize=4, linestyle="None", fmt='o', label='[Mg/Fe] Kroupa-IMF')





# plt.scatter(Dynamical_mass, Mass_weighted_stellar_Z_over_X, label='Mass_weighted_stellar_Z_over_X Kroupa-IMF')
# plt.scatter(Dynamical_mass, Mass_weighted_stellar_Mg_over_Fe, label='Mass_weighted_stellar_Mg_over_Fe Kroupa-IMF')
# # plt.scatter(Dynamical_mass, luminosity_weighted_stellar_Z_over_X, label='luminosity-weighted stellar [Z/X] Kroupa-IMF')
# # plt.scatter(Dynamical_mass, luminosity_weighted_stellar_Mg_over_Fe, label='luminosity-weighted stellar [Mg/Fe] Kroupa-IMF')
# plt.scatter(Dynamical_mass, gas_Z_over_X, label='gas [Z/X] Kroupa-IMF')
# plt.scatter(Dynamical_mass, gas_Mg_over_Fe, label='gas [Mg/Fe] Kroupa-IMF')
#
# plt.scatter(Dynamical_mass_igimf, Mass_weighted_stellar_Z_over_X_igimf, label='Mass-weighted stellar [Z/X] igimf')
# plt.scatter(Dynamical_mass_igimf, Mass_weighted_stellar_Mg_over_Fe_igimf, label='Mass-weighted stellar [Mg/Fe] igimf')
# # plt.scatter(Dynamical_mass_igimf, luminosity_weighted_stellar_Z_over_X_igimf, label='luminosity-weighted stellar [Z/X] igimf')
# # plt.scatter(Dynamical_mass_igimf, luminosity_weighted_stellar_Mg_over_Fe_igimf, label='luminosity-weighted stellar [Mg/Fe] igimf')
# plt.scatter(Dynamical_mass_igimf, gas_Z_over_X_igimf, label='gas [Z/X] igimf')
# plt.scatter(Dynamical_mass_igimf, gas_Mg_over_Fe_igimf, label='gas [Mg/Fe] igimf')




# # plot for STF = 0.2
#
# line = []
# line_igimf = []
# IMF = []
# IMF_igimf = []
# Log_SFR = []
# Log_SFR_igimf = []
# SFEN = []
# SFEN_igimf = []
# STF = []
# STF_igimf = []
# alive_stellar_mass = []
# alive_stellar_mass_igimf = []
# Dynamical_mass = []
# Dynamical_mass_igimf = []
# Mass_weighted_stellar_Mg_over_Fe = []
# Mass_weighted_stellar_Mg_over_Fe_igimf = []
# Mass_weighted_stellar_Z_over_X = []
# Mass_weighted_stellar_Z_over_X_igimf = []
# luminosity_weighted_stellar_Mg_over_Fe = []
# luminosity_weighted_stellar_Mg_over_Fe_igimf = []
# luminosity_weighted_stellar_Z_over_X = []
# luminosity_weighted_stellar_Z_over_X_igimf = []
# gas_Mg_over_Fe = []
# gas_Mg_over_Fe_igimf = []
# gas_Z_over_X = []
# gas_Z_over_X_igimf = []
# middle_Mg_over_Fe = []
# middle_Mg_over_Fe_igimf = []
# middle_Z_over_X = []
# middle_Z_over_X_igimf = []
# error_Mg_over_Fe = []
# error_Mg_over_Fe_igimf = []
# error_Z_over_X = []
# error_Z_over_X_igimf = []
# # points = []  # [Dynamical_mass, Mass_weighted_stellar_Mg_over_Fe, Mass_weighted_stellar_Z_over_X]
# # points_igimf = []  # [Dynamical_mass, Mass_weighted_stellar_Mg_over_Fe, Mass_weighted_stellar_Z_over_X]
# # values = []  # SFEN
#
# length_data = len(data)
# i = 1
# while i < length_data:
#     line_i = [float(x) for x in data[i].split()]
#     if True: #line_i[3] == 0.2 and line_i[1] > 0:
#         line.append(line_i)
#         IMF.append(line_i[0])
#         Log_SFR.append(line_i[1])
#         SFEN.append(line_i[2])
#         STF.append(line_i[3])
#         alive_stellar_mass.append(line_i[4])
#         Dynamical_mass.append(line_i[5])
#         Mass_weighted_stellar_Mg_over_Fe.append(line_i[6])
#         Mass_weighted_stellar_Z_over_X.append(line_i[7])
#         if len(line_i) > 8:
#             gas_Mg_over_Fe.append(line_i[8])
#             gas_Z_over_X.append(line_i[9])
#             middle_Mg_over_Fe.append((line_i[6] + line_i[8]) / 2)
#             middle_Z_over_X.append((line_i[7] + line_i[9]) / 2)
#             error_Mg_over_Fe.append(abs((line_i[6] - line_i[8]) / 2))
#             error_Z_over_X.append(abs((line_i[7] - line_i[9]) / 2))
#         else:
#             gas_Mg_over_Fe.append(None)
#             gas_Z_over_X.append(None)
#         if len(line_i) > 10:
#             luminosity_weighted_stellar_Mg_over_Fe.append(line_i[10])
#             luminosity_weighted_stellar_Z_over_X.append(line_i[11])
#         else:
#             luminosity_weighted_stellar_Mg_over_Fe.append(None)
#             luminosity_weighted_stellar_Z_over_X.append(None)
#     (i) = (i+1)
#
# length_data_igimf = len(data_igimf)
# i = 1
# while i < length_data_igimf:
#     line_i = [float(x) for x in data_igimf[i].split()]
#     if True: #line_i[3] == 0.2 and line_i[1] > 0:
#         line_igimf.append(line_i)
#         IMF_igimf.append(line_i[0])
#         Log_SFR_igimf.append(line_i[1])
#         SFEN_igimf.append(line_i[2])
#         STF_igimf.append(line_i[3])
#         alive_stellar_mass_igimf.append(line_i[4])
#         Dynamical_mass_igimf.append(line_i[5])
#         Mass_weighted_stellar_Mg_over_Fe_igimf.append(line_i[6])
#         Mass_weighted_stellar_Z_over_X_igimf.append(line_i[7])
#         gas_Mg_over_Fe_igimf.append(line_i[8])
#         gas_Z_over_X_igimf.append(line_i[9])
#         middle_Mg_over_Fe_igimf.append((line_i[6] + line_i[8]) / 2)
#         middle_Z_over_X_igimf.append((line_i[7] + line_i[9]) / 2)
#         error_Mg_over_Fe_igimf.append(abs((line_i[6] - line_i[8]) / 2))
#         error_Z_over_X_igimf.append(abs((line_i[7] - line_i[9]) / 2))
#         if len(line_i) > 10:
#             luminosity_weighted_stellar_Mg_over_Fe_igimf.append(line_i[10])
#             luminosity_weighted_stellar_Z_over_X_igimf.append(line_i[11])
#         else:
#             luminosity_weighted_stellar_Mg_over_Fe_igimf.append(None)
#             luminosity_weighted_stellar_Z_over_X_igimf.append(None)
#             # points_igimf.append([line_i[5], line_i[6], line_i[7]])
#     (i) = (i+1)
#
# # plt.scatter(Dynamical_mass, Mass_weighted_stellar_Z_over_X, label='Mass_weighted_stellar_Z_over_X Kroupa-IMF')
# # plt.scatter(Dynamical_mass, Mass_weighted_stellar_Mg_over_Fe, label='Mass_weighted_stellar_Mg_over_Fe Kroupa-IMF')
# # # plt.scatter(Dynamical_mass, luminosity_weighted_stellar_Z_over_X, label='luminosity-weighted stellar [Z/X] Kroupa-IMF')
# # # plt.scatter(Dynamical_mass, luminosity_weighted_stellar_Mg_over_Fe, label='luminosity-weighted stellar [Mg/Fe] Kroupa-IMF')
# # plt.scatter(Dynamical_mass, gas_Z_over_X, label='gas [Z/X] Kroupa-IMF')
# # plt.scatter(Dynamical_mass, gas_Mg_over_Fe, label='gas [Mg/Fe] Kroupa-IMF')
# #
# # plt.scatter(Dynamical_mass_igimf, Mass_weighted_stellar_Z_over_X_igimf, label='Mass-weighted stellar [Z/X] igimf')
# # plt.scatter(Dynamical_mass_igimf, Mass_weighted_stellar_Mg_over_Fe_igimf, label='Mass-weighted stellar [Mg/Fe] igimf')
# # # plt.scatter(Dynamical_mass_igimf, luminosity_weighted_stellar_Z_over_X_igimf, label='luminosity-weighted stellar [Z/X] igimf')
# # # plt.scatter(Dynamical_mass_igimf, luminosity_weighted_stellar_Mg_over_Fe_igimf, label='luminosity-weighted stellar [Mg/Fe] igimf')
# # plt.scatter(Dynamical_mass_igimf, gas_Z_over_X_igimf, label='gas [Z/X] igimf')
# # plt.scatter(Dynamical_mass_igimf, gas_Mg_over_Fe_igimf, label='gas [Mg/Fe] igimf')
#
# plt.errorbar(Dynamical_mass_igimf, middle_Mg_over_Fe_igimf, yerr=error_Mg_over_Fe_igimf, capsize=4, linestyle="None", fmt='o', label='[Mg/Fe] igimf')
# # plt.errorbar(Dynamical_mass, middle_Mg_over_Fe, yerr=error_Mg_over_Fe, capsize=4, linestyle="None", fmt='o', label='[Mg/Fe] Kroupa-IMF')
# plt.errorbar(Dynamical_mass_igimf, middle_Z_over_X_igimf, yerr=error_Z_over_X_igimf, capsize=4, linestyle="None", fmt='o', label='[Z/X] igimf')
# # plt.errorbar(Dynamical_mass, middle_Z_over_X, yerr=error_Z_over_X, capsize=4, linestyle="None", fmt='o', label='[Z/X] Kroupa-IMF')


#
# plt.xlabel(r'log$_{10}$(M$_{dyn}$ [M$_\odot$])')
# plt.ylabel(r'[Z/X] (blue) or [Mg/Fe] (orange)')
# plt.legend(prop={'size': 6}, loc='lower right')
# # plt.xlim(-1.2, 2.26)
# # plt.ylim(-1.05, 0.8)
# plt.tight_layout()
# # plt.savefig('Metal_mass_relation.png', dpi=250)




# # STF fixe for entire map:
# fig = plt.figure(2, figsize=(5, 5))
# SFEN_interpolate_list = []
# STF_interpolate_list = []
# Dynamical_mass_interpolate_list = []
# total_likelihood_list = []
#
# for SFEN_i in range(number_of_SFEN):
#     for Dynamical_mass_k in range(number_of_Dynamical_mass):
#         STF_j = 80
#         SFEN_interpolate, STF_interpolate, Dynamical_mass_interpolate, total_likelihood = \
#             calculate_likelihood(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)
#         SFEN_interpolate_list.append(SFEN_interpolate/100)
#         STF_interpolate_list.append(STF_interpolate)
#         Dynamical_mass_interpolate_list.append(Dynamical_mass_interpolate)
#         total_likelihood_list.append(total_likelihood)
# sc = plt.scatter(Dynamical_mass_interpolate_list, SFEN_interpolate_list, c=total_likelihood_list, s=30)
# plt.colorbar(sc)
# plt.plot(t, y1, c='r')
# plt.plot(t, y2, c='orange')


# # STF, the one with the best outcome likelihood, fixed for a given galaxy mass:
# fig = plt.figure(3, figsize=(5, 5))
# SFEN_interpolate_list = []
# STF_interpolate_list = []
# Dynamical_mass_interpolate_list = []
# total_likelihood_list = []
#
# for Dynamical_mass_k in range(number_of_Dynamical_mass):
#     total_likelihood_ = []
#     for STF_j in range(number_of_STF):
#         total_likelihood__ = []
#         for SFEN_i in range(number_of_SFEN):
#             SFEN_interpolate, STF_interpolate, Dynamical_mass_interpolate, total_likelihood = \
#                 calculate_likelihood(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)
#             if np.isnan(total_likelihood):
#                 total_likelihood__.append(0)
#             else:
#                 total_likelihood__.append(total_likelihood)
#         SFEN_i = total_likelihood__.index(max(total_likelihood__))
#         total_likelihood_.append(total_likelihood__[SFEN_i])
#     STF_j = total_likelihood_.index(max(total_likelihood_))
#     for SFEN_i in range(number_of_SFEN):
#         SFEN_interpolate, STF_interpolate, Dynamical_mass_interpolate, total_likelihood = \
#             calculate_likelihood(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)
#         SFEN_interpolate_list.append(SFEN_interpolate/100)
#         STF_interpolate_list.append(STF_interpolate) # STF_interpolate is the value while STF_j is the list index
#         Dynamical_mass_interpolate_list.append(Dynamical_mass_interpolate)
#         total_likelihood_list.append(total_likelihood)
# plt.scatter(Dynamical_mass_interpolate_list, SFEN_interpolate_list, c=total_likelihood_list, s=30)



# # STF is variable for all galaxy-mass and SFT, the plot shows the highest possible likelihood:
# fig = plt.figure(4, figsize=(5, 5))
# SFEN_interpolate_list = []
# STF_interpolate_list = []
# Dynamical_mass_interpolate_list = []
# total_likelihood_list = []
#
# for Dynamical_mass_k in range(number_of_Dynamical_mass):
#     for SFEN_i in range(number_of_SFEN):
#         total_likelihood__ = []
#         for STF_j in range(number_of_STF):
#             SFEN_interpolate, STF_interpolate, Dynamical_mass_interpolate, total_likelihood = \
#                 calculate_likelihood(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)
#             if np.isnan(total_likelihood):
#                 total_likelihood__.append(0)
#             else:
#                 total_likelihood__.append(total_likelihood)
#         STF_j = total_likelihood__.index(max(total_likelihood__))
#         SFEN_interpolate, STF_interpolate, Dynamical_mass_interpolate, total_likelihood = \
#             calculate_likelihood(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)
#         SFEN_interpolate_list.append(SFEN_interpolate/100)
#         STF_interpolate_list.append(STF_interpolate)
#         Dynamical_mass_interpolate_list.append(Dynamical_mass_interpolate)
#         total_likelihood_list.append(total_likelihood)
# plt.scatter(Dynamical_mass_interpolate_list, SFEN_interpolate_list, c=total_likelihood_list, s=30)
#



# Instead of fit both obseravtion simutanously,
# here first fit only Z to determine STF
# then fit Mg/Fe to determin SFT
# In this way, it mimics the Recchi or Thomas's work that only consider [Mg/Fe]
########################
###### Fig 2 & 6 #######

# fig = plt.figure(10, figsize=(4, 3.5))
# plt.xlabel(r'log$_{10}$(M$_{*}$ [M$_\odot$])')
# plt.ylabel(r't$_{\rm sf}$ [Gyr]')
# plt.xlim(9, 12)
# plt.ylim(0, 4)
# for error_metal, error_alpha in itertools.product(error_metal_list, error_alpha_list):
#     print(error_metal, error_alpha)
#     prior_likelihood = math.erfc(abs(error_metal) / standard_deviation_systematic_metal / 2 ** 0.5) * \
#                        math.erfc(abs(error_alpha) / standard_deviation_systematic_alpha / 2 ** 0.5)
#     SFEN_interpolate_list = []
#     STF_interpolate_list = []
#     Dynamical_mass_interpolate_list = []
#     total_likelihood_list = []
#     best_SFEN_list = []
#     dynamical_mass_list = []
#
#     for Dynamical_mass_k in range(number_of_Dynamical_mass):
#         total_likelihood__ = 0
#         print(Dynamical_mass_k+1, '%')
#         for STF_j__, SFEN_i__ in itertools.product(range(number_of_STF), range(number_of_SFEN)):
#             total_likelihood = calculate_likelihood_metal(SFEN_i__, STF_j__, Dynamical_mass_k, IMF, error_metal)
#             # total_likelihood = calculate_likelihood(SFEN_i__, STF_j__, Dynamical_mass_k, IMF, error_metal, error_alpha)[3]
#             if not np.isnan(total_likelihood):
#                 if total_likelihood > total_likelihood__:
#                     total_likelihood__ = total_likelihood
#                     STF_j = STF_j__
#         ### cut the plot ###
#         SFEN_test_1 = round(number_of_SFEN/10*8)
#         SFEN_test_2 = round(number_of_SFEN/10*2)
#         total_likelihood_test_1 = calculate_likelihood(SFEN_test_1, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)[3]
#         total_likelihood_test_2 = calculate_likelihood(SFEN_test_2, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)[3]
#         # total_likelihood_test_1 = True
#         # total_likelihood_test_2 = True
#         if not (np.isnan(total_likelihood_test_1) and np.isnan(total_likelihood_test_2)):
#             total_likelihood__ = 0
#             SFEN__ = 0
#             for SFEN_i in range(number_of_SFEN):
#                 SFEN_interpolate, STF_interpolate, Dynamical_mass_interpolate, total_likelihood = \
#                     calculate_likelihood(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)
#                 SFEN_interpolate_list.append(SFEN_interpolate / 100)
#                 STF_interpolate_list.append(STF_interpolate)
#                 Dynamical_mass_interpolate_list.append(Dynamical_mass_interpolate)
#                 total_likelihood_list.append(total_likelihood * prior_likelihood)
#                 if not np.isnan(total_likelihood):
#                     if total_likelihood > total_likelihood__:
#                         total_likelihood__ = total_likelihood
#                         SFEN__ = SFEN_interpolate
#             if SFEN__ == 0:
#                 best_SFEN_list.append(None)
#             else:
#                 best_SFEN_list.append(SFEN__ / 100)
#             dynamical_mass_list.append(Dynamical_mass_interpolate)
#     # print("Plotting...")
#     plt.scatter(Dynamical_mass_interpolate_list, SFEN_interpolate_list, c=total_likelihood_list, s=10)
#     plt.plot(dynamical_mass_list, best_SFEN_list, c='r',linewidth=0.8)
#     sc = plt.scatter([0, 0], [0, 0], c=[0, 1], s=1)
#     clb = plt.colorbar(sc)
#     clb.set_label('likelihood')
#     plt.title('[Z/X]_lit+({}); [Mg/Fe]_lit+({})'.format(error_metal, error_alpha), fontsize=8)
# plt.plot(t, y1, c='k', ls='dashed', label='Thomas05')
# mass_thomas_list = np.arange(9, 12.1, 0.1)
# SFT_thomas_list = []
# for i in mass_thomas_list:
#     alpha_interpolate_obs = fun_obs_alpha_mass(i)
#     SFT_thomas = 10**((1/5-alpha_interpolate_obs)*6)
#     SFT_thomas_list.append(SFT_thomas)
# plt.plot(mass_thomas_list, SFT_thomas_list, c='k', ls='dotted', label='Thomas05-new')
# plt.plot(t, y2, c='k', ls='-.', label='Recchi09')
# plt.legend(prop={'size': 7}, loc='upper right')
# plt.tight_layout()
# if IMF == "Kroupa":
#     plt.savefig('likelihood_Kroupa.pdf', dpi=250)
# if IMF == "igimf":
#     plt.savefig('likelihood_igimf.pdf', dpi=250)

########################
##### Fig A1 & A2 ######

# fig = plt.figure(10, figsize=(8, 6))
# ax = fig.add_subplot(111)
# ax.spines['top'].set_color('none')
# ax.spines['bottom'].set_color('none')
# ax.spines['left'].set_color('none')
# ax.spines['right'].set_color('none')
# ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
# ax.set_xlabel(r'log$_{10}$(M$_{*}$ [M$_\odot$])')
# ax.set_ylabel(r't$_{\rm sf}$ [Gyr]')
# # sc = ax.scatter([0, 0], [0, 0], c=[0, 1], s=1)
# # clb = plt.colorbar(sc, shrink=1, aspect=40)
# # clb.set_label('likelihood')
# plot_number = 0
# for error_metal, error_alpha in itertools.product(error_metal_list, error_alpha_list):
#     print(error_metal, error_alpha)
#     prior_likelihood = math.erfc(abs(error_metal) / standard_deviation_systematic_metal / 2 ** 0.5) * \
#                        math.erfc(abs(error_alpha) / standard_deviation_systematic_alpha / 2 ** 0.5)
#     plot_number = plot_number + 1
#     ax1 = fig.add_subplot(3, 3, plot_number)
#     SFEN_interpolate_list = []
#     STF_interpolate_list = []
#     Dynamical_mass_interpolate_list = []
#     total_likelihood_list = []
#     best_SFEN_list = []
#     dynamical_mass_list = []
#
#     for Dynamical_mass_k in range(number_of_Dynamical_mass):
#         total_likelihood__ = 0
#         # print(Dynamical_mass_k+1, '%')
#         for STF_j__, SFEN_i__ in itertools.product(range(number_of_STF), range(number_of_SFEN)):
#             total_likelihood = calculate_likelihood_metal(SFEN_i__, STF_j__, Dynamical_mass_k, IMF, error_metal)
#             # total_likelihood = calculate_likelihood(SFEN_i__, STF_j__, Dynamical_mass_k, IMF, error_metal, error_alpha)[3]
#             if not np.isnan(total_likelihood):
#                 if total_likelihood > total_likelihood__:
#                     total_likelihood__ = total_likelihood
#                     STF_j = STF_j__
#         SFEN_test_1 = round(number_of_SFEN/10*8)
#         SFEN_test_2 = round(number_of_SFEN/10*2)
#         total_likelihood_test_1 = calculate_likelihood(SFEN_test_1, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)[3]
#         total_likelihood_test_2 = calculate_likelihood(SFEN_test_2, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)[3]
#         # total_likelihood_test_1 = True
#         # total_likelihood_test_2 = True
#         if not (np.isnan(total_likelihood_test_1) and np.isnan(total_likelihood_test_2)):
#             total_likelihood__ = 0
#             SFEN__ = 0
#             for SFEN_i in range(number_of_SFEN):
#                 SFEN_interpolate, STF_interpolate, Dynamical_mass_interpolate, total_likelihood = \
#                     calculate_likelihood(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)
#                 SFEN_interpolate_list.append(SFEN_interpolate / 100)
#                 STF_interpolate_list.append(STF_interpolate)
#                 Dynamical_mass_interpolate_list.append(Dynamical_mass_interpolate)
#                 total_likelihood_list.append(total_likelihood * prior_likelihood)
#                 if not np.isnan(total_likelihood):
#                     if total_likelihood > total_likelihood__:
#                         total_likelihood__ = total_likelihood
#                         SFEN__ = SFEN_interpolate
#             if SFEN__ == 0:
#                 best_SFEN_list.append(None)
#             else:
#                 best_SFEN_list.append(SFEN__ / 100)
#             dynamical_mass_list.append(Dynamical_mass_interpolate)
#     # print("Plotting...")
#     ax1.scatter(Dynamical_mass_interpolate_list, SFEN_interpolate_list, c=total_likelihood_list, s=10)
#     ax1.plot(dynamical_mass_list, best_SFEN_list, c='r',linewidth=0.8)
#     ax1.plot(t, y1, c='k', ls='dashed')
#     ax1.plot(t, y2, c='k', ls='-.')
#     ax1.set_xlim(8.4, 12.1)
#     ax1.set_ylim(0, 8.1)
#     sc = ax1.scatter([0, 0], [0, 0], c=[0, 1], s=1)
#     # if plot_number == 3 or plot_number == 6 or plot_number == 9:
#     #     clb = plt.colorbar(sc)
#     #     clb.set_label('likelihood')
#     ax1.set_title('[Z/X]_lit+({}); [Mg/Fe]_lit+({})'.format(error_metal, error_alpha), fontsize=8)
# plt.tight_layout()
# if IMF == "Kroupa":
#     plt.savefig('grid_Kroupa.pdf', dpi=250)
# if IMF == "igimf":
#     plt.savefig('grid_igimf.pdf', dpi=250)

#################################
##### Fig 3, 4, 5; 7, 8, 9 ######
plot_number = 0
flag0 = 0
flag1 = 0
flag2 = 0
for error_metal, error_alpha in itertools.product(error_metal_list, error_alpha_list):
    print(error_metal, error_alpha)
    prior_likelihood = math.erfc(abs(error_metal) / standard_deviation_systematic_metal / 2 ** 0.5) * \
                       math.erfc(abs(error_alpha) / standard_deviation_systematic_alpha / 2 ** 0.5)
    plot_number = plot_number + 1
    # ax1 = fig.add_subplot(3, 3, plot_number)
    SFEN_interpolate_list = []
    STF_interpolate_list = []
    Dynamical_mass_interpolate_list = []
    total_likelihood_list = []
    best_SFEN_list = []
    dynamical_mass_list = []
    alpha_interpolate_obs_list = []

    for Dynamical_mass_k in range(number_of_Dynamical_mass):
        total_likelihood__ = 0
        # print(Dynamical_mass_k+1, '%')
        for STF_j__, SFEN_i__ in itertools.product(range(number_of_STF), range(number_of_SFEN)):
            total_likelihood = calculate_likelihood_metal(SFEN_i__, STF_j__, Dynamical_mass_k, IMF, error_metal)
            # total_likelihood = calculate_likelihood(SFEN_i__, STF_j__, Dynamical_mass_k, IMF, error_metal, error_alpha)[3]
            if not np.isnan(total_likelihood):
                if total_likelihood > total_likelihood__:
                    total_likelihood__ = total_likelihood
                    STF_j = STF_j__
        SFEN_test_1 = round(number_of_SFEN/10*8)
        SFEN_test_2 = round(number_of_SFEN/10*2)
        total_likelihood_test_1 = calculate_likelihood(SFEN_test_1, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)[3]
        total_likelihood_test_2 = calculate_likelihood(SFEN_test_2, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)[3]
        # total_likelihood_test_1 = True
        # total_likelihood_test_2 = True
        if not (np.isnan(total_likelihood_test_1) and np.isnan(total_likelihood_test_2)):
            total_likelihood__ = 0
            SFEN__ = 0
            for SFEN_i in range(number_of_SFEN):
                SFEN_interpolate, STF_interpolate, Dynamical_mass_interpolate, total_likelihood = \
                    calculate_likelihood(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)
                SFEN_interpolate_list.append(SFEN_interpolate / 100)
                STF_interpolate_list.append(STF_interpolate)
                Dynamical_mass_interpolate_list.append(Dynamical_mass_interpolate)
                total_likelihood_list.append(total_likelihood * prior_likelihood)
                if not np.isnan(total_likelihood):
                    if total_likelihood > total_likelihood__:
                        total_likelihood__ = total_likelihood
                        SFEN__ = SFEN_interpolate
            if SFEN__ == 0:
                best_SFEN_list.append(None)
            else:
                # alpha_interpolate_obs = fun_obs_alpha_mass(Dynamical_mass_interpolate) + error_alpha
                alpha_interpolate_obs = -0.459+0.062*Dynamical_mass_interpolate + error_alpha
                alpha_interpolate_obs_list.append(alpha_interpolate_obs)
                best_SFEN_list.append(SFEN__ / 100)
            dynamical_mass_list.append(Dynamical_mass_interpolate)
    fig = plt.figure(1001, figsize=(4, 3.5))
    if error_metal == error_metal_list[0] and flag0 == 0:
        flag0 = 1
        plt.plot(Dynamical_mass_interpolate_list, STF_interpolate_list, label='[Z/X] {}'.format(error_metal), c='k', ls='dashed')
    if error_metal == error_metal_list[1] and flag1 == 0:
        flag1 = 1
        plt.plot(Dynamical_mass_interpolate_list, STF_interpolate_list, label='[Z/X] {}'.format(error_metal), c='k')
    if error_metal == error_metal_list[2] and flag2 == 0:
        flag2 = 1
        plt.plot(Dynamical_mass_interpolate_list, STF_interpolate_list, label='[Z/X] {}'.format(error_metal), c='k', ls='dotted')

    fig = plt.figure(1002, figsize=(4, 3.5))
    if error_metal == error_metal_list[0] and error_alpha == error_alpha_list[0]:
        plt.plot(dynamical_mass_list, best_SFEN_list, c='r', ls='dashed')
    if error_metal == error_metal_list[0] and error_alpha == error_alpha_list[1]:
        plt.plot(dynamical_mass_list, best_SFEN_list, c='g', ls='dashed')
    if error_metal == error_metal_list[0] and error_alpha == error_alpha_list[2]:
        plt.plot(dynamical_mass_list, best_SFEN_list, c='b', ls='dashed')
    if error_metal == error_metal_list[1] and error_alpha == error_alpha_list[0]:
        plt.plot(dynamical_mass_list, best_SFEN_list, c='r')
    if error_metal == error_metal_list[1] and error_alpha == error_alpha_list[1]:
        plt.plot(dynamical_mass_list, best_SFEN_list, c='g')
    if error_metal == error_metal_list[1] and error_alpha == error_alpha_list[2]:
        plt.plot(dynamical_mass_list, best_SFEN_list, c='b')
    if error_metal == error_metal_list[2] and error_alpha == error_alpha_list[0]:
        plt.plot(dynamical_mass_list, best_SFEN_list, c='r', ls='dotted')
    if error_metal == error_metal_list[2] and error_alpha == error_alpha_list[1]:
        plt.plot(dynamical_mass_list, best_SFEN_list, c='g', ls='dotted')
    if error_metal == error_metal_list[2] and error_alpha == error_alpha_list[2]:
        plt.plot(dynamical_mass_list, best_SFEN_list, c='b', ls='dotted')

    fig = plt.figure(1003, figsize=(4, 3.5))
    if error_metal == error_metal_list[0] and error_alpha == error_alpha_list[0]:
        plt.plot(best_SFEN_list, alpha_interpolate_obs_list, c='r', ls='dashed')
    if error_metal == error_metal_list[0] and error_alpha == error_alpha_list[1]:
        plt.plot(best_SFEN_list, alpha_interpolate_obs_list, c='g', ls='dashed')
    if error_metal == error_metal_list[0] and error_alpha == error_alpha_list[2]:
        plt.plot(best_SFEN_list, alpha_interpolate_obs_list, c='b', ls='dashed')
    if error_metal == error_metal_list[1] and error_alpha == error_alpha_list[0]:
        plt.plot(best_SFEN_list, alpha_interpolate_obs_list, c='r')
    if error_metal == error_metal_list[1] and error_alpha == error_alpha_list[1]:
        plt.plot(best_SFEN_list, alpha_interpolate_obs_list, c='g')
    if error_metal == error_metal_list[1] and error_alpha == error_alpha_list[2]:
        plt.plot(best_SFEN_list, alpha_interpolate_obs_list, c='b')
    if error_metal == error_metal_list[2] and error_alpha == error_alpha_list[0]:
        plt.plot(best_SFEN_list, alpha_interpolate_obs_list, c='r', ls='dotted')
    if error_metal == error_metal_list[2] and error_alpha == error_alpha_list[1]:
        plt.plot(best_SFEN_list, alpha_interpolate_obs_list, c='g', ls='dotted')
    if error_metal == error_metal_list[2] and error_alpha == error_alpha_list[2]:
        plt.plot(best_SFEN_list, alpha_interpolate_obs_list, c='b', ls='dotted')
    # plt.scatter(best_SFEN_list, alpha_interpolate_obs_list, c='r', s=1)

fig = plt.figure(1001, figsize=(4, 3.5))
plt.xlabel(r'log$_{10}$(M$_{*}$ [M$_\odot$])')
plt.ylabel(r'f$_{\rm st}$')
plt.xlim(9, 12)
# plt.ylim(0, 1)
# plt.plot([], [], label='[Z/X], [Mg/Fe]:')
plt.legend(prop={'size': 7}, loc='best')
plt.tight_layout()
if IMF == "Kroupa":
    plt.savefig('best_STF_Kroupa.pdf', dpi=250)
if IMF == "igimf":
    plt.savefig('best_STF_igimf.pdf', dpi=250)


fig = plt.figure(1002, figsize=(4, 3.5))
plt.xlabel(r'log$_{10}$(M$_{*}$ [M$_\odot$])')
plt.ylabel(r't$_{\rm sf}$ [Gyr]')
plt.xlim(9, 13)
plt.ylim(0, 4)
plt.plot([], [], label='[Z/X] {}'.format(error_metal_list[0]), c='k', ls='dashed')
plt.plot([], [], label='[Z/X] {}'.format(error_metal_list[1]), c='k')
plt.plot([], [], label='[Z/X] {}'.format(error_metal_list[2]), c='k', ls='dotted')
plt.plot([], [], label='[Mg/Fe] {}'.format(error_alpha_list[0]), c='r')
plt.plot([], [], label='[Mg/Fe] {}'.format(error_alpha_list[1]), c='g')
plt.plot([], [], label='[Mg/Fe] {}'.format(error_alpha_list[2]), c='b')
plt.plot(t, y1, c='k', ls='dashed', label='Thomas05')
plt.plot(t, y2, c='k', ls='-.', label='Recchi09')
mass_thomas_list = np.arange(9, 12.1, 0.1)
SFT_thomas_list = []
for i in mass_thomas_list:
    # alpha_interpolate_obs = fun_obs_alpha_mass(i)
    alpha_interpolate_obs = -0.459+0.062*i
    SFT_thomas = 10**((1/5-alpha_interpolate_obs)*6)
    SFT_thomas_list.append(SFT_thomas)
plt.plot(mass_thomas_list, SFT_thomas_list, c='k', ls='dotted', label='Thomas05-new')
# plt.plot([], [], label='[Z/X], [Mg/Fe]:')
plt.legend(prop={'size': 7}, loc='upper right')
plt.tight_layout()
if IMF == "Kroupa":
    plt.savefig('best_SFEN_Kroupa.pdf', dpi=250)
if IMF == "igimf":
    plt.savefig('best_SFEN_igimf.pdf', dpi=250)


fig = plt.figure(1003, figsize=(4, 3.5))
plt.xlabel(r't$_{\rm sf}$ [Gyr]')
plt.ylabel('[Mg/Fe]')
plt.xlim(0, 4)
plt.ylim(0, 0.65)
plt.plot([], [], label='[Z/X] {}'.format(error_metal_list[0]), c='k', ls='dashed')
plt.plot([], [], label='[Z/X] {}'.format(error_metal_list[1]), c='k')
plt.plot([], [], label='[Z/X] {}'.format(error_metal_list[2]), c='k', ls='dotted')
plt.plot([], [], label='[Mg/Fe] {}'.format(error_alpha_list[0]), c='r')
plt.plot([], [], label='[Mg/Fe] {}'.format(error_alpha_list[1]), c='g')
plt.plot([], [], label='[Mg/Fe] {}'.format(error_alpha_list[2]), c='b')
SFT_thomas = np.arange(0.01, 4.1, 0.1)
alpha_thomas = []
for i in SFT_thomas:
    alpha_thomas.append(1/5-1/6*math.log(i, 10))
plt.plot(SFT_thomas, alpha_thomas, c='k', ls='dashed', label='Thomas05')
plt.legend(prop={'size': 7}, loc='best')
plt.tight_layout()
if IMF == "Kroupa":
    plt.savefig('alpha_SFEN_Kroupa.pdf', dpi=250)
if IMF == "igimf":
    plt.savefig('alpha_SFEN_igimf.pdf', dpi=250)




# # renormalize the color map for each galaxy-mass bin
# fig = plt.figure(6, figsize=(5, 5))
# for Dynamical_mass_k in range(number_of_Dynamical_mass):
#     SFEN_interpolate_list = []
#     STF_interpolate_list = []
#     Dynamical_mass_interpolate_list = []
#     total_likelihood_list = []
#
#     total_likelihood_ = []
#     for STF_j in range(number_of_STF):
#         total_likelihood__ = []
#         for SFEN_i in range(number_of_SFEN):
#             SFEN_interpolate, STF_interpolate, Dynamical_mass_interpolate, total_likelihood = \
#                 calculate_likelihood(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)
#             if np.isnan(total_likelihood):
#                 total_likelihood__.append(0)
#             else:
#                 total_likelihood__.append(total_likelihood)
#         SFEN_i = total_likelihood__.index(max(total_likelihood__))
#         total_likelihood_.append(total_likelihood__[SFEN_i])
#     STF_j = total_likelihood_.index(max(total_likelihood_))
#     for SFEN_i in range(number_of_SFEN):
#         SFEN_interpolate, STF_interpolate, Dynamical_mass_interpolate, total_likelihood = \
#             calculate_likelihood(SFEN_i, STF_j, Dynamical_mass_k, IMF, error_metal, error_alpha)
#         SFEN_interpolate_list.append(SFEN_interpolate/100)
#         STF_interpolate_list.append(STF_interpolate)
#         Dynamical_mass_interpolate_list.append(Dynamical_mass_interpolate)
#         total_likelihood_list.append(total_likelihood)
#     plt.scatter(Dynamical_mass_interpolate_list, SFEN_interpolate_list, c=total_likelihood_list, s=30)





plt.tight_layout()
plt.show()
