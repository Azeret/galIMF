# A python3 code
# This is a single-zone closed-box galaxy chemical evolution module.
# It is coupled with a variable galaxy-wide IMF that depends on the galactic property at the time of star formation.
# The stellar population forms at every 10 Myr (the shortest time step) over 10 Gyr;
# with each stellar population a different galaxy-wide IMF calculated using the IGIMF theory (the galimf.py model).
# The parameters assumed for the simulation are specified at the end of this file or imported from other files,
# i.e., element_weight_table.py, element_abundances_solar.py, element_abundances_primordial.py.


import numpy as np
import time
import math
import importlib
from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import gc
import sys
import warnings
import os
import errno
warnings.filterwarnings("ignore")

sys.path.insert(0, 'Generated_IGIMFs')
sys.path.insert(0, 'IMFs')
sys.path.insert(0, 'yield_tables')
import element_weight_table, element_abundances_solar, element_abundances_primordial, stellar_luminosity
from IMFs import Kroupa_IMF, diet_Salpeter_IMF
from yield_tables import SNIa_yield


def galaxy_evol(imf='igimf', STF=0.5, SFEN=1, Z_0=0.000000134, solar_mass_component='Anders1989_mass',
                str_yield_table='portinari98',
                IMF_name='Kroupa', steller_mass_upper_bound=150,
                time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
                SFH_model='provided', SFE=0.05,
                SNIa_ON=True, SNIa_yield_table='Thielemann1993', solar_abu_table='Anders1989',
                high_time_resolution=None, plot_show=None, plot_save=None, outflow=None, check_igimf=False):
    start_time = time.time()

    ######################
    # If imf='igimf', the model will use variable IMF, imf='Kroupa' will use Kroupa IMF
    # A 1 in SFH.txt stand for SFR = 1 [solar mass/year] in a 10 Myr epoch.
    # STF is the total stellar mass/total gas mass in 13Gyr, which determines the initial gas mass. See Yan et al. 2019
    # Z_0 is the initial metallicity
    ######################
    global igimf_mass_function, mass_grid_table, mass_grid_table2, Mfinal_table, Mmetal_table, M_element_table
    global time_axis, gas_Z_over_X_list, total_star_formed, \
        O_over_H_list, Mg_over_H_list, C_over_H_list, N_over_H_list, Ca_over_H_list, Fe_over_H_list, \
        S_over_H_list, Si_over_H_list, Ne_over_H_list, \
        S_over_H_list, Si_over_H_list, Ne_over_H_list, X_list, Y_list, Z_list, \
        ejected_O_mass_till_this_time_tot_list, ejected_O_mass_till_this_time_SNII_list, ejected_O_mass_till_this_time_SNIa_list, \
        ejected_Mg_mass_till_this_time_tot_list, ejected_Mg_mass_till_this_time_SNII_list, ejected_Mg_mass_till_this_time_SNIa_list, \
        ejected_Fe_mass_till_this_time_tot_list, ejected_Fe_mass_till_this_time_SNII_list, ejected_Fe_mass_till_this_time_SNIa_list, \
        ejected_S_mass_till_this_time_tot_list, ejected_S_mass_till_this_time_SNII_list, ejected_S_mass_till_this_time_SNIa_list, \
        ejected_Si_mass_till_this_time_tot_list, ejected_Si_mass_till_this_time_SNII_list, ejected_Si_mass_till_this_time_SNIa_list, \
        ejected_Ne_mass_till_this_time_tot_list, ejected_Ne_mass_till_this_time_SNII_list, ejected_Ne_mass_till_this_time_SNIa_list, \
        ejected_Ca_mass_till_this_time_tot_list, ejected_Ca_mass_till_this_time_SNII_list, ejected_Ca_mass_till_this_time_SNIa_list, \
        Mg_over_Fe_list, C_over_Fe_list, N_over_O_list, Ca_over_Fe_list, O_over_Fe_list, S_over_Fe_list, Si_over_Fe_list, Ne_over_Fe_list, \
        stellar_O_over_H_list, stellar_Mg_over_H_list, stellar_C_over_H_list, stellar_N_over_H_list, \
        stellar_Ca_over_H_list, stellar_Fe_over_H_list, stellar_Si_over_H_list, stellar_S_over_H_list, stellar_Ne_over_H_list, \
        stellar_X_list, stellar_Y_list, stellar_Z_list, \
        stellar_Mg_over_Fe_list, stellar_C_over_Fe_list, stellar_N_over_O_list, stellar_Ca_over_Fe_list, \
        stellar_S_over_Fe_list, stellar_Si_over_Fe_list, stellar_Ne_over_Fe_list, \
        stellar_O_over_Fe_list, stellar_Z_over_X_list, stellar_Z_over_H_list, \
        stellar_O_over_H_list_luminosity_weighted, stellar_Mg_over_H_list_luminosity_weighted, \
        stellar_C_over_H_list_luminosity_weighted, stellar_N_over_H_list_luminosity_weighted, \
        stellar_Ca_over_H_list_luminosity_weighted, stellar_Ne_over_H_list_luminosity_weighted, stellar_Si_over_H_list_luminosity_weighted, stellar_S_over_H_list_luminosity_weighted, \
        stellar_X_list_luminosity_weighted, stellar_Y_list_luminosity_weighted, stellar_Z_list_luminosity_weighted, \
        stellar_Fe_over_H_list_luminosity_weighted, stellar_Mg_over_Fe_list_luminosity_weighted, \
        stellar_C_over_Fe_list_luminosity_weighted, stellar_N_over_O_list_luminosity_weighted, \
        stellar_Ca_over_Fe_list_luminosity_weighted, stellar_O_over_Fe_list_luminosity_weighted, \
        stellar_S_over_Fe_list_luminosity_weighted, stellar_Si_over_Fe_list_luminosity_weighted, stellar_Ne_over_Fe_list_luminosity_weighted, \
        stellar_Z_over_X_list_luminosity_weighted, stellar_Z_over_H_list_luminosity_weighted, \
        remnant_mass_list, total_gas_mass_list, stellar_mass_list, \
        ejected_gas_mass_list, ejected_metal_mass_list, \
        ejected_gas_Mg_over_Fe_list, instant_ejected_gas_Mg_over_Fe_list, expansion_factor_list, \
        expansion_factor_instantaneous_list, expansion_factor_adiabat_list
    global SNIa_energy_release_list, SNIa_number_list, SNIa_number_per_century, \
        SNII_energy_release_list, SNII_number_list, SNII_number_per_century, \
        total_energy_release_list, SN_number_per_century, total_gas_kinetic_energy_list, original_gas_mass  # , binding_energy_list
    global BH_mass_list, NS_mass_list, WD_mass_list, all_sf_imf, all_sfr
    global times_calculate_igimf, instantaneous_recycling, primary_He_mass_fraction
    global X_solar, Y_solar, Z_solar, log_binding_energy_initial

    instantaneous_recycling = False
    times_calculate_igimf = 0
    ###################
    ### preparation ###
    ###################

    # Warning flags:
    Warning_ejected_gas_mass_of_this_epoch = False
    Warning_WD_mass_till_this_time = False
    Warning_galaxy_mass_ejected_gas_mass = False

    # get all avaliable metallicity from stellar evolution table
    (Z_table_list, Z_table_list_2, Z_table_list_3) = function_get_avaliable_Z(str_yield_table)

    # read in SFH
    SFH_input = np.loadtxt('SFH.txt')
    length_list_SFH_input = len(SFH_input)

    i = 0
    total_SF = 0
    while i < length_list_SFH_input:
        total_SF += SFH_input[i]
        (i) = (i + 1)

    # Star Trasnformation fraction (STF)
    total_star_formed = 10 ** 7 * total_SF
    original_gas_mass = total_star_formed / STF  # in solar mass unit
    print("original_gas_mass =", math.log(original_gas_mass, 10))
    ratio_gas_over_DM_radii = 0.3
    log_binding_energy_initial = round(
        math.log(6.674 * 1.989 ** 2 / 3.086 / 242, 10) + 40 +
        math.log(0.5 * original_gas_mass ** 2 + ratio_gas_over_DM_radii * (
        1 + 1.37 * ratio_gas_over_DM_radii) / 2 / 3.1415926 * 3 * 10 ** 6 * original_gas_mass, 10), 3)
    print("log_binding_energy_initial =", log_binding_energy_initial)
    # Create the time steps (x axis) for final output
    time_axis = [10**6]
    time_resolution = time_resolution_in_Myr * 10 ** 5 * 10
    for i in range(10 ** 9, 15 * 10 ** 9, time_resolution * 1000):
        time_axis += [i]
    if high_time_resolution==True:
        for i in range(10 ** 7, 10 ** 8, time_resolution * 10):
            time_axis += [i]
        for i in range(10 ** 8, 10 ** 9, time_resolution * 100):
            time_axis += [i]
    else:
        plot_at_age = [1 * 10 ** 7, 2 * 10 ** 7, 3 * 10 ** 7, 4 * 10 ** 7, 5 * 10 ** 7,
                       6 * 10 ** 7, 7 * 10 ** 7, 8 * 10 ** 7, 9 * 10 ** 7,
                       1 * 10 ** 8, 2 * 10 ** 8, 3 * 10 ** 8, 4 * 10 ** 8, 5 * 10 ** 8,
                       6 * 10 ** 8, 7 * 10 ** 8, 8 * 10 ** 8, 9 * 10 ** 8,
                       1 * 10 ** 9, 101 * 10 ** 7, 102 * 10 ** 7, 103 * 10 ** 7,
                       104 * 10 ** 7, 105 * 10 ** 7, 106 * 10 ** 7, 107 * 10 ** 7, 108 * 10 ** 7, 11 * 10 ** 8,
                       12 * 10 ** 8, 14 * 10 ** 8, 16 * 10 ** 8, 18 * 10 ** 8,
                       2 * 10 ** 9, 23 * 10 ** 8, 26 * 10 ** 8, 29 * 10 ** 8,
                       32 * 10 ** 8, 35 * 10 ** 8, 38 * 10 ** 8, 41 * 10 ** 8, 45 * 10 ** 8,
                       5 * 10 ** 9, 6 * 10 ** 9,
                       7 * 10 ** 9, 8 * 10 ** 9, 9 * 10 ** 9, 10 * 10 ** 9, 11 * 10 ** 9]
        # plot_at_age = [1 * 10 ** 8, 1 * 10 ** 9, 10.8 * 10 ** 9]
        time_axis += plot_at_age
        for i in range(10 ** 9, 15 * 10 ** 9, time_resolution * 1000):
            time_axis += [i]

    # consider also all star formation event happend times
    # where the time resolution should be temporarily increased.
    time_axis_for_SFH_input = []
    time_axis_for_SFH_input_D = []
    i = 0
    while i < length_list_SFH_input:
        if SFH_input[i] > 0:
            if high_time_resolution == True:
                add_time = 1
                while add_time < 70:
                    add_time_step = round(10**(add_time/20)) * 10 ** 7
                    add_time_step += i * 10 ** 7 + add_time_step
                    if add_time_step < 14 * 10**9:
                        time_axis_for_SFH_input += [add_time_step]
                    (add_time) = (add_time+1)

                # time_axis_for_SFH_input += [i * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 1 * 10 ** 4]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 5 * 10 ** 4]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 7 * 10 ** 4]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 8 * 10 ** 4]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 9 * 10 ** 4]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 1 * 10 ** 5]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 2 * 10 ** 5]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 5 * 10 ** 5]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 1 * 10 ** 6]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 2 * 10 ** 6]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 5 * 10 ** 6]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 1 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 2 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 3 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 4 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 5 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 6 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 7 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 8 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 9 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 10 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 11 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 12 * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 2 * 10 ** 8]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 5 * 10 ** 8]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 1 * 10 ** 9]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 2 * 10 ** 9]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 5 * 10 ** 9]
            else:
                time_axis_for_SFH_input += [i * 10 ** 7]
                # time_axis_for_SFH_input += [i * 10 ** 7 + 9 * 10 ** 6]
        (i) = (i + 1)

    # the final time axis is the sorted combination of the two
    time_axis = sorted(list(set(time_axis + time_axis_for_SFH_input)))
    # print("\nSimulation results will be give at galactic age [yr] =\n", time_axis)
    length_list_time_step = len(time_axis)

    print('time_axis:', length_list_time_step, time_axis)
    ###################
    ###  main loop  ###
    ###################
    S_F_R_of_this_epoch_list = []
    S_F_F = 1
    # define an array save SF event informations that will be used in every latter time steps
    all_sf_imf = []
    all_sfr = []
    epoch_info = []  # This array saves the properties,
    # i.e., [S_F_R_of_this_epoch, M_tot_of_this_epoch, igimf_mass_function, igimf_normalization],
    # of all the stellar populations formed at different time step (the so-called star formation epoch),
    # Thus that at any later given time (aging), the effects (yield) of these populations can be easily computed.
    BH_mass_list = []
    NS_mass_list = []
    WD_mass_list = []
    gas_Z_over_X_list = []
    O_over_H_list = []
    ejected_O_mass_till_this_time_tot_list = []
    ejected_O_mass_till_this_time_SNII_list = []
    ejected_O_mass_till_this_time_SNIa_list = []
    ejected_Mg_mass_till_this_time_tot_list = []
    ejected_Mg_mass_till_this_time_SNII_list = []
    ejected_Mg_mass_till_this_time_SNIa_list = []
    ejected_Fe_mass_till_this_time_tot_list = []
    ejected_Fe_mass_till_this_time_SNII_list = []
    ejected_Fe_mass_till_this_time_SNIa_list = []
    ejected_Ca_mass_till_this_time_tot_list = []
    ejected_Ca_mass_till_this_time_SNII_list = []
    ejected_Ca_mass_till_this_time_SNIa_list = []
    ejected_S_mass_till_this_time_tot_list = []
    ejected_S_mass_till_this_time_SNII_list = []
    ejected_S_mass_till_this_time_SNIa_list = []
    ejected_Si_mass_till_this_time_tot_list = []
    ejected_Si_mass_till_this_time_SNII_list = []
    ejected_Si_mass_till_this_time_SNIa_list = []
    ejected_Ne_mass_till_this_time_tot_list = []
    ejected_Ne_mass_till_this_time_SNII_list = []
    ejected_Ne_mass_till_this_time_SNIa_list = []
    X_list = []
    Y_list = []
    Z_list = []
    Mg_over_H_list = []
    C_over_H_list = []
    N_over_H_list = []
    Ca_over_H_list = []
    Ne_over_H_list = []
    Si_over_H_list = []
    S_over_H_list = []
    Fe_over_H_list = []
    Mg_over_Fe_list = []
    C_over_Fe_list = []
    N_over_O_list = []
    Ca_over_Fe_list = []
    Ne_over_Fe_list = []
    Si_over_Fe_list = []
    S_over_Fe_list = []
    O_over_Fe_list = []
    stellar_O_over_H_list = []
    stellar_Mg_over_H_list = []
    stellar_C_over_H_list = []
    stellar_N_over_H_list = []
    stellar_Ca_over_H_list = []
    stellar_Ne_over_H_list = []
    stellar_Si_over_H_list = []
    stellar_S_over_H_list = []
    stellar_Fe_over_H_list = []
    stellar_X_list = []
    stellar_Y_list = []
    stellar_Z_list = []
    stellar_Mg_over_Fe_list = []
    stellar_C_over_Fe_list = []
    stellar_N_over_O_list = []
    stellar_Ca_over_Fe_list = []
    stellar_Ne_over_Fe_list = []
    stellar_Si_over_Fe_list = []
    stellar_S_over_Fe_list = []
    stellar_O_over_Fe_list = []
    stellar_Z_over_X_list = []
    stellar_Z_over_H_list = []
    stellar_O_over_H_list_luminosity_weighted = []
    stellar_Mg_over_H_list_luminosity_weighted = []
    stellar_C_over_H_list_luminosity_weighted = []
    stellar_N_over_H_list_luminosity_weighted = []
    stellar_Ca_over_H_list_luminosity_weighted = []
    stellar_Ne_over_H_list_luminosity_weighted = []
    stellar_Si_over_H_list_luminosity_weighted = []
    stellar_S_over_H_list_luminosity_weighted = []
    stellar_X_list_luminosity_weighted = []
    stellar_Y_list_luminosity_weighted = []
    stellar_Z_list_luminosity_weighted = []
    stellar_Fe_over_H_list_luminosity_weighted = []
    stellar_Mg_over_Fe_list_luminosity_weighted = []
    stellar_C_over_Fe_list_luminosity_weighted = []
    stellar_N_over_O_list_luminosity_weighted = []
    stellar_Ca_over_Fe_list_luminosity_weighted = []
    stellar_Ne_over_Fe_list_luminosity_weighted = []
    stellar_Si_over_Fe_list_luminosity_weighted = []
    stellar_S_over_Fe_list_luminosity_weighted = []
    stellar_O_over_Fe_list_luminosity_weighted = []
    stellar_Z_over_X_list_luminosity_weighted = []
    stellar_Z_over_H_list_luminosity_weighted = []
    remnant_mass_list = []
    total_gas_mass_list = []
    ejected_gas_mass_list = []
    ejected_gas_Mg_over_Fe_list = []
    instant_ejected_gas_Mg_over_Fe_list = []
    ejected_metal_mass_list = []
    expansion_factor_instantaneous_list = []
    expansion_factor_adiabat_list = []
    expansion_factor_list = []
    stellar_mass_list = []
    total_energy_release_list = []
    SN_number_per_century = []
    # binding_energy_list = []
    total_gas_kinetic_energy_list = []
    SNIa_energy_release_list = []
    SNIa_number_list = []
    SNIa_number_per_century = []
    SNII_energy_release_list = []
    SNII_number_list = []
    SNII_number_per_century = []
    Z_solar = element_abundances_solar.function_solar_element_abundances(solar_mass_component, 'Metal')
    Y_solar = element_abundances_solar.function_solar_element_abundances(solar_mass_component, 'He')
    X_solar = element_abundances_solar.function_solar_element_abundances(solar_mass_component, 'H')
    primary_H_mass_fraction = element_abundances_primordial.function_element_mass_primary_fraction(
        solar_abu_table, "H", Z_0, Z_solar)
    primary_He_mass_fraction = element_abundances_primordial.function_element_mass_primary_fraction(
        solar_abu_table, "He", Z_0, Z_solar)
    Z_over_X = math.log(Z_0 / primary_H_mass_fraction, 10) - math.log(Z_solar / X_solar, 10)
    # do calculation for each time start from time 0
    time_step = 0
    gc_collect_check = 1
    # do calculation for each time to the end time
    while time_step < length_list_time_step:
        # get time
        this_time = time_axis[time_step]
        # calculated the array index (line number in SFH.txt) this_time has reached
        epoch_index_limit = (this_time + 1) / 10 ** 7
        if epoch_index_limit > length_list_SFH_input:
            epoch_index_limit = length_list_SFH_input
        last_time_age = 0
        age_of_this_epoch = 0
        number_in_SNIa_boundary = 0
        # get masses and metallicity at this time (values are calculated by the end of last time step)
        # initialize values
        total_energy_release = 0
        SNIa_energy_release = 0
        SNIa_number_from_all_epoch = 0
        SNII_energy_release = 0
        SNII_number = 0
        if time_step == 0:
            eject_H_mass = 0
            eject_C_mass = 0
            eject_N_mass = 0
            eject_O_mass = 0
            eject_Mg_mass = 0
            eject_Ca_mass = 0
            eject_Ne_mass = 0
            eject_Si_mass = 0
            eject_S_mass = 0
            eject_Fe_mass = 0
            eject_metal_mass = 0

            total_gas_mass_at_this_time = original_gas_mass
            ejected_gas_mass_at_this_time = 0
            ejected_metal_mass_at_last_time = 0

            M_tot_up_to_last_time = 0
            M_tot_up_to_this_time = 0
            stellar_mass_at_last_time = 0
            stellar_mass_at_this_time = 0
            stellar_luminosity_at_this_time = 0

            ejected_gas_mass_till_last_time = 0
            ejected_metal_mass_till_last_time = 0
            ejected_H_mass_till_last_time = 0
            ejected_He_mass_till_last_time = 0
            ejected_C_mass_till_last_time = 0
            ejected_N_mass_till_last_time = 0
            ejected_O_mass_till_last_time = 0
            ejected_Mg_mass_till_last_time = 0
            ejected_Ca_mass_till_last_time = 0
            ejected_Ne_mass_till_last_time = 0
            ejected_Si_mass_till_last_time = 0
            ejected_S_mass_till_last_time = 0
            ejected_Fe_mass_till_last_time = 0

            ejected_gas_mass_till_this_time = 0
            ejected_metal_mass_till_this_time = 0
            ejected_H_mass_till_this_time = 0
            ejected_He_mass_till_this_time = 0
            ejected_C_mass_till_this_time = 0
            ejected_N_mass_till_this_time = 0
            ejected_O_mass_till_this_time = 0
            ejected_Mg_mass_till_this_time = 0
            ejected_Ca_mass_till_this_time = 0
            ejected_Ne_mass_till_this_time = 0
            ejected_Si_mass_till_this_time = 0
            ejected_S_mass_till_this_time = 0
            ejected_Fe_mass_till_this_time = 0
            BH_mass_till_this_time = 0
            NS_mass_till_this_time = 0
            WD_mass_till_this_time = 0
            remnant_mass_at_this_time = 0

            # Fe_H_mass_ratio_at_last_time = 0 #################################
            Z_gas_this_time_step = Z_0
            total_metal_mass_at_this_time = total_gas_mass_at_this_time * Z_gas_this_time_step
            total_H_mass_at_this_time = 0
            total_He_mass_at_this_time = 0
            total_C_mass_at_this_time = 0
            total_N_mass_at_this_time = 0
            total_O_mass_at_this_time = 0
            total_Mg_mass_at_this_time = 0
            total_Ca_mass_at_this_time = 0
            total_Ne_mass_at_this_time = 0
            total_Si_mass_at_this_time = 0
            total_S_mass_at_this_time = 0
            total_Fe_mass_at_this_time = 0

            total_H_mass_at_last_time = original_gas_mass * primary_H_mass_fraction
            H_weight = element_weight_table.function_element_weight("H")
            total_He_mass_at_last_time = original_gas_mass * primary_He_mass_fraction
            total_C_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(
                solar_abu_table, "C", Z_0, Z_solar)
            total_N_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(
                solar_abu_table, "N", Z_0, Z_solar)
            total_O_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(
                solar_abu_table, "O", Z_0, Z_solar)
            total_Mg_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(
                solar_abu_table, "Mg", Z_0, Z_solar)
            total_Ca_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(
                solar_abu_table, "Ca", Z_0, Z_solar)
            total_Ne_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(
                solar_abu_table, "Ne", Z_0, Z_solar)
            total_Si_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(
                solar_abu_table, "Si", Z_0, Z_solar)
            total_S_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(
                solar_abu_table, "S", Z_0, Z_solar)
            total_Fe_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(
                solar_abu_table, "Fe", Z_0, Z_solar)
            total_metal_mass_in_gas_at_last_time = original_gas_mass * Z_0
            total_gas_mass_at_last_time = original_gas_mass

            stellar_metal_mass_at_this_time = 0
            stellar_H_mass_at_this_time = 0
            stellar_He_mass_at_this_time = 0
            stellar_C_mass_at_this_time = 0
            stellar_N_mass_at_this_time = 0
            stellar_O_mass_at_this_time = 0
            stellar_Mg_mass_at_this_time = 0
            stellar_Ca_mass_at_this_time = 0
            stellar_Fe_mass_at_this_time = 0
            stellar_Ne_mass_at_this_time = 0
            stellar_Si_mass_at_this_time = 0
            stellar_S_mass_at_this_time = 0

            stellar_metal_luminosity_at_this_time = 0
            stellar_H_luminosity_at_this_time = 0
            stellar_He_luminosity_at_this_time = 0
            stellar_C_luminosity_at_this_time = 0
            stellar_N_luminosity_at_this_time = 0
            stellar_O_luminosity_at_this_time = 0
            stellar_Mg_luminosity_at_this_time = 0
            stellar_Ca_luminosity_at_this_time = 0
            stellar_Fe_luminosity_at_this_time = 0
            stellar_Ne_luminosity_at_this_time = 0
            stellar_Si_luminosity_at_this_time = 0
            stellar_S_luminosity_at_this_time = 0

            metal_mass_fraction_in_gas = [Z_gas_this_time_step, primary_H_mass_fraction, primary_He_mass_fraction,
                                          total_C_mass_at_last_time / original_gas_mass,
                                          total_N_mass_at_last_time / original_gas_mass,
                                          total_O_mass_at_last_time / original_gas_mass,
                                          total_Mg_mass_at_last_time / original_gas_mass,
                                          total_Ca_mass_at_last_time / original_gas_mass,
                                          total_Fe_mass_at_last_time / original_gas_mass,
                                          total_Ne_mass_at_last_time / original_gas_mass,
                                          total_Si_mass_at_last_time / original_gas_mass,
                                          total_S_mass_at_last_time / original_gas_mass]
                                          # H He C N O Mg Ca Fe Ne Si S
        else:
            total_gas_mass_at_last_time = total_gas_mass_at_this_time
            # total_gas_mass_at_this_time is set in below
            ejected_gas_mass_at_this_time = 0
            total_metal_mass_in_gas_at_last_time = total_metal_mass_at_this_time
            total_metal_mass_at_this_time = 0
            total_H_mass_at_last_time = total_H_mass_at_this_time
            total_H_mass_at_this_time = 0
            total_He_mass_at_last_time = total_He_mass_at_this_time
            total_He_mass_at_this_time = 0
            total_C_mass_at_last_time = total_C_mass_at_this_time
            total_C_mass_at_this_time = 0
            total_N_mass_at_last_time = total_N_mass_at_this_time
            total_N_mass_at_this_time = 0
            total_O_mass_at_last_time = total_O_mass_at_this_time
            total_O_mass_at_this_time = 0
            total_Mg_mass_at_last_time = total_Mg_mass_at_this_time
            total_Mg_mass_at_this_time = 0
            total_Ca_mass_at_last_time = total_Ca_mass_at_this_time
            total_Ca_mass_at_this_time = 0
            total_Si_mass_at_last_time = total_Si_mass_at_this_time
            total_Si_mass_at_this_time = 0
            total_S_mass_at_last_time = total_S_mass_at_this_time
            total_S_mass_at_this_time = 0
            total_Ne_mass_at_last_time = total_Ne_mass_at_this_time
            total_Ne_mass_at_this_time = 0
            total_Fe_mass_at_last_time = total_Fe_mass_at_this_time
            total_Fe_mass_at_this_time = 0
            M_tot_up_to_last_time = M_tot_up_to_this_time
            M_tot_up_to_this_time = 0
            stellar_mass_at_last_time = stellar_mass_at_this_time
            stellar_mass_at_this_time = 0
            stellar_luminosity_at_this_time = 0
            BH_mass_till_this_time = 0
            NS_mass_till_this_time = 0
            WD_mass_till_this_time = 0
            remnant_mass_at_this_time = 0
            ejected_gas_mass_till_last_time = ejected_gas_mass_till_this_time
            ejected_metal_mass_till_last_time = ejected_metal_mass_till_this_time
            ejected_H_mass_till_last_time = ejected_H_mass_till_this_time
            ejected_He_mass_till_last_time = ejected_He_mass_till_this_time
            ejected_C_mass_till_last_time = ejected_C_mass_till_this_time
            ejected_N_mass_till_last_time = ejected_N_mass_till_this_time
            ejected_O_mass_till_last_time = ejected_O_mass_till_this_time
            ejected_Mg_mass_till_last_time = ejected_Mg_mass_till_this_time
            ejected_Ca_mass_till_last_time = ejected_Ca_mass_till_this_time
            ejected_Ne_mass_till_last_time = ejected_Ne_mass_till_this_time
            ejected_Si_mass_till_last_time = ejected_Si_mass_till_this_time
            ejected_S_mass_till_last_time = ejected_S_mass_till_this_time
            ejected_Fe_mass_till_last_time = ejected_Fe_mass_till_this_time
            ejected_gas_mass_till_this_time = 0
            ejected_metal_mass_till_this_time = 0
            ejected_H_mass_till_this_time = 0
            ejected_He_mass_till_this_time = 0
            ejected_C_mass_till_this_time = 0
            ejected_N_mass_till_this_time = 0
            ejected_O_mass_till_this_time = 0
            ejected_Mg_mass_till_this_time = 0
            ejected_Ca_mass_till_this_time = 0
            ejected_Ne_mass_till_this_time = 0
            ejected_Si_mass_till_this_time = 0
            ejected_S_mass_till_this_time = 0
            ejected_Fe_mass_till_this_time = 0
            ejected_metal_mass_at_last_time = ejected_metal_mass_at_this_time
            # Fe_H_mass_ratio_at_last_time = Fe_H_mass_ratio_at_this_time
            Z_gas_this_time_step = total_metal_mass_in_gas_at_last_time / total_gas_mass_at_last_time
            metal_mass_fraction_in_gas = [Z_gas_this_time_step,
                                          total_H_mass_at_last_time / total_gas_mass_at_last_time,
                                          total_He_mass_at_last_time / total_gas_mass_at_last_time,
                                          total_C_mass_at_last_time / total_gas_mass_at_last_time,
                                          total_N_mass_at_last_time / total_gas_mass_at_last_time,
                                          total_O_mass_at_last_time / total_gas_mass_at_last_time,
                                          total_Mg_mass_at_last_time / total_gas_mass_at_last_time,
                                          total_Ca_mass_at_last_time / total_gas_mass_at_last_time,
                                          total_Fe_mass_at_last_time / total_gas_mass_at_last_time,
                                          total_Ne_mass_at_last_time / total_gas_mass_at_last_time,
                                          total_Si_mass_at_last_time / total_gas_mass_at_last_time,
                                          total_S_mass_at_last_time / total_gas_mass_at_last_time]
            stellar_metal_mass_at_this_time = 0
            stellar_H_mass_at_this_time = 0
            stellar_He_mass_at_this_time = 0
            stellar_C_mass_at_this_time = 0
            stellar_N_mass_at_this_time = 0
            stellar_O_mass_at_this_time = 0
            stellar_Mg_mass_at_this_time = 0
            stellar_Ca_mass_at_this_time = 0
            stellar_Ne_mass_at_this_time = 0
            stellar_Si_mass_at_this_time = 0
            stellar_S_mass_at_this_time = 0
            stellar_Fe_mass_at_this_time = 0

            stellar_metal_luminosity_at_this_time = 0
            stellar_H_luminosity_at_this_time = 0
            stellar_He_luminosity_at_this_time = 0
            stellar_C_luminosity_at_this_time = 0
            stellar_N_luminosity_at_this_time = 0
            stellar_O_luminosity_at_this_time = 0
            stellar_Mg_luminosity_at_this_time = 0
            stellar_Ca_luminosity_at_this_time = 0
            stellar_Fe_luminosity_at_this_time = 0
            stellar_Ne_luminosity_at_this_time = 0
            stellar_Si_luminosity_at_this_time = 0
            stellar_S_luminosity_at_this_time = 0
        # add up metals contributed by SSP from each SF epoch
        # consider only the SF event (epoch) that had happend
        Fe_production_SNII = 0
        Mg_production_SNII = 0
        Ca_production_SNII = 0
        Ne_production_SNII = 0
        Si_production_SNII = 0
        S_production_SNII = 0
        O_production_SNII = 0

        epoch_index = 0
        while epoch_index < epoch_index_limit:
            # get age
            age_of_this_epoch = this_time - epoch_index * 10 ** 7
            # get SFR, M_tot, igimf, integrated igimf, stellar lifetime and stellar remnant mass for this metallicity
            # check if the info of this epoch has been recorded in previous time steps...
            if epoch_index == len(epoch_info):  # if not:
                # SFR
                if SFH_model == 'provided':
                    # This model apply the SFH specified by the SFH.txt
                    S_F_R_of_this_epoch = SFH_input[epoch_index]
                elif SFH_model == 'gas_mass_dependent':
                    # In this model, the SFR is determined by the current gas mass
                    # if the current time is shorter than SFEN * 10^7 yr.
                    S_F_R_of_this_epoch = total_gas_mass_at_this_time * SFE / 10 ** 7
                    if SFH_input[epoch_index] == 0 or S_F_R_of_this_epoch < 3.5*1e-6 or epoch_index>99:
                        S_F_R_of_this_epoch = 0
                    # print(epoch_index, '*10 Myr  SFR:', S_F_R_of_this_epoch)
                    print(S_F_R_of_this_epoch)
                else:
                    print("Wrong input parameter for 'SFH_model'.")

                # M_tot
                # if total_gas_mass_at_last_time > 10**12:
                #     M_tot_of_this_epoch = max((min(((total_gas_mass_at_last_time - 10 * stellar_mass_at_last_time) / 5), 10**12)), 0)
                # else:
                #     M_tot_of_this_epoch = 0
                M_tot_of_this_epoch = S_F_R_of_this_epoch * 10 ** 7

                # if S_F_F == 1:
                #     S_F_R_of_this_epoch = total_gas_mass_at_last_time**(0.99) * 3.97 * 10**(-10)  # Pflamm-Altenburg & Kroupa 2009
                #     S_F_R_of_this_epoch_list += [S_F_R_of_this_epoch]
                #     if S_F_R_of_this_epoch < S_F_R_of_this_epoch_list[0] * 0.8:
                #         S_F_F = 0
                # else:
                #     S_F_R_of_this_epoch = 0
                #
                #
                # print(S_F_R_of_this_epoch)
                # M_tot_of_this_epoch = S_F_R_of_this_epoch * 10 ** 7
                # if S_F_R_of_this_epoch > 0:
                #     if high_time_resolution == True:
                #         time_axis_for_SFH_input_D += [epoch_index * 10 ** 7]
                #         # time_axis_for_SFH_input_D += [epoch_index * 10 ** 7 + 5 * 10 ** 5]
                #         # time_axis_for_SFH_input_D += [epoch_index * 10 ** 7 + 1 * 10 ** 6]
                #         # time_axis_for_SFH_input_D += [epoch_index * 10 ** 7 + 2 * 10 ** 6]
                #         # time_axis_for_SFH_input_D += [epoch_index * 10 ** 7 + 4 * 10 ** 6]
                #         # time_axis_for_SFH_input_D += [epoch_index * 10 ** 7 + 6 * 10 ** 6]
                #         # time_axis_for_SFH_input_D += [epoch_index * 10 ** 7 + 8 * 10 ** 6]
                #         time_axis_for_SFH_input_D += [epoch_index * 10 ** 7 + 10 * 10 ** 6]
                #     else:
                #         time_axis_for_SFH_input_D += [epoch_index * 10 ** 7]
                #     time_axis = sorted(list(set(time_axis + time_axis_for_SFH_input_D)))
                #     length_list_time_step = len(time_axis)


                if S_F_R_of_this_epoch > 0:
                    # Total mass normalized IGIMF and unnormalized other IMFs
                    if imf == 'igimf':
                        igimf_of_this_epoch = function_get_igimf_for_this_epoch(S_F_R_of_this_epoch, Z_over_X,
                                                                                this_time, epoch_index,
                                                                                check_igimf)  # Fe_over_H_number_ratio)
                    elif imf == 'Kroupa':
                        igimf_of_this_epoch = Kroupa_IMF
                    elif imf == 'Salpeter':
                        from IMFs import Salpeter_IMF
                        igimf_of_this_epoch = Salpeter_IMF
                    elif imf == 'diet_Salpeter':
                        igimf_of_this_epoch = diet_Salpeter_IMF
                    elif imf == 'given':
                        from IMFs import given_IMF
                        igimf_of_this_epoch = given_IMF
                    igimf = igimf_of_this_epoch

                    #
                    def igimf_xi_function(mass):
                        return igimf_of_this_epoch.custom_imf(mass, this_time)

                    def igimf_mass_function(mass):
                        return igimf_of_this_epoch.custom_imf(mass, this_time) * mass

                    def igimf_luminous_function(mass):
                        return igimf_of_this_epoch.custom_imf(mass, this_time) * \
                               stellar_luminosity.stellar_luminosity_function(mass)

                    # integrated igimf_mass_function from 0.08 to steller_mass_upper_bound
                    if imf == 'diet_Salpeter':
                        integrate_igimf_mass = quad(igimf_mass_function, 0.1, 100, limit=50)[0]
                    else:
                        integrate_igimf_mass = quad(igimf_mass_function, 0.08, steller_mass_upper_bound, limit=50)[0]
                    # as the integration of the IGIMF always has a small (at least for low SFRs) computational error,
                    # it need to be fixed by mutiplying a calibration factor which is close to 1:
                    mass_calibration_factor = M_tot_of_this_epoch / integrate_igimf_mass
                    # print("mass_calibration_factor:", mass_calibration_factor)  # the calibration factor is about 1%

                    # integrate_igimf_mass_l = quad(igimf_mass_function, 0.08, 3, limit=40)[0]
                    # integrate_igimf_mass_h = quad(igimf_mass_function, 8, steller_mass_upper_bound, limit=50)[0]
                    # integrate_igimf_mass_m = quad(igimf_mass_function, 3, 8, limit=40)[0]
                    # print("high mass star mass ratio:", integrate_igimf_mass_h/integrate_igimf_mass)
                    # print("middle mass star mass ratio:", integrate_igimf_mass_m/integrate_igimf_mass)
                    # print("Low mass star mass ratio:", integrate_igimf_mass_l/integrate_igimf_mass)
                    # integrate_igimf_number = quad(igimf_xi_function, 0.08, steller_mass_upper_bound, limit=50)[0]
                    # integrate_igimf_number_l = quad(igimf_xi_function, 0.08, 3, limit=40)[0]
                    # integrate_igimf_number_h = quad(igimf_xi_function, 8, steller_mass_upper_bound, limit=50)[0]
                    # integrate_igimf_number_m = quad(igimf_xi_function, 3, 8, limit=40)[0]
                    # print("high mass star number ratio:", integrate_igimf_number_h/integrate_igimf_number)
                    # print("middle mass star number ratio:", integrate_igimf_number_m/integrate_igimf_number)
                    # print("Low mass star number ratio:", integrate_igimf_number_l/integrate_igimf_number)

                    # Choose the closest metallicity
                    Z_select_in_table = function_select_metal(Z_gas_this_time_step, Z_table_list)
                    # Z_select_in_table = ('in/out', Z_select__low, Z_gas_this_time_step, Z_select__high)
                    Z_select_in_table_2 = function_select_metal(Z_gas_this_time_step, Z_table_list_2)
                    if str_yield_table != "portinari98":
                        Z_select_in_table_3 = function_select_metal(Z_gas_this_time_step, Z_table_list_3)
                    else:
                        Z_select_in_table_3 = None
                    # read in interpolated stellar lifetime table
                    (mass_1, mass, lifetime_table) = function_read_lifetime(str_yield_table, Z_select_in_table)
                    # read in interpolated stellar final mass
                    (mass_12, Mfinal_table) = function_read_Mfinal(str_yield_table, Z_select_in_table)
                    # read in interpolated stellar ejected metal mass
                    (mass_2, mass2, Mmetal_table) = function_read_Mmetal(str_yield_table, Z_select_in_table_2,
                                                                         Z_select_in_table_3)
                    # read in interpolated stellar ejected elements mass
                    MH_table = function_read_M_element("H", str_yield_table, Z_select_in_table_2, Z_select_in_table_3)
                    MHe_table = function_read_M_element("He", str_yield_table, Z_select_in_table_2, Z_select_in_table_3)
                    MC_table = function_read_M_element("C", str_yield_table, Z_select_in_table_2, Z_select_in_table_3)
                    MN_table = function_read_M_element("N", str_yield_table, Z_select_in_table_2, Z_select_in_table_3)
                    MO_table = function_read_M_element("O", str_yield_table, Z_select_in_table_2, Z_select_in_table_3)
                    MMg_table = function_read_M_element("Mg", str_yield_table, Z_select_in_table_2, Z_select_in_table_3)
                    MNe_table = function_read_M_element("Ne", str_yield_table, Z_select_in_table_2, Z_select_in_table_3)
                    MSi_table = function_read_M_element("Si", str_yield_table, Z_select_in_table_2, Z_select_in_table_3)
                    MS_table = function_read_M_element("S", str_yield_table, Z_select_in_table_2, Z_select_in_table_3)
                    MCa_table = function_read_M_element("Ca", str_yield_table, Z_select_in_table_2, Z_select_in_table_3)
                    MFe_table = function_read_M_element("Fe", str_yield_table, Z_select_in_table_2, Z_select_in_table_3)
                    M_element_table = [MH_table, MHe_table, MC_table, MN_table, MO_table, MMg_table, MNe_table,
                                       MSi_table, MS_table, MCa_table, MFe_table]

                    # check if the in put lifetime and final mass table used the same mass grid
                    # if mass_1 != mass_12:
                    #     print('Error! Stellar lifetime and final mass input data do not match.\n'
                    #           'Check the table file: yield_tables/rearranged___/setllar_final_mass_from_portinari98/portinari98_Z={}.txt\n'
                    #           'and table file: yield_tables/rearranged___/setllar_lifetime_from_portinari98/portinari98_Z={}.txt'.format(
                    #                                                                            Z_select_in_table,
                    #                                                                            Z_select_in_table))
                    # else:
                    #     mass_grid_table = mass
                    #     mass_grid_table2 = mass2
                    mass_grid_table = mass
                    mass_grid_table2 = mass2

                    last_time_age = age_of_this_epoch
                    number_in_SNIa_boundary = mass_calibration_factor * quad(igimf_xi_function, 3, 8, limit=50)[
                        0]  # see function_number_SNIa below
                    if imf == 'diet_Salpeter' or imf == 'Salpeter':
                        number_all = mass_calibration_factor * quad(igimf_xi_function, 0.1, 100, limit=50)[0]  # see function_number_SNIa below
                    else:
                        number_all = mass_calibration_factor * quad(igimf_xi_function, 0.08, steller_mass_upper_bound, limit=50)[0]  # see function_number_SNIa below
                    # number_low = quad(igimf_xi_function, 0.08, 2, limit=40)[0]  # see function_number_SNIa below
                    # number_up = quad(igimf_xi_function, 8, steller_mass_upper_bound, limit=50)[0]  # see function_number_SNIa below
                    # print("up", number_up/number_all)

                    SNIa_number_prob = number_in_SNIa_boundary ** 2 / number_all / M_tot_of_this_epoch

                    # SNIa_number_prob = number_in_SNIa_boundary**2 / number_all * 10**2 * 0.61
                    # number_in_SNIa_boundary = SNIa_number_prob
                    # SNIa_number_prob = number_in_SNIa_boundary / integrate_igimf_mass
                    # print("SNIa SNIa_number_prob:", SNIa_number_prob)
                    # print("total star number", number_all)
                    # print("low", number_low/number_all)

                    age_of_this_epoch_at_end = (length_list_SFH_input - epoch_index - 1) * 10 ** 7
                    mass_boundary_at_end = fucntion_mass_boundary(age_of_this_epoch_at_end, mass_grid_table,
                                                                  lifetime_table)
                    all_sf_imf.append([igimf, mass_boundary_at_end, this_time])
                    time_of_the_epoch_in_Gyr = epoch_index / 100
                    all_sfr.append([S_F_R_of_this_epoch, time_of_the_epoch_in_Gyr])
                    epoch_info.append(
                        [S_F_R_of_this_epoch, M_tot_of_this_epoch, igimf_of_this_epoch, integrate_igimf_mass,
                         mass_grid_table, lifetime_table, Mfinal_table, mass_grid_table2, Mmetal_table, M_element_table,
                         last_time_age, SNIa_number_prob, metal_mass_fraction_in_gas, mass_calibration_factor])
                    metal_in_gas = metal_mass_fraction_in_gas
                else:  # if SFR == 0
                    time_of_the_epoch_in_Gyr = epoch_index / 100
                    all_sfr.append([10 ** -22, time_of_the_epoch_in_Gyr])
                    epoch_info.append(
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 0, 0, [0, 0, 0, 0, 0], 0])
            else:  # if epoch_index =! len(epoch_info)
                S_F_R_of_this_epoch = epoch_info[epoch_index][0]
                M_tot_of_this_epoch = epoch_info[epoch_index][1]
                igimf_of_this_epoch = epoch_info[epoch_index][2]
                integrate_igimf_mass = epoch_info[epoch_index][3]
                mass_grid_table = epoch_info[epoch_index][4]
                lifetime_table = epoch_info[epoch_index][5]
                Mfinal_table = epoch_info[epoch_index][6]
                mass_grid_table2 = epoch_info[epoch_index][7]
                Mmetal_table = epoch_info[epoch_index][8]
                M_element_table = epoch_info[epoch_index][9]
                last_time_age = epoch_info[epoch_index][10]
                epoch_info[epoch_index][10] = age_of_this_epoch
                SNIa_number_prob = epoch_info[epoch_index][11]
                metal_in_gas = epoch_info[epoch_index][12]
                mass_calibration_factor = epoch_info[epoch_index][13]
                def igimf_xi_function(mass):
                    return igimf_of_this_epoch.custom_imf(mass, this_time)
                def igimf_mass_function(mass):
                    return igimf_of_this_epoch.custom_imf(mass, this_time) * mass
                def igimf_luminous_function(mass):
                    return igimf_of_this_epoch.custom_imf(mass, this_time) * \
                           stellar_luminosity.stellar_luminosity_function(mass)

            if S_F_R_of_this_epoch > 0:
                # get M_tot (total initial mass of all star ever formed)
                M_tot_up_to_this_time += M_tot_of_this_epoch
                # calculate stellar initial mass that is still alive (dead star mass boundary)
                mass_boundary = fucntion_mass_boundary(age_of_this_epoch, mass_grid_table, lifetime_table)
                # output of this epoch
                # Mtarget_table_number:
                # 1: Mfinal_table
                # 2: Mmetal_table
                # 3: MH_table
                # 4: M_element_table
                # ...
                if integrate_igimf_mass != 0:

                    # m1 = quad(igimf_mass_function, 0.08, 10, limit=40)[0]
                    # m2 = quad(igimf_mass_function, 10, 150, limit=40)[0]
                    # print(m1)
                    # print(m2)
                    # print(m1 / m2)

                    inte_limit = max(round((math.log(mass_boundary, 10)+1) / (math.log(steller_mass_upper_bound, 10)+1) * 50), 20)

                    if imf == 'diet_Salpeter':
                        integrate_star_mass = quad(igimf_mass_function, 0.1, 100, limit=inte_limit)[0]  # normalized mass
                        stellar_luminosity_of_a_epoch_at_a_time_step = \
                        quad(igimf_luminous_function, 0.1, 100, limit=inte_limit)[0]
                    else:
                        integrate_star_mass = quad(igimf_mass_function, 0.08, mass_boundary, limit=inte_limit)[0]  # normalized mass
                        stellar_luminosity_of_a_epoch_at_a_time_step = \
                        quad(igimf_luminous_function, 0.08, mass_boundary, limit=inte_limit)[0]
                    stellar_mass_of_a_epoch_at_a_time_step = mass_calibration_factor * integrate_star_mass  # real mass

                    # apprent metal mass (neglect stellar evolution, only account for the initial metal mass when SF):
                    # the stellar metal abandance is the gas abdandance at the time of star foramtion (metal_in_gas).
                    stellar_metal_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[0]
                    stellar_H_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[1]
                    stellar_He_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[2]
                    stellar_C_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[3]
                    stellar_N_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[4]
                    stellar_O_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[5]
                    stellar_Mg_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[6]
                    stellar_Ca_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[7]
                    stellar_Fe_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[8]
                    stellar_Ne_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[9]
                    stellar_Si_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[10]
                    stellar_S_mass_of_this_epoch = stellar_mass_of_a_epoch_at_a_time_step * metal_in_gas[11]

                    # The luminosity-weighted metallicity is in its exact form. However,
                    # the luminosity-weighted element abundance, e.g., weighted-with-luminosity([Fe/H]) is approximated
                    # by [the-number-of(weighted-with-luminosity(mass-fraction-of(Fe)))/the-number-of(weighted-with-luminosity(mass-fraction-of(H)))]
                    # below is the first step to calculate the weighted-with-luminosity(mass-fraction-of(An-element))
                    stellar_metal_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[0]
                    stellar_H_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[1]
                    stellar_He_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[2]
                    stellar_C_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[3]
                    stellar_N_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[4]
                    stellar_O_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[5]
                    stellar_Mg_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[6]
                    stellar_Ca_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[7]
                    stellar_Fe_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[8]
                    stellar_Ne_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[9]
                    stellar_Si_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[10]
                    stellar_S_luminosity_of_this_epoch = stellar_luminosity_of_a_epoch_at_a_time_step * metal_in_gas[11]
                    #
                    BH_mass_of_this_epoch = get_BH_mass(mass_boundary, 1, 1, mass_calibration_factor,
                                                        steller_mass_upper_bound)
                    NS_mass_of_this_epoch = get_NS_mass(mass_boundary, 1, 1, mass_calibration_factor)
                    WD_mass_of_this_epoch = get_WD_mass(mass_boundary, 1, 1, mass_calibration_factor)
                    remnant_mass_of_this_epoch = WD_mass_of_this_epoch + NS_mass_of_this_epoch + BH_mass_of_this_epoch
                    # Note: M_tot_of_this_epoch =! ejected_gas_mass_of_this_epoch +
                    # stellar_mass_of_a_epoch_at_a_time_step + remnant_mass_of_this_epoch
                    # because the remnant_mass is a spline fitted value
                    # while metall mass ejection is calculated with M_metal = M_ini - M_final - M_H - M_He,
                    # where M_final is the remnant mass given by the stellar yield table.

                    #
                    # # consider direct black hole as in Heger et al. (2003) (maybe not self-consistant with the stellar evolution table)
                    # if mass_boundary > 100:
                    #     SNII_number_of_this_epoch_1 = quad(igimf_xi_function, mass_boundary, steller_mass_upper_bound, limit=50)[0]
                    #     SNII_number_of_this_epoch_2 = 0
                    # elif mass_boundary > 40:
                    #     SNII_number_of_this_epoch_1 = quad(igimf_xi_function, 100, steller_mass_upper_bound, limit=50)[0]
                    #     SNII_number_of_this_epoch_2 = 0
                    # elif mass_boundary > 8:
                    #     SNII_number_of_this_epoch_1 = quad(igimf_xi_function, 100, steller_mass_upper_bound, limit=50)[0]
                    #     SNII_number_of_this_epoch_2 = quad(igimf_xi_function, mass_boundary, 40, limit=40)[0]
                    # else:
                    #     SNII_number_of_this_epoch_1 = quad(igimf_xi_function, 100, steller_mass_upper_bound, limit=50)[0]
                    #     SNII_number_of_this_epoch_2 = quad(igimf_xi_function, 8, 40, limit=40)[0]
                    # SNII_number_of_this_epoch = (SNII_number_of_this_epoch_1 + SNII_number_of_this_epoch_2) * mass_calibration_factor
                    if mass_boundary > 8:
                        SNII_number_of_this_epoch = \
                        quad(igimf_xi_function, mass_boundary, steller_mass_upper_bound, limit=50)[0]
                        SNII_ejected_mass_of_this_epoch = \
                        quad(igimf_xi_function, mass_boundary, steller_mass_upper_bound, limit=50)[0]
                    else:
                        SNII_number_of_this_epoch = quad(igimf_xi_function, 8, steller_mass_upper_bound, limit=50)[0]
                    SNII_number_of_this_epoch = SNII_number_of_this_epoch * mass_calibration_factor
                    SNII_energy_release_per_event = 0.03* 10 ** 51  # Bradamante 1998
                    SNII_number += SNII_number_of_this_epoch
                    SNII_energy_release += SNII_energy_release_per_event * SNII_number_of_this_epoch
                    # ejected_ :
                    metal_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary,
                                                                                 steller_mass_upper_bound, 2, 2,
                                                                                 mass_calibration_factor)
                    H_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound, 2,
                                                                             "H", mass_calibration_factor)
                    He_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound,
                                                                              2, "He", mass_calibration_factor)
                    ejected_gas_mass_of_this_epoch = H_mass_of_this_epoch + He_mass_of_this_epoch + \
                                                     metal_mass_of_this_epoch
                    C_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound,
                                                                             2, "C", mass_calibration_factor)
                    N_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound,
                                                                             2, "N", mass_calibration_factor)
                    O_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound,
                                                                             2, "O", mass_calibration_factor)
                    Mg_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound,
                                                                              2, "Mg", mass_calibration_factor)
                    Ca_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound,
                                                                              2, "Ca", mass_calibration_factor)
                    S_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound,
                                                                              2, "S", mass_calibration_factor)
                    Si_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound,
                                                                              2, "Si", mass_calibration_factor)
                    Ne_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound,
                                                                              2, "Ne", mass_calibration_factor)
                    Fe_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound,
                                                                              2, "Fe", mass_calibration_factor)
                    Fe_production_SNII += Fe_mass_of_this_epoch
                    Ca_production_SNII += Ca_mass_of_this_epoch
                    Ne_production_SNII += Ne_mass_of_this_epoch
                    Si_production_SNII += Si_mass_of_this_epoch
                    S_production_SNII += S_mass_of_this_epoch
                    Mg_production_SNII += Mg_mass_of_this_epoch
                    O_production_SNII += O_mass_of_this_epoch
                    # if age_of_this_epoch == 1 * 10 ** 9:
                    #     print("Fe_production_SNII", Fe_production_SNII)
                    #     print("O_production_SNII", O_production_SNII)
                    #     print("Mg_production_SNII", Mg_production_SNII)
                    # _mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound, 2, "",
                    #                                                           mass_calibration_factor)


                else:
                    print("Error: integrate_igimf_mass == 0 while S_F_R_of_this_epoch != 0.")
                    stellar_mass_of_a_epoch_at_a_time_step = 0
                    BH_mass_of_this_epoch = 0
                    NS_mass_of_this_epoch = 0
                    WD_mass_of_this_epoch = 0
                    remnant_mass_of_this_epoch = 0
                    ejected_gas_mass_of_this_epoch = 0
                    metal_mass_of_this_epoch = 0
                    H_mass_of_this_epoch = 0
                    He_mass_of_this_epoch = 0
                    C_mass_of_this_epoch = 0
                    N_mass_of_this_epoch = 0
                    O_mass_of_this_epoch = 0
                    Mg_mass_of_this_epoch = 0
                    Ca_mass_of_this_epoch = 0
                    Si_mass_of_this_epoch = 0
                    S_mass_of_this_epoch = 0
                    Ne_mass_of_this_epoch = 0
                    Fe_mass_of_this_epoch = 0
                # if consider SNIa
                if SNIa_ON == True or 'power-law' or 'SD':
                    # read in SNIa yield table
                    # (here only account for the most abandant element yields)
                    # (but should account as long as the SNIa yield is comparable with SNII yield)
                    Fe_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Fe')
                    Si_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Si')
                    O_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'O')
                    S_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'S')
                    Mg_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Mg')
                    Ne_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Ne')
                    if SNIa_yield_table=='Seitenzahl2013' or SNIa_yield_table=='Iwamoto1999':
                        Ca_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Ca')
                        Ne_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Ne')
                        S_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'S')
                        Si_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Si')
                        C_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'C')
                    total_mass_eject_per_SNIa = Fe_mass_eject + Si_mass_eject + O_mass_eject + S_mass_eject + Mg_mass_eject + Ne_mass_eject
                    Chandrasekhar_mass = 1.44
                    pre_SNIa_NS_mass = 1
                    SNIa_energy_release_per_event = 0.8*10**51  # in the unit of 10^51 erg
                    # # Yin J., Matteucci F., Vladilo G., 2011, A&A, 531, A136
                    # integrate SNIa number from last_delay_time to this_delay_time contributed by this SF epoch
                    if SNIa_ON == 'SD':
                        if age_of_this_epoch == 0:
                            mass_boundary_SNIa_Padovani93 = 100
                            # mass_boundary_SNIa_Greggio83 = 8
                        else:
                            mass_boundary_SNIa_Padovani93 = 10 ** (7.764 - (1.79 - (1.338 - math.log(age_of_this_epoch, 10) * 0.1116) ** 2) / 0.2232)
                            # if age_of_this_epoch <  3 * 1e7:
                            #     mass_boundary_SNIa_Greggio83 = 8
                            # else:
                            #     mass_boundary_SNIa_Greggio83 = fucntion_mass_boundary_SNIa_Greggio83(age_of_this_epoch, 8)
                        # SNIa_number_from_this_epoch_till_this_time = function_number_SNIa_SD(mass_boundary_SNIa_Greggio83, igimf_xi_function, mass_calibration_factor)
                        # SNIa_number_from_this_epoch_till_this_time = function_number_SNIa_SD(mass_boundary_SNIa_Padovani93, igimf_xi_function, mass_calibration_factor)
                        SNIa_number_from_this_epoch_till_this_time = function_number_SNIa_SD(mass_boundary, igimf_xi_function, mass_calibration_factor)
                        # print('...', time_axis[time_step]/1e6, mass_boundary)
                        # print(SNIa_number_from_this_epoch_till_this_time)
                    elif SNIa_ON == True or 'power-law':
                        SNIa_number_from_this_epoch_till_this_time = function_number_SNIa_power_law(0, age_of_this_epoch,
                                                                                          SNIa_number_prob,
                                                                                          S_F_R_of_this_epoch)

                    # the following should result in 0.0022+-50% for a SSP,
                    # but now calibrate to a different value to fit with galaxy [Fe/H] observation
                    if age_of_this_epoch == 10 * 10 ** 9 - 1 * 10 ** 7:
                        # print(function_number_SNIa_power_law(0, 10 * 10 ** 9, 1, 0))
                        # print("SN number per star in range:", SNIa_number_from_this_epoch_till_this_time/number_in_SNIa_boundary)
                        print("\nType Ia supernova activated. "
                              "Total SNIa number per solar mass of star formed at t = 10Gyr:",
                              SNIa_number_from_this_epoch_till_this_time / M_tot_of_this_epoch)
                        # print("Total SNIa number between t = 9 and 10 Gyr:",
                        #       function_number_SNIa_power_law(9 * 10 ** 9 - 1 * 10 ** 7, age_of_this_epoch,
                        #                            number_in_SNIa_boundary,
                        #                            S_F_R_of_this_epoch))
                    # update the element masses
                    ejected_gas_mass_of_this_epoch += total_mass_eject_per_SNIa * SNIa_number_from_this_epoch_till_this_time
                    metal_mass_of_this_epoch += (Chandrasekhar_mass - (Chandrasekhar_mass - pre_SNIa_NS_mass) *
                                                 Z_gas_this_time_step) * SNIa_number_from_this_epoch_till_this_time
                    O_mass_of_SNIa = O_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    Mg_mass_of_SNIa = Mg_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    Fe_mass_of_SNIa = (Fe_mass_eject
                                       # - (Chandrasekhar_mass - pre_SNIa_NS_mass) * Fe_H_mass_ratio_at_last_time * 0.7057 # this term is small and can be neglected
                                       ) * SNIa_number_from_this_epoch_till_this_time
                    # Si_mass_of_SNIa = Si_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    # S_mass_of_SNIa = S_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    # Ne_mass_of_SNIa = Ne_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    if SNIa_yield_table == 'Seitenzahl2013' or SNIa_yield_table=='Iwamoto1999':
                        Ca_mass_of_SNIa = Ca_mass_eject * SNIa_number_from_this_epoch_till_this_time
                        Si_mass_of_SNIa = Si_mass_eject * SNIa_number_from_this_epoch_till_this_time
                        S_mass_of_SNIa = S_mass_eject * SNIa_number_from_this_epoch_till_this_time
                        Ne_mass_of_SNIa = Ne_mass_eject * SNIa_number_from_this_epoch_till_this_time
                        C_mass_of_SNIa = C_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    O_mass_of_this_epoch += O_mass_of_SNIa
                    Mg_mass_of_this_epoch += Mg_mass_of_SNIa
                    Fe_mass_of_this_epoch += Fe_mass_of_SNIa
                    # Si_mass_of_this_epoch += Si_mass_of_SNIa
                    # S_mass_of_this_epoch += S_mass_of_SNIa
                    # Ne_mass_of_this_epoch += Ne_mass_of_SNIa
                    if SNIa_yield_table == 'Seitenzahl2013' or SNIa_yield_table=='Iwamoto1999':
                        Ca_mass_of_this_epoch += Ca_mass_of_SNIa
                        Ne_mass_of_this_epoch += Ne_mass_of_SNIa
                        S_mass_of_this_epoch += S_mass_of_SNIa
                        Si_mass_of_this_epoch += Si_mass_of_SNIa
                        C_mass_of_this_epoch += C_mass_of_SNIa
                    remnant_mass_of_this_epoch -= pre_SNIa_NS_mass * SNIa_number_from_this_epoch_till_this_time
                    WD_mass_of_this_epoch -= pre_SNIa_NS_mass * SNIa_number_from_this_epoch_till_this_time
                    SNIa_number_from_all_epoch += SNIa_number_from_this_epoch_till_this_time
                    SNIa_energy_release += SNIa_energy_release_per_event * SNIa_number_from_this_epoch_till_this_time
                #
                stellar_mass_at_this_time += stellar_mass_of_a_epoch_at_a_time_step
                stellar_metal_mass_at_this_time += stellar_metal_mass_of_this_epoch
                stellar_H_mass_at_this_time += stellar_H_mass_of_this_epoch
                stellar_He_mass_at_this_time += stellar_He_mass_of_this_epoch
                stellar_O_mass_at_this_time += stellar_O_mass_of_this_epoch
                stellar_C_mass_at_this_time += stellar_C_mass_of_this_epoch
                stellar_N_mass_at_this_time += stellar_N_mass_of_this_epoch
                stellar_Ca_mass_at_this_time += stellar_Ca_mass_of_this_epoch
                stellar_Si_mass_at_this_time += stellar_Si_mass_of_this_epoch
                stellar_S_mass_at_this_time += stellar_S_mass_of_this_epoch
                stellar_Ne_mass_at_this_time += stellar_Ne_mass_of_this_epoch
                stellar_Mg_mass_at_this_time += stellar_Mg_mass_of_this_epoch
                stellar_Fe_mass_at_this_time += stellar_Fe_mass_of_this_epoch
                #
                # The luminosity-weighted element mass fraction is,
                # e.g., stellar_Fe_luminosity_at_this_time / stellar_luminosity_at_this_time
                stellar_luminosity_at_this_time += stellar_luminosity_of_a_epoch_at_a_time_step
                stellar_metal_luminosity_at_this_time += stellar_metal_luminosity_of_this_epoch
                stellar_H_luminosity_at_this_time += stellar_H_luminosity_of_this_epoch
                stellar_He_luminosity_at_this_time += stellar_He_luminosity_of_this_epoch
                stellar_O_luminosity_at_this_time += stellar_O_luminosity_of_this_epoch
                stellar_C_luminosity_at_this_time += stellar_C_luminosity_of_this_epoch
                stellar_N_luminosity_at_this_time += stellar_N_luminosity_of_this_epoch
                stellar_Ca_luminosity_at_this_time += stellar_Ca_luminosity_of_this_epoch
                stellar_Ne_luminosity_at_this_time += stellar_Ne_luminosity_of_this_epoch
                stellar_S_luminosity_at_this_time += stellar_S_luminosity_of_this_epoch
                stellar_Si_luminosity_at_this_time += stellar_Si_luminosity_of_this_epoch
                stellar_Mg_luminosity_at_this_time += stellar_Mg_luminosity_of_this_epoch
                stellar_Fe_luminosity_at_this_time += stellar_Fe_luminosity_of_this_epoch

                BH_mass_till_this_time += BH_mass_of_this_epoch
                NS_mass_till_this_time += NS_mass_of_this_epoch
                WD_mass_till_this_time += WD_mass_of_this_epoch
                remnant_mass_at_this_time += remnant_mass_of_this_epoch
                ejected_gas_mass_till_this_time += ejected_gas_mass_of_this_epoch
                ejected_metal_mass_till_this_time += metal_mass_of_this_epoch
                ejected_H_mass_till_this_time += H_mass_of_this_epoch
                ejected_He_mass_till_this_time += He_mass_of_this_epoch
                ejected_O_mass_till_this_time += O_mass_of_this_epoch
                ejected_C_mass_till_this_time += C_mass_of_this_epoch
                ejected_N_mass_till_this_time += N_mass_of_this_epoch
                ejected_Ca_mass_till_this_time += Ca_mass_of_this_epoch
                ejected_Si_mass_till_this_time += Si_mass_of_this_epoch
                ejected_S_mass_till_this_time += S_mass_of_this_epoch
                ejected_Ne_mass_till_this_time += Ne_mass_of_this_epoch
                ejected_Mg_mass_till_this_time += Mg_mass_of_this_epoch
                ejected_Fe_mass_till_this_time += Fe_mass_of_this_epoch
            # Goes to the next SF epoch until all SF event before this time step is accounted:
            (epoch_index) = (epoch_index + 1)
        # output of this time step
        total_energy_release = SNIa_energy_release + SNII_energy_release

        ### yeilds at this time step from all SF epoch:
        ejected_gas_mass_at_this_time = ejected_gas_mass_till_this_time - ejected_gas_mass_till_last_time
        ejected_metal_mass_at_this_time = ejected_metal_mass_till_this_time - ejected_metal_mass_till_last_time
        ejected_H_mass_at_this_time = ejected_H_mass_till_this_time - ejected_H_mass_till_last_time
        ejected_He_mass_at_this_time = ejected_He_mass_till_this_time - ejected_He_mass_till_last_time
        ejected_C_mass_at_this_time = ejected_C_mass_till_this_time - ejected_C_mass_till_last_time
        ejected_N_mass_at_this_time = ejected_N_mass_till_this_time - ejected_N_mass_till_last_time
        ejected_O_mass_at_this_time = ejected_O_mass_till_this_time - ejected_O_mass_till_last_time
        ejected_Mg_mass_at_this_time = ejected_Mg_mass_till_this_time - ejected_Mg_mass_till_last_time
        ejected_Ca_mass_at_this_time = ejected_Ca_mass_till_this_time - ejected_Ca_mass_till_last_time
        ejected_Ne_mass_at_this_time = ejected_Ne_mass_till_this_time - ejected_Ne_mass_till_last_time
        ejected_S_mass_at_this_time = ejected_S_mass_till_this_time - ejected_S_mass_till_last_time
        ejected_Si_mass_at_this_time = ejected_Si_mass_till_this_time - ejected_Si_mass_till_last_time
        ejected_Fe_mass_at_this_time = ejected_Fe_mass_till_this_time - ejected_Fe_mass_till_last_time
        ejected_gas_Mg_over_Fe_till_this_time = function_element_abundunce(solar_abu_table, "Mg", "Fe",
                                                                           ejected_Mg_mass_till_this_time,
                                                                           ejected_Fe_mass_till_this_time, False)
        ejected_gas_Mg_over_Fe_at_this_time = function_element_abundunce(solar_abu_table, "Mg", "Fe",
                                                                         ejected_Mg_mass_at_this_time,
                                                                         ejected_Fe_mass_at_this_time, True)
        M_tot_of_this_time = M_tot_up_to_this_time - M_tot_up_to_last_time  # new SF mass added at this time step
        #
        galaxy_mass_without_gas_at_this_time = stellar_mass_at_this_time + remnant_mass_at_this_time
        if galaxy_mass_without_gas_at_this_time == 0 or ejected_gas_mass_at_this_time == 0:
            expansion_factor_instantaneous = 1
            expansion_factor_adiabat = 1
        elif galaxy_mass_without_gas_at_this_time < ejected_gas_mass_at_this_time:
            Warning_galaxy_mass_ejected_gas_mass = True
            # Warning: galaxy_mass < ejected_gas_mass.
            # This is due to too large a timestep.
            # It is easy to aviod this issue by applying the "high_time_resolution=True"
            # but the simulation will take much longer time.
            expansion_factor_instantaneous = 10
            expansion_factor_adiabat = (
                                       galaxy_mass_without_gas_at_this_time + ejected_gas_mass_at_this_time) / galaxy_mass_without_gas_at_this_time
        else:
            expansion_factor_instantaneous = galaxy_mass_without_gas_at_this_time / (
            galaxy_mass_without_gas_at_this_time - ejected_gas_mass_at_this_time)
            expansion_factor_adiabat = (
                                       galaxy_mass_without_gas_at_this_time + ejected_gas_mass_at_this_time) / galaxy_mass_without_gas_at_this_time

        expansion_factor = 10 ** (
            (math.log(expansion_factor_instantaneous, 10) + math.log(expansion_factor_adiabat, 10)) / 2)


        # calculate the gravitational binding engergy:
        # gravitational_constant = 6.674
        # # galaxy_mass_without_gas_at_this_time, original_gas_mass, total_gas_mass_at_this_time, ejected_gas_mass_at_this_time
        # # gas_mass = max(ejected_gas_mass_at_this_time, 1)
        # # galaxy mass--radii relation adopted from Dabringhausen 2008 eq.4
        # Dabringhausen_2008_a = 2.95
        # Dabringhausen_2008_b = 0.596
        # initial_expansion_factor = 10000  # need change for every simulation. use the expansion_factor at final time
        # # initial_expansion_factor = expansion_factor_list[-1]
        # if expansion_factor_list == []:
        #     current_expansion_factor = initial_expansion_factor
        # else:
        #     current_expansion_factor = initial_expansion_factor - expansion_factor_list[-1]
        # # log_binding_energy = round(
        # #     math.log(3 / 5 * gravitational_constant * 1.989**2 / 3.086, 10) + 40 + (2 - Dabringhausen_2008_b) *
        # #     math.log(original_gas_mass, 10) - math.log(Dabringhausen_2008_a, 10) +
        # #     6 * Dabringhausen_2008_b + math.log(current_expansion_factor, 10), 3)
        # # # 40 = 30 (solar mass) * 2 - 11 (Gravitational constant) - 16 (pc to meter) + 7 (J to erg)
        # # # binding_energy = 10 ** log_binding_energy  # [erg]

        ### Element abundances in the gas phase (in solar unit):
        # log_binding_energy = 53.7-(7-math.log(original_gas_mass, 10))*2 # 52.8 # 52.74 #  #

        # print('log_binding_energy', log_binding_energy)
        # if total_energy_release > 0:
        #     print(epoch_index, "log total_energy_release =", math.log(total_energy_release, 10))
        if outflow is not None and total_energy_release > 0 and math.log(total_energy_release, 10) > log_binding_energy_initial:
            lockup_and_outflow_mass = M_tot_of_this_epoch * outflow  # lockup gas mass in BDs is about 4% thus neglected while the uniform outflow is often assumed to be the same value as the formed stellar mass.
            # print("lockup_and_outflow_mass", lockup_and_outflow_mass)
        else:
            lockup_and_outflow_mass = M_tot_of_this_epoch

        total_gas_mass_at_this_time = total_gas_mass_at_last_time - lockup_and_outflow_mass + ejected_gas_mass_at_this_time
        if total_gas_mass_at_this_time < 0.0001:
            total_gas_mass_at_this_time = 0.0001
        total_metal_mass_at_this_time = total_metal_mass_in_gas_at_last_time - lockup_and_outflow_mass * \
                                                                               Z_gas_this_time_step + ejected_metal_mass_at_this_time
        if total_metal_mass_at_this_time < 0.0001:
            total_metal_mass_at_this_time = 0.0001
        total_H_mass_at_this_time = total_H_mass_at_last_time - lockup_and_outflow_mass * (
            total_H_mass_at_last_time / total_gas_mass_at_last_time) + ejected_H_mass_at_this_time
        if total_H_mass_at_this_time < 0.0001:
            total_H_mass_at_this_time = 0.0001
        total_He_mass_at_this_time = total_He_mass_at_last_time - lockup_and_outflow_mass * (
            total_He_mass_at_last_time / total_gas_mass_at_last_time) + ejected_He_mass_at_this_time
        if total_He_mass_at_this_time < 0.0001:
            total_He_mass_at_this_time = 0.0001
        total_C_mass_at_this_time = total_C_mass_at_last_time - lockup_and_outflow_mass * (
            total_C_mass_at_last_time / total_gas_mass_at_last_time) + ejected_C_mass_at_this_time
        if total_C_mass_at_this_time < 0.0001:
            total_C_mass_at_this_time = 0.0001
        total_N_mass_at_this_time = total_N_mass_at_last_time - lockup_and_outflow_mass * (
            total_N_mass_at_last_time / total_gas_mass_at_last_time) + ejected_N_mass_at_this_time
        if total_N_mass_at_this_time < 0.0001:
            total_N_mass_at_this_time = 0.0001
        total_O_mass_at_this_time = total_O_mass_at_last_time - lockup_and_outflow_mass * (
            total_O_mass_at_last_time / total_gas_mass_at_last_time) + ejected_O_mass_at_this_time
        if total_O_mass_at_this_time < 0.0001:
            total_O_mass_at_this_time = 0.0001
        total_Mg_mass_at_this_time = total_Mg_mass_at_last_time - lockup_and_outflow_mass * (
            total_Mg_mass_at_last_time / total_gas_mass_at_last_time) + ejected_Mg_mass_at_this_time
        if total_Mg_mass_at_this_time < 0.0001:
            total_Mg_mass_at_this_time = 0.0001
        total_Ca_mass_at_this_time = total_Ca_mass_at_last_time - lockup_and_outflow_mass * (
            total_Ca_mass_at_last_time / total_gas_mass_at_last_time) + ejected_Ca_mass_at_this_time
        if total_Ca_mass_at_this_time < 0.0001:
            total_Ca_mass_at_this_time = 0.0001
        total_Ne_mass_at_this_time = total_Ne_mass_at_last_time - lockup_and_outflow_mass * (
            total_Ne_mass_at_last_time / total_gas_mass_at_last_time) + ejected_Ne_mass_at_this_time
        if total_Ne_mass_at_this_time < 0.0001:
            total_Ne_mass_at_this_time = 0.0001
        total_Si_mass_at_this_time = total_Si_mass_at_last_time - lockup_and_outflow_mass * (
            total_Si_mass_at_last_time / total_gas_mass_at_last_time) + ejected_Si_mass_at_this_time
        if total_Si_mass_at_this_time < 0.0001:
            total_Si_mass_at_this_time = 0.0001
        total_S_mass_at_this_time = total_S_mass_at_last_time - lockup_and_outflow_mass * (
            total_S_mass_at_last_time / total_gas_mass_at_last_time) + ejected_S_mass_at_this_time
        if total_S_mass_at_this_time < 0.0001:
            total_S_mass_at_this_time = 0.0001
        total_Fe_mass_at_this_time = total_Fe_mass_at_last_time - lockup_and_outflow_mass * (
            total_Fe_mass_at_last_time / total_gas_mass_at_last_time) + ejected_Fe_mass_at_this_time
        if total_Fe_mass_at_this_time < 0.0001:
            total_Fe_mass_at_this_time = 0.0001



        # calculate the kinetic energy of the gas if they have a uniform temperature of 2 keV:
        X_for_H = total_H_mass_at_this_time / total_gas_mass_at_this_time
        Y_for_He = total_He_mass_at_this_time / total_gas_mass_at_this_time
        Z_for_metal = total_metal_mass_at_this_time / total_gas_mass_at_this_time
        mean_molecular_weight = 1 / (2 * X_for_H + 3 / 4 * Y_for_He + Z_for_metal / 2) * \
                                element_weight_table.function_element_weight("H") / 6.022140857 / 1.9891
        # / 10**23 / 10**33 (i.e., 10**56) mean_molecular_weight in solar mass unit.
        log_mean_molecular_weight = math.log(mean_molecular_weight, 10) - 56  # log [M_sun]
        log_total_number_of_molecule = math.log(total_gas_mass_at_this_time,
                                                10) - log_mean_molecular_weight  # log [Number]
        # 1 [keV] = 1.60217662 * 10**(-9) [erg]
        log_energy_per_molecule = math.log(2 * 1.60217662, 10) - 9  # [erg]
        log_total_gas_kinetic_energy = log_total_number_of_molecule + log_energy_per_molecule  # log [erg]
        total_gas_kinetic_energy = 10 ** log_total_gas_kinetic_energy

        # if outflow is None:
        #     if total_energy_release == 0:
        #         outflow = None
        #     elif math.log(total_energy_release, 10) + 51 > log_binding_energy:
        #         outflow = True
        # elif outflow == True:
        #     if total_energy_release == 0:
        #         outflow = None
        #     elif math.log(total_energy_release, 10) + 51 < log_binding_energy:
        #         outflow = None
        #
        # if gas_infall == True:
        #     function_update_element_gas_infall()

        # gas metallicity_at_this_time = total_metal_mass_at_this_time (in gas) / total_gas_mass_at_this_time
        Z_over_X = math.log(total_metal_mass_at_this_time / total_H_mass_at_this_time, 10) - math.log(Z_solar / X_solar,
                                                                                                      10)
        # Fe_H_mass_ratio_at_this_time = total_Fe_mass_at_this_time / total_H_mass_at_this_time
        gas_X_at_this_time = X_for_H
        gas_Y_at_this_time = Y_for_He
        gas_Z_at_this_time = Z_for_metal
        O_over_H_number_ratio = function_element_abundunce(solar_abu_table, "O", "H",
                                                           total_O_mass_at_this_time, total_H_mass_at_this_time, False)
        Mg_over_H_number_ratio = function_element_abundunce(solar_abu_table, "Mg", "H",
                                                            total_Mg_mass_at_this_time, total_H_mass_at_this_time, False)
        C_over_H_number_ratio = function_element_abundunce(solar_abu_table, "C", "H",
                                                            total_C_mass_at_this_time, total_H_mass_at_this_time, False)
        N_over_H_number_ratio = function_element_abundunce(solar_abu_table, "N", "H",
                                                            total_N_mass_at_this_time, total_H_mass_at_this_time, False)
        Ca_over_H_number_ratio = function_element_abundunce(solar_abu_table, "Ca", "H",
                                                            total_Ca_mass_at_this_time, total_H_mass_at_this_time, False)
        S_over_H_number_ratio = function_element_abundunce(solar_abu_table, "S", "H",
                                                            total_S_mass_at_this_time, total_H_mass_at_this_time, False)
        Si_over_H_number_ratio = function_element_abundunce(solar_abu_table, "Si", "H",
                                                            total_Si_mass_at_this_time, total_H_mass_at_this_time, False)
        Ne_over_H_number_ratio = function_element_abundunce(solar_abu_table, "Ne", "H",
                                                            total_Ne_mass_at_this_time, total_H_mass_at_this_time, False)
        Fe_over_H_number_ratio = function_element_abundunce(solar_abu_table, "Fe", "H",
                                                            total_Fe_mass_at_this_time, total_H_mass_at_this_time, False)
        C_over_Fe_number_ratio = function_element_abundunce(solar_abu_table, "C", "Fe",
                                                            total_C_mass_at_this_time, total_Fe_mass_at_this_time, False)
        N_over_O_number_ratio = function_element_abundunce(solar_abu_table, "N", "O",
                                                            total_N_mass_at_this_time, total_Fe_mass_at_this_time, False)
        O_over_Fe_number_ratio = function_element_abundunce(solar_abu_table, "O", "Fe",
                                                            total_O_mass_at_this_time, total_Fe_mass_at_this_time, False)
        Mg_over_Fe_number_ratio = function_element_abundunce(solar_abu_table, "Mg", "Fe",
                                                             total_Mg_mass_at_this_time, total_Fe_mass_at_this_time, False)
        Ca_over_Fe_number_ratio = function_element_abundunce(solar_abu_table, "Ca", "Fe",
                                                             total_Ca_mass_at_this_time, total_Fe_mass_at_this_time, False)
        Ne_over_Fe_number_ratio = function_element_abundunce(solar_abu_table, "Ne", "Fe",
                                                             total_Ne_mass_at_this_time, total_Fe_mass_at_this_time, False)
        S_over_Fe_number_ratio = function_element_abundunce(solar_abu_table, "S", "Fe",
                                                             total_S_mass_at_this_time, total_Fe_mass_at_this_time, False)
        Si_over_Fe_number_ratio = function_element_abundunce(solar_abu_table, "Si", "Fe",
                                                             total_Si_mass_at_this_time, total_Fe_mass_at_this_time, False)

        ### Element abundances in of stars (consider only the metal of stellar surface, i.e., neglect stellar evolution
        # This raises errors from very low mass stars which are fully convective but may not be observationally important):
        ##### mass weighted abundances
        # (total metal in stars / total H in stars):
        if stellar_mass_at_this_time > 0:
            mass_weighted_stellar_X = stellar_H_mass_at_this_time / stellar_mass_at_this_time
        else:
            mass_weighted_stellar_X = primary_H_mass_fraction
        if stellar_mass_at_this_time > 0:
            mass_weighted_stellar_Y = stellar_He_mass_at_this_time / stellar_mass_at_this_time
        else:
            mass_weighted_stellar_Y = primary_He_mass_fraction
        if stellar_mass_at_this_time > 0:
            mass_weighted_stellar_Z = stellar_metal_mass_at_this_time / stellar_mass_at_this_time
        else:
            mass_weighted_stellar_Z = 1-primary_H_mass_fraction-primary_He_mass_fraction
        mass_weighted_stellar_O_over_H = function_element_abundunce(solar_abu_table, "O", "H",
                                                                    stellar_O_mass_at_this_time,
                                                                    stellar_H_mass_at_this_time, False)
        mass_weighted_stellar_Mg_over_H = function_element_abundunce(solar_abu_table, "Mg", "H",
                                                                     stellar_Mg_mass_at_this_time,
                                                                     stellar_H_mass_at_this_time, False)
        mass_weighted_stellar_C_over_H = function_element_abundunce(solar_abu_table, "C", "H",
                                                                     stellar_C_mass_at_this_time,
                                                                     stellar_H_mass_at_this_time, False)
        mass_weighted_stellar_N_over_H = function_element_abundunce(solar_abu_table, "N", "H",
                                                                     stellar_N_mass_at_this_time,
                                                                     stellar_H_mass_at_this_time, False)
        mass_weighted_stellar_Ca_over_H = function_element_abundunce(solar_abu_table, "Ca", "H",
                                                                     stellar_Ca_mass_at_this_time,
                                                                     stellar_H_mass_at_this_time, False)
        mass_weighted_stellar_Si_over_H = function_element_abundunce(solar_abu_table, "Si", "H",
                                                                     stellar_Si_mass_at_this_time,
                                                                     stellar_H_mass_at_this_time, False)
        mass_weighted_stellar_S_over_H = function_element_abundunce(solar_abu_table, "S", "H",
                                                                     stellar_S_mass_at_this_time,
                                                                     stellar_H_mass_at_this_time, False)
        mass_weighted_stellar_Ne_over_H = function_element_abundunce(solar_abu_table, "Ne", "H",
                                                                     stellar_Ne_mass_at_this_time,
                                                                     stellar_H_mass_at_this_time, False)
        mass_weighted_stellar_Fe_over_H = function_element_abundunce(solar_abu_table, "Fe", "H",
                                                                     stellar_Fe_mass_at_this_time,
                                                                     stellar_H_mass_at_this_time, False)
        mass_weighted_stellar_C_over_Fe = function_element_abundunce(solar_abu_table, "C", "Fe",
                                                                     stellar_C_mass_at_this_time,
                                                                     stellar_Fe_mass_at_this_time, False)
        mass_weighted_stellar_N_over_O = function_element_abundunce(solar_abu_table, "N", "O",
                                                                     stellar_N_mass_at_this_time,
                                                                     stellar_Fe_mass_at_this_time, False)
        mass_weighted_stellar_O_over_Fe = function_element_abundunce(solar_abu_table, "O", "Fe",
                                                                     stellar_O_mass_at_this_time,
                                                                     stellar_Fe_mass_at_this_time, False)
        mass_weighted_stellar_Mg_over_Fe = function_element_abundunce(solar_abu_table, "Mg", "Fe",
                                                                      stellar_Mg_mass_at_this_time,
                                                                      stellar_Fe_mass_at_this_time, False)
        mass_weighted_stellar_Ca_over_Fe = function_element_abundunce(solar_abu_table, "Ca", "Fe",
                                                                      stellar_Ca_mass_at_this_time,
                                                                      stellar_Fe_mass_at_this_time, False)
        mass_weighted_stellar_Ne_over_Fe = function_element_abundunce(solar_abu_table, "Ne", "Fe",
                                                                      stellar_Ne_mass_at_this_time,
                                                                      stellar_Fe_mass_at_this_time, False)
        mass_weighted_stellar_S_over_Fe = function_element_abundunce(solar_abu_table, "S", "Fe",
                                                                      stellar_S_mass_at_this_time,
                                                                      stellar_Fe_mass_at_this_time, False)
        mass_weighted_stellar_Si_over_Fe = function_element_abundunce(solar_abu_table, "Si", "Fe",
                                                                      stellar_Si_mass_at_this_time,
                                                                      stellar_Fe_mass_at_this_time, False)
        ##### luminosity weighted abundances
        # (total metal in stars / total H in stars):
        if stellar_luminosity_at_this_time > 0:
            luminosity_weighted_stellar_X = stellar_H_luminosity_at_this_time / stellar_luminosity_at_this_time
        else:
            luminosity_weighted_stellar_X = primary_H_mass_fraction
        if stellar_luminosity_at_this_time > 0:
            luminosity_weighted_stellar_Y = stellar_He_luminosity_at_this_time / stellar_luminosity_at_this_time
        else:
            luminosity_weighted_stellar_Y = primary_He_mass_fraction
        if stellar_luminosity_at_this_time > 0:
            luminosity_weighted_stellar_Z = stellar_metal_luminosity_at_this_time / stellar_luminosity_at_this_time
        else:
            luminosity_weighted_stellar_Z = 1-primary_H_mass_fraction-primary_He_mass_fraction
            # below the input shall be the luminosity-weighted element mass,
        # e.g., stellar_O_luminosity_at_this_time / stellar_luminosity_at_this_time * total-stellar-mass-at-this-time,
        # but since stellar_luminosity_at_this_time and total-stellar-mass-at-this-time are the same for both element,
        # the constants cancel in function_element_abundunce.
        luminosity_weighted_stellar_O_over_H = function_element_abundunce(solar_abu_table, "O", "H",
                                                                          stellar_O_luminosity_at_this_time,
                                                                          stellar_H_luminosity_at_this_time, False)
        luminosity_weighted_stellar_Mg_over_H = function_element_abundunce(solar_abu_table, "Mg", "H",
                                                                           stellar_Mg_luminosity_at_this_time,
                                                                           stellar_H_luminosity_at_this_time, False)
        luminosity_weighted_stellar_C_over_H = function_element_abundunce(solar_abu_table, "C", "H",
                                                                           stellar_C_luminosity_at_this_time,
                                                                           stellar_H_luminosity_at_this_time, False)
        luminosity_weighted_stellar_N_over_H = function_element_abundunce(solar_abu_table, "N", "H",
                                                                           stellar_N_luminosity_at_this_time,
                                                                           stellar_H_luminosity_at_this_time, False)
        luminosity_weighted_stellar_Ca_over_H = function_element_abundunce(solar_abu_table, "Ca", "H",
                                                                           stellar_Ca_luminosity_at_this_time,
                                                                           stellar_H_luminosity_at_this_time, False)
        luminosity_weighted_stellar_Si_over_H = function_element_abundunce(solar_abu_table, "Si", "H",
                                                                           stellar_Si_luminosity_at_this_time,
                                                                           stellar_H_luminosity_at_this_time, False)
        luminosity_weighted_stellar_S_over_H = function_element_abundunce(solar_abu_table, "S", "H",
                                                                           stellar_S_luminosity_at_this_time,
                                                                           stellar_H_luminosity_at_this_time, False)
        luminosity_weighted_stellar_Ne_over_H = function_element_abundunce(solar_abu_table, "Ne", "H",
                                                                           stellar_Ne_luminosity_at_this_time,
                                                                           stellar_H_luminosity_at_this_time, False)
        luminosity_weighted_stellar_Fe_over_H = function_element_abundunce(solar_abu_table, "Fe", "H",
                                                                           stellar_Fe_luminosity_at_this_time,
                                                                           stellar_H_luminosity_at_this_time, False)
        luminosity_weighted_stellar_C_over_Fe = function_element_abundunce(solar_abu_table, "C", "Fe",
                                                                           stellar_C_luminosity_at_this_time,
                                                                           stellar_Fe_luminosity_at_this_time, False)
        luminosity_weighted_stellar_N_over_O = function_element_abundunce(solar_abu_table, "N", "O",
                                                                           stellar_N_luminosity_at_this_time,
                                                                           stellar_Fe_luminosity_at_this_time, False)
        luminosity_weighted_stellar_O_over_Fe = function_element_abundunce(solar_abu_table, "O", "Fe",
                                                                           stellar_O_luminosity_at_this_time,
                                                                           stellar_Fe_luminosity_at_this_time, False)
        luminosity_weighted_stellar_Mg_over_Fe = function_element_abundunce(solar_abu_table, "Mg", "Fe",
                                                                            stellar_Mg_luminosity_at_this_time,
                                                                            stellar_Fe_luminosity_at_this_time, False)
        luminosity_weighted_stellar_Ca_over_Fe = function_element_abundunce(solar_abu_table, "Ca", "Fe",
                                                                            stellar_Ca_luminosity_at_this_time,
                                                                            stellar_Fe_luminosity_at_this_time, False)
        luminosity_weighted_stellar_Ne_over_Fe = function_element_abundunce(solar_abu_table, "Ne", "Fe",
                                                                            stellar_Ne_luminosity_at_this_time,
                                                                            stellar_Fe_luminosity_at_this_time, False)
        luminosity_weighted_stellar_S_over_Fe = function_element_abundunce(solar_abu_table, "S", "Fe",
                                                                            stellar_S_luminosity_at_this_time,
                                                                            stellar_Fe_luminosity_at_this_time, False)
        luminosity_weighted_stellar_Si_over_Fe = function_element_abundunce(solar_abu_table, "Si", "Fe",
                                                                            stellar_Si_luminosity_at_this_time,
                                                                            stellar_Fe_luminosity_at_this_time, False)


        if stellar_H_mass_at_this_time == 0:
            mass_weighted_stellar_Z_over_X = math.log(Z_0 / Z_solar, 10)  # approximated with [Z]
            luminosity_weighted_stellar_Z_over_X = mass_weighted_stellar_Z_over_X
        else:
            mass_weighted_stellar_Z_over_X = math.log(stellar_metal_mass_at_this_time / stellar_H_mass_at_this_time, 10) \
                                             - math.log(Z_solar / X_solar, 10)
            luminosity_weighted_stellar_Z_over_X = math.log(
                stellar_metal_luminosity_at_this_time / stellar_H_luminosity_at_this_time, 10) \
                                                   - math.log(Z_solar / X_solar, 10)

        # the so called [Z/H] as determined by observation with equation: [Z/H] = [Fe/H] + A[Mg/Fe] where A=0.94 (Thomas 2003)
        if mass_weighted_stellar_Fe_over_H is None or mass_weighted_stellar_Mg_over_Fe is None:
            mass_weighted_stellar_Z_over_H = None
        else:
            mass_weighted_stellar_Z_over_H = mass_weighted_stellar_Fe_over_H + 0.94 * mass_weighted_stellar_Mg_over_Fe
            luminosity_weighted_stellar_Z_over_H = luminosity_weighted_stellar_Fe_over_H + 0.94 * luminosity_weighted_stellar_Mg_over_Fe

        ##########################################
        ##### luminosity weighted abundances #####
        ##########################################

        if BH_mass_till_this_time == 0:
            BH_mass_list += [10 ** (-10)]
        else:
            BH_mass_list += [BH_mass_till_this_time]

        if NS_mass_till_this_time == 0:
            NS_mass_list += [10 ** (-10)]
        else:
            NS_mass_list += [NS_mass_till_this_time]

        if WD_mass_till_this_time == 0:
            WD_mass_list += [10 ** (-10)]
        elif WD_mass_till_this_time < 0:
            Warning_WD_mass_till_this_time = True
            # Warning: more SNIa formed than WD avaliable. Please modify the SNIa rate assumption
            WD_mass_list += [10 ** (-10)]
        else:
            WD_mass_list += [WD_mass_till_this_time]

        if Z_over_X == 0:
            print("Warning: Z_over_X == 0")
            gas_Z_over_X_list += [math.log(Z_0 / Z_solar, 10)]  # approximated with [Z]
        else:
            gas_Z_over_X_list += [Z_over_X]

        ejected_O_mass_till_this_time_tot_list += [ejected_O_mass_till_this_time]
        ejected_O_mass_till_this_time_SNII_list += [O_production_SNII]
        ejected_O_mass_till_this_time_SNIa_list += [ejected_O_mass_till_this_time - O_production_SNII]

        ejected_Mg_mass_till_this_time_tot_list += [ejected_Mg_mass_till_this_time]
        ejected_Mg_mass_till_this_time_SNII_list += [Mg_production_SNII]
        ejected_Mg_mass_till_this_time_SNIa_list += [ejected_Mg_mass_till_this_time - Mg_production_SNII]

        ejected_Fe_mass_till_this_time_tot_list += [ejected_Fe_mass_till_this_time]
        ejected_Fe_mass_till_this_time_SNII_list += [Fe_production_SNII]
        ejected_Fe_mass_till_this_time_SNIa_list += [ejected_Fe_mass_till_this_time - Fe_production_SNII]

        ejected_Ca_mass_till_this_time_tot_list += [ejected_Ca_mass_till_this_time]
        ejected_Ca_mass_till_this_time_SNII_list += [Ca_production_SNII]
        ejected_Ca_mass_till_this_time_SNIa_list += [ejected_Ca_mass_till_this_time - Ca_production_SNII]

        ejected_S_mass_till_this_time_tot_list += [ejected_S_mass_till_this_time]
        ejected_S_mass_till_this_time_SNII_list += [S_production_SNII]
        ejected_S_mass_till_this_time_SNIa_list += [ejected_S_mass_till_this_time - S_production_SNII]

        ejected_Si_mass_till_this_time_tot_list += [ejected_Si_mass_till_this_time]
        ejected_Si_mass_till_this_time_SNII_list += [Si_production_SNII]
        ejected_Si_mass_till_this_time_SNIa_list += [ejected_Si_mass_till_this_time - Si_production_SNII]

        ejected_Ne_mass_till_this_time_tot_list += [ejected_Ne_mass_till_this_time]
        ejected_Ne_mass_till_this_time_SNII_list += [Ne_production_SNII]
        ejected_Ne_mass_till_this_time_SNIa_list += [ejected_Ne_mass_till_this_time - Ne_production_SNII]

        X_list += [gas_X_at_this_time]
        Y_list += [gas_Y_at_this_time]
        Z_list += [gas_Z_at_this_time]
        O_over_H_list += [O_over_H_number_ratio]
        Mg_over_H_list += [Mg_over_H_number_ratio]
        C_over_H_list += [C_over_H_number_ratio]
        N_over_H_list += [N_over_H_number_ratio]
        Ca_over_H_list += [Ca_over_H_number_ratio]
        Si_over_H_list += [Si_over_H_number_ratio]
        S_over_H_list += [S_over_H_number_ratio]
        Ne_over_H_list += [Ne_over_H_number_ratio]
        Fe_over_H_list += [Fe_over_H_number_ratio]
        C_over_Fe_list += [C_over_Fe_number_ratio]
        N_over_O_list += [N_over_O_number_ratio]
        O_over_Fe_list += [O_over_Fe_number_ratio]
        Mg_over_Fe_list += [Mg_over_Fe_number_ratio]
        Ca_over_Fe_list += [Ca_over_Fe_number_ratio]
        Ne_over_Fe_list += [Ne_over_Fe_number_ratio]
        S_over_Fe_list += [S_over_Fe_number_ratio]
        Si_over_Fe_list += [Si_over_Fe_number_ratio]

        stellar_O_over_H_list += [mass_weighted_stellar_O_over_H]
        stellar_Mg_over_H_list += [mass_weighted_stellar_Mg_over_H]
        stellar_C_over_H_list += [mass_weighted_stellar_C_over_H]
        stellar_N_over_H_list += [mass_weighted_stellar_N_over_H]
        stellar_Ca_over_H_list += [mass_weighted_stellar_Ca_over_H]
        stellar_Si_over_H_list += [mass_weighted_stellar_Si_over_H]
        stellar_S_over_H_list += [mass_weighted_stellar_S_over_H]
        stellar_Ne_over_H_list += [mass_weighted_stellar_Ne_over_H]
        stellar_X_list += [mass_weighted_stellar_X]
        stellar_Y_list += [mass_weighted_stellar_Y]
        stellar_Z_list += [mass_weighted_stellar_Z]
        stellar_Fe_over_H_list += [mass_weighted_stellar_Fe_over_H]
        stellar_C_over_Fe_list += [mass_weighted_stellar_C_over_Fe]
        stellar_N_over_O_list += [mass_weighted_stellar_N_over_O]
        stellar_O_over_Fe_list += [mass_weighted_stellar_O_over_Fe]
        stellar_Mg_over_Fe_list += [mass_weighted_stellar_Mg_over_Fe]
        stellar_Ca_over_Fe_list += [mass_weighted_stellar_Ca_over_Fe]
        stellar_Ne_over_Fe_list += [mass_weighted_stellar_Ne_over_Fe]
        stellar_S_over_Fe_list += [mass_weighted_stellar_S_over_Fe]
        stellar_Si_over_Fe_list += [mass_weighted_stellar_Si_over_Fe]
        stellar_Z_over_X_list += [mass_weighted_stellar_Z_over_X]
        stellar_Z_over_H_list += [mass_weighted_stellar_Z_over_H]

        stellar_O_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_O_over_H]
        stellar_Mg_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_Mg_over_H]
        stellar_C_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_C_over_H]
        stellar_N_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_N_over_H]
        stellar_Ca_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_Ca_over_H]
        stellar_Si_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_Si_over_H]
        stellar_S_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_S_over_H]
        stellar_Ne_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_Ne_over_H]
        stellar_X_list_luminosity_weighted += [luminosity_weighted_stellar_X]
        stellar_Y_list_luminosity_weighted += [luminosity_weighted_stellar_Y]
        stellar_Z_list_luminosity_weighted += [luminosity_weighted_stellar_Z]
        stellar_Fe_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_Fe_over_H]
        stellar_C_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_C_over_Fe]
        stellar_N_over_O_list_luminosity_weighted += [luminosity_weighted_stellar_N_over_O]
        stellar_O_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_O_over_Fe]
        stellar_Mg_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_Mg_over_Fe]
        stellar_Ca_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_Ca_over_Fe]
        stellar_Ne_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_Ne_over_Fe]
        stellar_S_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_S_over_Fe]
        stellar_Si_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_Si_over_Fe]
        stellar_Z_over_X_list_luminosity_weighted += [luminosity_weighted_stellar_Z_over_X]
        stellar_Z_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_Z_over_H]

        if remnant_mass_at_this_time == 0:
            remnant_mass_list += [10 ** (-10)]
        else:
            remnant_mass_list += [remnant_mass_at_this_time]

        if total_gas_mass_at_this_time == 0:
            total_gas_mass_list += [10 ** (-10)]
        else:
            total_gas_mass_list += [total_gas_mass_at_this_time]

        if ejected_gas_mass_till_this_time == 0 or ejected_gas_mass_till_this_time < 0:
            ejected_gas_mass_list += [10 ** (-10)]
        else:
            ejected_gas_mass_list += [ejected_gas_mass_till_this_time]

        ejected_gas_Mg_over_Fe_list += [ejected_gas_Mg_over_Fe_till_this_time]
        instant_ejected_gas_Mg_over_Fe_list += [ejected_gas_Mg_over_Fe_at_this_time]

        if M_tot_up_to_this_time > 0:
            ejected_metal_mass_list += [ejected_metal_mass_till_this_time / M_tot_up_to_this_time]
        else:
            ejected_metal_mass_list += [0]

        if expansion_factor_instantaneous_list == []:
            expansion_factor_instantaneous_list += [1]
        else:
            expansion_factor_instantaneous_list += [
                expansion_factor_instantaneous * expansion_factor_instantaneous_list[-1]]

        if expansion_factor_adiabat_list == []:
            expansion_factor_adiabat_list += [1]
        else:
            expansion_factor_adiabat_list += [expansion_factor_adiabat * expansion_factor_adiabat_list[-1]]

        if expansion_factor_list == []:
            expansion_factor_list += [1]
        else:
            expansion_factor_list += [expansion_factor * expansion_factor_list[-1]]

        if stellar_mass_at_this_time == 0:
            stellar_mass_list += [10 ** (-10)]
        else:
            stellar_mass_list += [stellar_mass_at_this_time]

        SNIa_energy_release_list += [SNIa_energy_release]
        SNIa_number_list += [SNIa_number_from_all_epoch]

        if len(SNIa_number_per_century) == 0:
            SNIa_number_per_century += [SNIa_number_list[0]]
        else:
            SNIa_number_per_century += [
                (SNIa_number_list[-1] - SNIa_number_list[-2]) / (time_axis[time_step] - time_axis[time_step - 1]) * 100]

        SNII_energy_release_list += [SNII_energy_release]
        SNII_number_list += [SNII_number]

        if len(SNII_number_per_century) == 0:
            SNII_number_per_century += [SNII_number_list[0]]
        else:
            SNII_number_per_century += [
                (SNII_number_list[-1] - SNII_number_list[-2]) / (time_axis[time_step] - time_axis[time_step - 1]) * 100]

        if total_energy_release == 0:
            total_energy_release_list += [0]
        else:
            total_energy_release_list += [total_energy_release]  # [math.log((total_energy_release), 10)]

        # if binding_energy == 0:
        #     binding_energy_list += [0]
        # else:
        #     binding_energy_list += [binding_energy]#[math.log((binding_energy), 10)]

        if total_gas_kinetic_energy_list == 0:
            total_gas_kinetic_energy_list += [0]
        else:
            total_gas_kinetic_energy_list += [total_gas_kinetic_energy]  # [math.log((total_gas_kinetic_energy), 10)]

        SN_number_per_century += [SNIa_number_per_century[-1] + SNII_number_per_century[-1]]

        # go to next time step
        if time_step / 50 > gc_collect_check:
            gc_collect_check += 1
            print("gc_collect:", gc_collect_check)
            gc.collect()
        (time_step) = (time_step + 1)

    ######################
    ### Show Warnings ###
    ######################

    if Warning_ejected_gas_mass_of_this_epoch == True:
        print('Warning: ejected_gas_mass_of_this_epoch < 0. See comments in galevo.py')

    if Warning_WD_mass_till_this_time == True:
        print("Warning: WD_mass_till_this_time < 0. See comments in galevo.py")

    if Warning_galaxy_mass_ejected_gas_mass == True:
        print('Warning: galaxy_mass < ejected_gas_mass. See comments in galevo.py.')
        # Warning: galaxy_mass < ejected_gas_mass.
        # This is due to too large a timestep.
        # It is easy to aviod this issue by applying the "high_time_resolution=True"
        # but the simulation will take much longer time.

    computation_time_seconds = round((time.time() - start_time), 2)
    minutes, seconds = divmod(computation_time_seconds, 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    days = int(days)
    hours = int(hours)
    minutes = int(minutes)
    seconds = round(seconds, 4)
    print("- Simulation complete. Computation time: {} d {} h {} m {} s -".format(days, hours, minutes, seconds))

    ###################
    ### output data ###
    ###################

    # Remnant_Star_ratio = [0]*len(stellar_mass_list)
    # for i in range(len(remnant_mass_list)):
    #     Remnant_Star_ratio[i] = remnant_mass_list[i]/stellar_mass_list[i]
    # import csv
    # with open('GalEvo_time.txt', 'w') as f:
    #     writer = csv.writer(f, delimiter=' ')
    #     f.write("# galevo.py output file.\n# time\n")
    #     writer.writerows(
    #         zip(time_axis))
    # with open('GalEvo_ratio.txt', 'w') as f:
    #     writer = csv.writer(f, delimiter=' ')
    #     f.write("# galevo.py output file.\n# Remnant_Star_ratio\n")
    #     writer.writerows(
    #         zip(Remnant_Star_ratio))


    ###################
    ###    output   ###
    ###################
    log_Z_0 = round(math.log(Z_0 / Z_solar, 10), 2)

    text_output(imf, STF, round(math.log(max(SFH_input), 10), 1), SFEN, original_gas_mass, log_Z_0)

    # if output plot applies
    plot_output(plot_show, plot_save, imf, igimf, round(math.log(max(SFH_input), 10), 1), SFEN, log_Z_0, STF)

    ###################
    ###     end     ###
    ###################

    return


# def function_update_element_gas_infall():
#     return


# # # calculate the diet_Salpeter_mass_to_number_ratio:
# # Bell & de Jong (2001). Salpeter IMF x = 1.35 with a flat x = 0 slope below 0.35
# def function_xi_diet_Salpeter_IMF(mass):
#     # integrate this function's output xi result in the number of stars in mass limits.
#     xi = diet_Salpeter_IMF.custom_imf(mass, 0)
#     return xi


# def function_mass_diet_Salpeter_IMF(mass):
#     # integrate this function's output m result in the total stellar mass for stars in mass limits.
#     m = mass * diet_Salpeter_IMF.custom_imf(mass, 0)
#     return m


# integrate_all_for_function_mass_SNIa = quad(function_mass_diet_Salpeter_IMF, 0.1, 100, limit=50)[0]
# integrate_28_for_function_number_SNIa = quad(function_xi_diet_Salpeter_IMF, 3, 8, limit=50)[0]
# diet_Salpeter_mass_to_number_ratio = integrate_all_for_function_mass_SNIa / integrate_28_for_function_number_SNIa

def function_xi_Kroupa_IMF(mass):  # there is no time dependence for Kroupa IMF
    if mass < 0.08:
        return 0
    elif mass < 0.5:
        return 2*mass**(-1.3)
    elif mass < 150:
        return mass**(-2.3)
    else:
        return 0

def function_mass_Kroupa_IMF(mass):
    # integrate this function's output m result in the total stellar mass for stars in mass limits.
    m = mass * function_xi_Kroupa_IMF(mass)
    return m

integrate_all_for_function_mass_SNIa = quad(function_mass_Kroupa_IMF, 0.08, 150, limit=50)[0]
integrate_total_number_SNIa = quad(function_xi_Kroupa_IMF, 0.08, 150, limit=50)[0]
integrate_28_for_function_number_SNIa = quad(function_xi_Kroupa_IMF, 3, 8, limit=50)[0]
SNIa_number_prob_Kroupa = integrate_28_for_function_number_SNIa ** 2 / integrate_total_number_SNIa /integrate_all_for_function_mass_SNIa


def function_number_SNIa_power_law(last_delay_time, this_delay_time, SNIa_number_prob__, S_F_R_of_this_epoch):
    # This function calculate the number of SNIa between last_delay_time and this_delay_time

    # It is commonly assumed that the maximum stellar mass able to produce a degenerate C–O white dwarf is 8 M⊙,
    # The minimum possible binary mass is assumed to be 3 M⊙ in order to ensure that the
    # smallest possible white dwarf can accrete enough mass from the secondary star to reach the Chandrasekhar mass.
    # see Greggio, L., & Renzini, A. 1983, A & A, 118, 217
    # Thus we should normalize the DTD according to the number (but currently, mass) of stars between 1.5 and 8 solar mass
    # normalized with a SNIa assuming fixed diet-Salpeter IMF (Bell et al. 149:289–312, 2003)
    # See Dan Maoz and Filippo Mannucci 2012 review
    global diet_Salpeter_mass_to_number_ratio
    SNIa_normalization_parameter = funtion_SNIa_DTD_normalization_parameter(S_F_R_of_this_epoch)
    # SNIa_normalization_parameter considers the possible variation of binary encounter rate in different system density
    # integrate SNIa number from last_delay_time to this_delay_time using observationally determined DTD assuming diet-Salpeter IMF
    SNIa_number_per_solar_mass = quad(function_SNIa_DTD, last_delay_time, this_delay_time, limit=40)[0]
    # calculate actual SNIa event number
    # SNIa_number = stellar_number_in_SNIa_boundary * SNIa_normalization_parameter * SNIa_number_per_solar_mass * diet_Salpeter_mass_to_number_ratio
    SNIa_number = S_F_R_of_this_epoch * 10**7 * SNIa_number_per_solar_mass / SNIa_number_prob_Kroupa * SNIa_number_prob__ * SNIa_normalization_parameter
    return SNIa_number

def function_number_SNIa_SD(mass_boundary, igimf_xi_function, mass_calibration_factor):
    # this function calculate the total number of SNIa events orignicated from a single 10 Myr star formation epoch
    # (i.e. a burst) since the birth of this stellar population till the current time step
    # The age of the population corresponds to a stellar mass (mass_boundary) with a lifetime equal to this age
    if mass_boundary > 8:
        SNIa_number = 0
    else:
        M2min = max(mass_boundary, 0.8)
        M2max = 8
        A_SNIa_Matteucci01 = 0.006 * 14 / 0.0024077124353644908 * 0.002
        SNIa_number = mass_calibration_factor * A_SNIa_Matteucci01 * quad(function_SNIa_SD_xi_2, M2min, M2max, args=(igimf_xi_function))[0]
    # print(SNIa_number)
    return SNIa_number

def function_SNIa_SD_xi_2(M_2, igimf_xi_function):
    # M_B_min = 3.9
    # M_1_min = 3
    # M_B_inf = max((2*M_2), M_B_min, M_1_min+M_2)
    # M_up = min((8 + M_2), 9.5)
    M_B_min = 3
    M_B_inf = max((2*M_2), M_B_min)
    M_up = 8 + M_2
    xi_2 = quad(function_SD_f_times_xi, M_B_inf, M_up, args=(M_2, igimf_xi_function))[0]
    return xi_2

def function_SD_f_times_xi(M_B, M_2, igimf_xi_function):
    gamma_f = 2
    mu = M_2 / M_B
    if mu > 0.5 or mu == 0 or mu < 0:
        f_times_xi = 0  # since f = 0
    else:
        # f_times_xi = (mu**gamma_f) * igimf_xi_function(M_B)
        f_times_xi = (2**(1+gamma_f) * (1+gamma_f) * mu**gamma_f) * igimf_xi_function(M_B) / M_B
    return f_times_xi

def funtion_SNIa_DTD_normalization_parameter(SFR):
    # this modification on the SNIa rate is to honor the fact that the number of SNIa should
    # not only depends on the number of potential progenitor but also the density of the stellar system
    # as is expected by the dynamical encounter rate.
    # x = 0.95
    # xplusSFR = SFR / 16 + x
    # gamma_DTD = 0.8 / (0.5 + math.log(xplusSFR, 10))
    # output = xplusSFR ** gamma_DTD
    output = 1
    # This is a toy model relation. Roughly keep the SNIa-rate/stellar-mass-formed unchanged for different SFRs (IMFs).
    # When SFR is high, top-heavy IMF reduce the number possible SNIa progenitor within mass range 1.5-8 solar mass.
    # the dense stellar region pump the SNIa rate as the stars have a larger chance to meet with another star.
    # The lower SFR do not further reduce the SNIa number as the low SFR happens after the star burst epoch
    # thus the newly formed star can still meet with stars formed at ealier epochs regredless of its current SFR.
    return output

####### the following code is for test, showing the renormalization function of SNIa# #######

# xxx = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
# y0 = funtion_SNIa_DTD_normalization_parameter(0.0001)
# yyy = [1, funtion_SNIa_DTD_normalization_parameter(0.001)/y0,
#        funtion_SNIa_DTD_normalization_parameter(0.01)/y0, funtion_SNIa_DTD_normalization_parameter(0.1)/y0,
#        funtion_SNIa_DTD_normalization_parameter(1)/y0, funtion_SNIa_DTD_normalization_parameter(10)/y0,
#        funtion_SNIa_DTD_normalization_parameter(100)/y0, funtion_SNIa_DTD_normalization_parameter(1000)/y0,
#        funtion_SNIa_DTD_normalization_parameter(10000)/y0]
#
# plt.rc('font', family='serif')
# plt.rc('xtick', labelsize='x-small')
# plt.rc('ytick', labelsize='x-small')
# fig = plt.figure(200, figsize=(3, 2.5))
# fig.add_subplot(1, 1, 1)
# plt.plot(xxx, yyy)
# plt.xlabel(r'log$_{10}$(SFR [$M_\odot$/yr])')
# plt.ylabel(r'SNIa Number renormalization')
# plt.legend(prop={'size': 7})
# plt.tight_layout()
# plt.show()




def function_SNIa_DTD(delay_time):
    # The delay time distribution (DTD) in the unit of per year per total stellar mass [solar]
    # DTD for SNIa is adopted from Maoz & Mannucci 2012, 29, 447–465, their equation 13
    # with a consistent assumed IMF – the Bell et al. 2003 diet-Salpeter IMF
    # gaptime = 6 * 10 ** 7
    gaptime = 4 * 10 ** 7
    powerindex = -1
    if delay_time < gaptime:  # [yr] #  2.3 * 10 ** 7 for a burst of star formation from Greggio 1983
        number = 0
    else:
        number = 4 * 10 ** (-4) * delay_time ** (powerindex) / (gaptime)**(powerindex) * (gaptime)**(-1)  # DTD of Maoz2012
        # Normalized such that the DTD integral over 10Gyr for diet-Salpeter is N_SN/M_sun = 2 * 10^-3 (M_sun^-1)
        # number = 0.16288551211 * 10 ** (-4) / 10**(delay_time/(5*10**9)) / 10**(gaptime/(5*10**9)) * (gaptime)**(-1)  # DTD of Lacchin19
        # Normalized such that the DTD integral over 10Gyr for diet-Salpeter is N_SN/M_sun = 2 * 10^-3 (M_sun^-1)

        # number = 0.50986652089 * 10 ** (-4) / 10**(delay_time/(5*10**9)) / 10**(gaptime/(5*10**9)) * (gaptime)**(-1)  # DTD of Lacchin19
        # Normalized such that the SNIa rate at 10 Gyr is the same as the Maoz2012's DTD.

        # Normalized such that the DTD integral over 10Gyr for diet-Salpeter is N_SN/M_sun = 2 * 10^-3 (M_sun^-1)
        # This value changes with igimf where top-heavy and bottom-heavy IGIMF will have lower number of SNIa
        # as the number of stars within the mass range 3.0001 to 8 solar mass is smaller.
        # The observational uncertainty being +-50%. See Maoz & Mannucci 2012 their Table 1
    return number


def function_read_lifetime(str_yield_table, Z_select_in_table):
    #### if apply instantaneous recycling approximation ####
    global instantaneous_recycling
    if instantaneous_recycling == True:
        mass_1 = 0
        mass = [0.08, 1, 1.00001, 150]
        lifetime_table = [1e12, 1e12, 0.1, 0.1]
    elif Z_select_in_table[0] == 'out':
        file_lifetime = open(
            'yield_tables/rearranged/setllar_lifetime_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table[1]),
            'r')
        data = file_lifetime.readlines()
        metallicity = data[1]
        mass_1 = data[3]
        lifetime_ = data[5]
        file_lifetime.close()
        mass = [float(x) for x in mass_1.split()]
        lifetime_table = [float(x) for x in lifetime_.split()]
    else:
        file_lifetime_low = open(
            'yield_tables/rearranged/setllar_lifetime_from_portinari98/portinari98_Z={}.txt'.format(
                Z_select_in_table[1]),
            'r')
        data_low = file_lifetime_low.readlines()
        metallicity = data_low[1]
        mass_1 = data_low[3]
        lifetime_low = data_low[5]
        file_lifetime_low.close()

        file_lifetime_high = open(
            'yield_tables/rearranged/setllar_lifetime_from_portinari98/portinari98_Z={}.txt'.format(
                Z_select_in_table[3]),
            'r')
        data_high = file_lifetime_high.readlines()
        lifetime_high = data_high[5]
        file_lifetime_high.close()
        mass = [float(x) for x in mass_1.split()]
        lifetime_table_low = [float(x) for x in lifetime_low.split()]
        lifetime_table_high = [float(x) for x in lifetime_high.split()]
        x1 = Z_select_in_table[1]
        x2 = Z_select_in_table[2]
        x3 = Z_select_in_table[3]
        lifetime_table = [y1+(y3-y1)*(x2-x1)/(x3-x1) for y1, y3 in zip(lifetime_table_low, lifetime_table_high)]
    return (mass_1, mass, lifetime_table)


def function_read_Mfinal(str_yield_table, Z_select_in_table):
    if Z_select_in_table[0] == 'out':
        file_final_mass = open(
            "yield_tables/rearranged___/setllar_final_mass_from_portinari98/portinari98_Z={}.txt".format(
                Z_select_in_table[1]),
            'r')
        data = file_final_mass.readlines()
        metallicity2 = data[1]
        mass_2 = data[3]
        Mfinal_ = data[5]
        file_final_mass.close()
        Mfinal_table = [float(x) for x in Mfinal_.split()]
    else:
        file_final_mass = open(
            "yield_tables/rearranged___/setllar_final_mass_from_portinari98/portinari98_Z={}.txt".format(
                Z_select_in_table[1]),
            'r')
        data = file_final_mass.readlines()
        metallicity2 = data[1]
        mass_2 = data[3]
        Mfinal_low = data[5]
        file_final_mass.close()
        Mfinal_table_low = [float(x) for x in Mfinal_low.split()]
        file_final_mass = open(
            "yield_tables/rearranged___/setllar_final_mass_from_portinari98/portinari98_Z={}.txt".format(
                Z_select_in_table[3]),
            'r')
        data = file_final_mass.readlines()
        Mfinal_high = data[5]
        file_final_mass.close()
        Mfinal_table_high = [float(x) for x in Mfinal_high.split()]
        x1 = Z_select_in_table[1]
        x2 = Z_select_in_table[2]
        x3 = Z_select_in_table[3]
        Mfinal_table = [y1+(y3-y1)*(x2-x1)/(x3-x1) for y1, y3 in zip(Mfinal_table_low, Mfinal_table_high)]
    return (mass_2, Mfinal_table)


def lindexsplit(List, *lindex):
    index = list(lindex)
    index.sort()
    templist1 = []
    templist2 = []
    templist3 = []
    breakcounter = 0
    itemcounter = 0
    finalcounter = 0
    numberofbreaks = len(index)
    totalitems = len(List)
    lastindexval = index[(len(index) - 1)]
    finalcounttrigger = (totalitems - (lastindexval + 1))
    for item in List:
        itemcounter += 1
        indexofitem = itemcounter - 1
        nextbreakindex = index[breakcounter]
        # Less than the last cut
        if breakcounter <= numberofbreaks:
            if indexofitem < nextbreakindex:
                templist1.append(item)
            elif breakcounter < (numberofbreaks - 1):
                templist1.append(item)
                templist2.append(templist1)
                templist1 = []
                breakcounter += 1
            else:
                if indexofitem <= lastindexval and indexofitem <= totalitems:
                    templist1.append(item)
                    templist2.append(templist1)
                    templist1 = []
                else:
                    if indexofitem >= lastindexval and indexofitem < totalitems + 1:
                        finalcounter += 1
                        templist3.append(item)
                        if finalcounter == finalcounttrigger:
                            templist2.append(templist3)
    return templist2


def function_read_Mmetal(str_yield_table, Z_select_in_table_2, Z_select_in_table_3):
    global mm, zz
    if str_yield_table == "Kobayashi06" or str_yield_table == "portinari98":
        if Z_select_in_table_2[0] == 'out':
            file_Metal_eject = open(
                'yield_tables/rearranged___/setllar_Metal_eject_mass_from_{}/{}_Z={}.txt'.format(str_yield_table,
                                                                                                 str_yield_table,
                                                                                                 Z_select_in_table_2[1]),
                'r')
            data = file_Metal_eject.readlines()
            metallicity = data[1]
            mass_2 = data[3]
            Metal_eject_ = data[5]
            file_Metal_eject.close()
            mass = [float(x) for x in mass_2.split()]
            Metal_eject_table = [float(x) for x in Metal_eject_.split()]
        else:
            file_Metal_eject = open(
                'yield_tables/rearranged___/setllar_Metal_eject_mass_from_{}/{}_Z={}.txt'.format(str_yield_table,
                                                                                                 str_yield_table,
                                                                                                 Z_select_in_table_2[
                                                                                                     1]),
                'r')
            data = file_Metal_eject.readlines()
            metallicity = data[1]
            mass_2 = data[3]
            Metal_eject_low = data[5]
            file_Metal_eject.close()
            mass = [float(x) for x in mass_2.split()]
            Metal_eject_table_low = [float(x) for x in Metal_eject_low.split()]
            file_Metal_eject = open(
                'yield_tables/rearranged___/setllar_Metal_eject_mass_from_{}/{}_Z={}.txt'.format(str_yield_table,
                                                                                                 str_yield_table,
                                                                                                 Z_select_in_table_2[
                                                                                                     3]),
                'r')
            data = file_Metal_eject.readlines()
            Metal_eject_high = data[5]
            file_Metal_eject.close()
            Metal_eject_table_high = [float(x) for x in Metal_eject_high.split()]
            x1 = Z_select_in_table_2[1]
            x2 = Z_select_in_table_2[2]
            x3 = Z_select_in_table_2[3]
            Metal_eject_table = [y1 + (y3 - y1) * (x2 - x1) / (x3 - x1) for y1, y3 in
                            zip(Metal_eject_table_low, Metal_eject_table_high)]

    elif str_yield_table == "WW95":
        if Z_select_in_table_2[2] < (Z_select_in_table_2[1]+Z_select_in_table_2[3])/2:
            Z_select_in_table_2 = Z_select_in_table_2[1]
        else:
            Z_select_in_table_2 = Z_select_in_table_2[3]
        if Z_select_in_table_3[2] < (Z_select_in_table_3[1]+Z_select_in_table_3[3])/2:
            Z_select_in_table_3 = Z_select_in_table_3[1]
        else:
            Z_select_in_table_3 = Z_select_in_table_3[3]

        file_Metal_eject = open(
            'yield_tables/rearranged___/setllar_Metal_eject_mass_from_{}/{}_Z={}.txt'.format(str_yield_table,
                                                                                             str_yield_table,
                                                                                             Z_select_in_table_2),
            'r')
        data = file_Metal_eject.readlines()
        mass_2 = data[3]
        Metal_eject_ = data[5]
        file_Metal_eject.close()
        mass = [float(x) for x in mass_2.split()]
        mass = lindexsplit(mass, 153)[1]
        Metal_eject_table = [float(x) for x in Metal_eject_.split()]
        Metal_eject_table = lindexsplit(Metal_eject_table, 153)[1]

        file_Metal_eject = open(
            'yield_tables/rearranged___/setllar_Metal_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                Z_select_in_table_3),
            'r')
        data = file_Metal_eject.readlines()
        mass_2 = data[3]
        Metal_eject_ = data[5]
        file_Metal_eject.close()
        mass_agb = [float(x) for x in mass_2.split()]
        mass_agb = lindexsplit(mass_agb, 153)[0]
        Metal_eject_table_agb = [float(x) for x in Metal_eject_.split()]
        Metal_eject_table_agb = lindexsplit(Metal_eject_table_agb, 153)[0]

        mass = mass_agb + mass
        Metal_eject_table = Metal_eject_table_agb + Metal_eject_table
    else:
        print('Input str_yield_table does not exist.')
    return (mass_2, mass, Metal_eject_table)


def function_read_M_element(element, str_yield_table, Z_select_in_table_2, Z_select_in_table_3):
    if str_yield_table == "portinari98" or str_yield_table == "Kobayashi06":
        if element == "H" or element == "He" or element == "C" or element == "N" or element == "O" or element == "Mg"\
                or element == "Ne" or element == "Si" or element == "S" or element == "Ca" or element == "Fe":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(element, str_yield_table, str_yield_table,
                    Z_select_in_table_2[1]),
                'r')
            data = file_M_eject.readlines()
            M_eject_low = data[5]
            file_M_eject.close()
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(element, str_yield_table, str_yield_table,
                    Z_select_in_table_2[3]),
                'r')
            data = file_M_eject.readlines()
            M_eject_high = data[5]
            file_M_eject.close()
        else:
            print("Error: element parameter for function_read_M_element do not exsit.")
        M_eject_table_low = [float(x) for x in M_eject_low.split()]
        M_eject_table_high = [float(x) for x in M_eject_high.split()]
        x1 = Z_select_in_table_2[1]
        x2 = Z_select_in_table_2[2]
        x3 = Z_select_in_table_2[3]
        if x3 == x1:
            M_eject_table = M_eject_table_high
        else:
            M_eject_table = [y1 + (y3 - y1) * (x2 - x1) / (x3 - x1) for y1, y3 in
                             zip(M_eject_table_low, M_eject_table_high)]
    elif str_yield_table == "WW95":
        if Z_select_in_table_2[2] < (Z_select_in_table_2[1]+Z_select_in_table_2[3])/2:
            Z_select_in_table_2 = Z_select_in_table_2[1]
        else:
            Z_select_in_table_2 = Z_select_in_table_2[3]
        if Z_select_in_table_3[2] < (Z_select_in_table_3[1]+Z_select_in_table_3[3])/2:
            Z_select_in_table_3 = Z_select_in_table_3[1]
        else:
            Z_select_in_table_3 = Z_select_in_table_3[3]

        if element == "H":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_H_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "He":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_He_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "C":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_C_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "N":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_N_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "O":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_O_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Mg":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_Mg_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Ne":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_Ne_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Si":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_Si_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "S":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_S_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Ca":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_Ca_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Fe":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_Fe_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        else:
            file_M_eject = 0
            print("Error: element parameter for function_read_M_element do not exsit.")
        data = file_M_eject.readlines()
        M_eject_ = data[5]
        file_M_eject.close()
        M_eject_table = [float(x) for x in M_eject_.split()]
        M_eject_table = lindexsplit(M_eject_table, 153)[1]

        if element == "H":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_H_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                    Z_select_in_table_3),
                'r')
        elif element == "He":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_He_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                    Z_select_in_table_3),
                'r')
        elif element == "C":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_C_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                    Z_select_in_table_3),
                'r')
        elif element == "N":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_N_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                    Z_select_in_table_3),
                'r')
        elif element == "O":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_O_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                    Z_select_in_table_3),
                'r')
        elif element == "Mg":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_Mg_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                    Z_select_in_table_3),
                'r')
        elif element == "Ne":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_Ne_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                    Z_select_in_table_3),
                'r')
        elif element == "Si":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_Si_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                    Z_select_in_table_3),
                'r')
        elif element == "S":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_S_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                    Z_select_in_table_3),
                'r')
        elif element == "Ca":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_Ca_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                    Z_select_in_table_3),
                'r')
        elif element == "Fe":
            file_M_eject = open(
                'yield_tables/rearranged___/setllar_Fe_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(
                    Z_select_in_table_3),
                'r')
        else:
            file_M_eject = 0
            print("Error: element parameter for function_read_M_element do not exsit.")
        data = file_M_eject.readlines()
        M_eject_ = data[5]
        file_M_eject.close()
        M_eject_table_agb = [float(x) for x in M_eject_.split()]
        M_eject_table_agb = lindexsplit(M_eject_table_agb, 153)[0]

        M_eject_table = M_eject_table_agb + M_eject_table
    return M_eject_table


def get_BH_mass(mass_boundary, mass_grid_table_number, Mtarget_table_number, mass_calibration_factor,
                steller_mass_upper_bound):
    if mass_boundary < steller_mass_upper_bound:
        BH_mass = function_get_target_mass_in_range(max(mass_boundary, 40), steller_mass_upper_bound,
                                                    mass_grid_table_number,
                                                    Mtarget_table_number, mass_calibration_factor)
    else:
        BH_mass = 0
    return BH_mass


def get_NS_mass(mass_boundary, mass_grid_table_number, Mtarget_table_number, mass_calibration_factor):
    if mass_boundary < 40:
        NS_mass = function_get_target_mass_in_range(max(mass_boundary, 8), 40, mass_grid_table_number,
                                                    Mtarget_table_number, mass_calibration_factor)
    else:
        NS_mass = 0
    return NS_mass


def get_WD_mass(mass_boundary, mass_grid_table_number, Mtarget_table_number, mass_calibration_factor):
    if mass_boundary < 8:
        WD_mass = function_get_target_mass_in_range(max(mass_boundary, 0.08), 8, mass_grid_table_number,
                                                    Mtarget_table_number, mass_calibration_factor)
    else:
        WD_mass = 0
    return WD_mass


def function_get_target_mass_in_range(lower_mass_limit, upper_mass_limit, mass_grid_table_number, Mtarget_table_number,
                                      mass_calibration_factor):
    integrate_in_range = quad(integrator_for_function_get_target_mass_in_range, lower_mass_limit, upper_mass_limit,
                              (mass_grid_table_number, Mtarget_table_number), limit=40)[
        0]  ####################################
    target_mass_in_range = mass_calibration_factor * integrate_in_range
    return target_mass_in_range


def integrator_for_function_get_target_mass_in_range(initial_mass, mass_grid_table_number, Mtarget_table_number):
    global igimf_mass_function
    mass = igimf_mass_function(initial_mass)
    mass_fraction = function_get_target_mass(initial_mass, mass_grid_table_number, Mtarget_table_number) / initial_mass
    integrator = mass * mass_fraction
    return integrator


def function_get_target_mass(initial_mass, mass_grid_table_number, Mtarget_table_number):
    global mass_grid_table, mass_grid_table2, Mfinal_table, Mmetal_table, M_element_table
    if Mtarget_table_number == 1:
        Mtarget_table = Mfinal_table
    if Mtarget_table_number == 2:
        Mtarget_table = Mmetal_table
    if Mtarget_table_number == "H":
        Mtarget_table = M_element_table[0]
    if Mtarget_table_number == "He":
        Mtarget_table = M_element_table[1]
    if Mtarget_table_number == "C":
        Mtarget_table = M_element_table[2]
    if Mtarget_table_number == "N":
        Mtarget_table = M_element_table[3]
    if Mtarget_table_number == "O":
        Mtarget_table = M_element_table[4]
    if Mtarget_table_number == "Mg":
        Mtarget_table = M_element_table[5]
    if Mtarget_table_number == "Ne":
        Mtarget_table = M_element_table[6]
    if Mtarget_table_number == "Si":
        Mtarget_table = M_element_table[7]
    if Mtarget_table_number == "S":
        Mtarget_table = M_element_table[8]
    if Mtarget_table_number == "Ca":
        Mtarget_table = M_element_table[9]
    if Mtarget_table_number == "Fe":
        Mtarget_table = M_element_table[10]
    if mass_grid_table_number == 1:
        mass_grid_table_n = mass_grid_table
    if mass_grid_table_number == 2:
        mass_grid_table_n = mass_grid_table2
    if initial_mass < mass_grid_table_n[0] or initial_mass > mass_grid_table_n[-1]:
        print('Warning: function_get_remnant_mass initial_mass out of range')
        print("initial_mass=", initial_mass, "< mass grid lower boundary =", mass_grid_table_n[0])
    length_list_mass = len(mass_grid_table_n)
    x = round(length_list_mass / 2)
    i = 0
    low = 0
    high = length_list_mass
    if initial_mass == mass_grid_table_n[0]:
        x = 0
    elif initial_mass == mass_grid_table_n[-1]:
        x = -1
    else:
        while i < math.ceil(math.log(length_list_mass, 2)):
            if initial_mass == mass_grid_table_n[x]:
                break
            elif initial_mass > mass_grid_table_n[x]:
                low = x
                x = x + round((high - x) / 2)
            else:
                high = x
                x = x - round((x - low) / 2)
            (i) = (i + 1)
    if mass_grid_table_n[x - 1] < initial_mass < mass_grid_table_n[x]:
        x = x - 1
    target_mass = round(
        (Mtarget_table[x] + (Mtarget_table[x + 1] - Mtarget_table[x]) * (initial_mass - mass_grid_table_n[x]) /
         (mass_grid_table_n[x + 1] - mass_grid_table_n[x])), 5)
    return target_mass


    # ### Define initial stellar mass boundary for WD, NS, and BH.
    # mass_boundary_WD_to_NS = 8  # [solar mass]
    # mass_boundary_NS_to_BH = 40  # [solar mass]
    #
    # # Define the observational sensitive mass range for galaxy total mass estimation
    # mass_boundary_observe = [mass_boundary_observe_low, mass_boundary_observe_up]


    # ### Calculate total mass at each time ###
    # M_tot = 0
    # M_tot_time_list = []
    # new_time = 1
    # M_tot_list = []
    # for SFH in SFH_input:
    #     formed_mass = SFH * 10 ** 7
    #     M_tot += formed_mass
    #     M_tot_time_list += [new_time]
    #     if M_tot == 0:
    #         M_tot_list += [1, 1]
    #     else:
    #         M_tot_list += [M_tot, M_tot]
    #     new_time += 10 ** 7
    #     M_tot_time_list += [new_time]
    #
    # Log_M_tot = math.log(M_tot, 10)
    # M_tot_time_list += [time_axis[-1]]
    # M_tot_list += [M_tot_list[-1]]
    #
    #
    # ### Calculate the observational estimated total mass of the galaxy ###
    # # Assuming the estimation done everything right, e.g., stellar evolution module, SFH, dust extinction, metallicity,
    # # excepet assumed an universal Kroupa IMF that is not what really happend
    # # (although this assumption contradict itself because it is impossible to get everything else right with a wrong IMF).
    # # We using the stellar population with mass in 0.08 - 3 solar mass to estimate the total stellar mass with Kroupa IMF
    # # and compare it with the real total mass
    #
    # imf_file_name = "{}_IMF".format(IMF_name)
    #
    # # estimated total mass with Kroupa IMF =
    # M_tot_est_list = []
    # IMF = __import__(imf_file_name)
    # a = quad(IMF.imf, 0.08, steller_mass_upper_bound, limit=50)[0]
    # b = quad(IMF.imf, mass_boundary_observe[0], mass_boundary_observe[1], limit=40)[0]
    # for mass_in_range in M_in_range_list:
    #     est_mass = mass_in_range * a / b
    #     if est_mass == 0:
    #         M_tot_est_list += [1]
    #     else:
    #         M_tot_est_list += [est_mass]


def function_get_igimf_for_this_epoch(SFR_input, Z_over_X, this_time, this_epoch, check_igimf):
    # this function calculate igimf, write them in directory Generated_IGIMFs, and import the file
    # with igimf = function_get_igimf_for_every_epoch(SFH_input, Z, Z_solar),
    # the igimf can be called by: igimf.custom_imf(stellar_mass, this_time).
    function_generate_igimf_file(SFR=SFR_input, Z_over_X=Z_over_X, printout=None, sf_epoch=this_epoch,
                                 check_igimf=check_igimf)
    if SFR_input == 0:
        igimf_file_name = "igimf_SFR_Zero"
    else:
        igimf_file_name = "igimf_SFR_{}_Fe_over_H_{}".format(round(math.log(SFR_input, 10) * 100000),
                                                             round(Z_over_X * 100000))
    igimf = __import__(igimf_file_name)
    # import os
    # if os.path.isfile('Generated_IGIMFs/' + igimf_file_name + '.py'):
    #     igimf = __import__(igimf_file_name)
    # else:
    #     cwd = os.getcwd()
    #     igimf = __import__(cwd + '/galIMF/Generated_IGIMFs/' + igimf_file_name)

    # if shows ModuleNotFoundError:
    # No module named 'igimf_SFR_..._Fe_over_H_...',
    # then try clear all (except for the first and last) lines in the file Generated_IGIMFs/all_igimf_list.txt.
    # This will force the program to generate new IGIMF functions for future use,
    # instead of looking for the IGIMF in the old generated ones.
    return igimf


def function_generate_igimf_file(SFR=None, Z_over_X=None, printout=False, sf_epoch=0, check_igimf=False):
    # This funtion check if the parameter for generating a new IGIMF match an old one,
    # if not, the function generate a new IGIMF and add it to the generated-IGIMF list.

    # --------------------------------------------------------------------------------------------------------------------------------
    # import modules and libraries
    # --------------------------------------------------------------------------------------------------------------------------------

    import galimf  # galIMF containing IGIMF function and OSGIMF function and additional computational modules
    import numpy as np
    import math
    import time
    import os

    Generated_IGIMFs_path = 'Generated_IGIMFs'
    if os.path.isdir(Generated_IGIMFs_path) == False:
        Generated_IGIMFs_path = '/galIMF/Generated_IGIMFs'
        if os.path.isdir(Generated_IGIMFs_path) == False:
            cwd = os.getcwd()
            Generated_IGIMFs_path = cwd + '/galIMF/Generated_IGIMFs'
    file_name = '/igimf_SFR_{}_Fe_over_H_{}.py'.format(round(math.log(SFR, 10) * 100000),
                                                       round(Z_over_X * 100000))
    file_path_and_name = Generated_IGIMFs_path + file_name

    # --------------------------------------------------------------------------------------------------------------------------------
    # check if the required IGIMF has already been generated
    # --------------------------------------------------------------------------------------------------------------------------------

    exist = 0

    if check_igimf == True:

        if os.path.isfile(file_path_and_name):
            igimf_file_name = "igimf_SFR_{}_Fe_over_H_{}".format(round(math.log(SFR, 10) * 100000),
                                                                 round(Z_over_X * 100000))
            igimf_____ = __import__(igimf_file_name)
            if hasattr(igimf_____, "custom_imf"):
                # print("find IGIMF file '{}' for a galaxy with [Z/X]={}, SFR={}".format(file_path_and_name, round(Z_over_X, 2), SFR))
                exist = 1
        # else:
        #     print("{} is not a file".format(file_path_and_name))

        # check_file = open(Generated_IGIMFs_path + '/all_igimf_list.txt', 'r')
        # igimf_list_line = check_file.readlines()
        # check_file.close()
        # i = 0
        # while i < len(igimf_list_line):
        #     data = [float(a_block) for a_block in igimf_list_line[i].split()]
        #     if SFR == data[0] and Z_over_X == data[1]:
        #         exist = 1
        #         break
        #     (i) = (i + 1)

    if exist == 0 and SFR != 0:
        # print("Generating new IGIMF file '{}' for a galaxy with [Z/X]={}, SFR={}".format(file_path_and_name, Z_over_X, SFR))

        # # --------------------------------------------------------------------------------------------------------------------------------
        # # add new headline into the list file -- all_igimf_list.txt:
        # # --------------------------------------------------------------------------------------------------------------------------------
        #
        # check_file = open('Generated_IGIMFs/all_igimf_list.txt', 'r')
        # igimf_list = check_file.read()
        # check_file.close()
        #
        # check_file = open('Generated_IGIMFs/all_igimf_list.txt', 'w')
        # new_headline = igimf_list + '{} {}\n'.format(SFR, Z_over_X)
        # check_file.write(new_headline)
        # check_file.close()

        # --------------------------------------------------------------------------------------------------------------------------------
        # Define code parameters necesarry for the computations:
        # --------------------------------------------------------------------------------------------------------------------------------

        # the most crutial ones, which you most likely might want to change

        if SFR is None:
            SFR = float(
                input(
                    "Please input the galaxy-wide SFR in solar mass per year and ended the input with the return key. "
                    "(A typical input SFR is from 0.0001 to 10000. "
                    "We recommed a value smallar than 0.01 for the first run as high SFR calculations take more time.)\n"
                    "You can input 1e-4 as 0.0001\n"
                    "\nSFR [Msolar/yr] = "))
            # Star Formation Rate [solar mass / yr]
        if SFR != 0:
            bindw = galimf.resolution_histogram_relative = 10 ** (max((0 - math.log(SFR, 10)), 0) ** (0.2) - 1.9)
        # will change the resolution of histogram for optimall sampling automatically addjusted with SFR value.

        gwIMF_model = "IGIMF_Z"

        if gwIMF_model == "IGIMF3":
            alpha3_model = 2  # 1  # IMF high-mass-end power-index model, see Function_alpha_3_change in file 'galimf.py'
            alpha_2 = 2.3  # IMF middle-mass power-index
            alpha_1 = 1.3  # IMF low-mass-end power-index
            alpha2_model = 1  # 1  # see file 'galimf.py'
            alpha1_model = 1  # 0 # see file 'galimf.py'
            beta_model = 1
            R14orNOT = False
        elif gwIMF_model == "IGIMF_Z":
            alpha3_model = 2  # 1  # IMF high-mass-end power-index model, see Function_alpha_3_change in file 'galimf.py'
            alpha_2 = 2.3  # IMF middle-mass power-index
            alpha_1 = 1.3  # IMF low-mass-end power-index
            alpha2_model = 'Z'  # 1  # see file 'galimf.py'
            alpha1_model = 'Z'  # 0 # see file 'galimf.py'
            beta_model = 1
            R14orNOT = False
        elif gwIMF_model == "IGIMF2d5":
            alpha3_model = 2  # 1  # IMF high-mass-end power-index model, see Function_alpha_3_change in file 'galimf.py'
            alpha_2 = 2.3  # IMF middle-mass power-index
            alpha_1 = 1.3  # IMF low-mass-end power-index
            alpha2_model = 'IGIMF2.5'  # 1  # see file 'galimf.py'
            alpha1_model = 'IGIMF2.5' # 0 # see file 'galimf.py'
            beta_model = 1
            R14orNOT = False
        elif gwIMF_model == "IGIMF2":
            alpha3_model = 2  # 1  # IMF high-mass-end power-index model, see Function_alpha_3_change in file 'galimf.py'
            alpha_2 = 2.3  # IMF middle-mass power-index
            alpha_1 = 1.3  # IMF low-mass-end power-index
            alpha2_model = 0  # 1  # see file 'galimf.py'
            alpha1_model = 0  # 0 # see file 'galimf.py'
            beta_model = 1
            R14orNOT = False
        elif gwIMF_model == "IGIMF_R14":
            alpha3_model = 'R14' # 'R14'  # 2  # 1  # IMF high-mass-end power-index model, see Function_alpha_3_change in file 'galimf.py'
            alpha_2 = 2.3  # IMF middle-mass power-index
            alpha_1 = 1.3  # IMF low-mass-end power-index
            alpha2_model = 'R14' # 'R14'  # 1  # see file 'galimf.py'
            alpha1_model = 0 # 0 # see file 'galimf.py'
            beta_model = 0
            R14orNOT = True

        # ----------------------------------------------------------------

        # Parameters below are internal parameters of the theory.
        # Read Yan et al. 2017 carefully before change them!

        delta_t = 10.  # star formation epoch [Myr]
        I_ecl = 1.  # normalization factor in the Optimal Sampling condition equation
        M_ecl_U = 10 ** 9  # 10**(0.75 * math.log(SFR, 10) + 4.8269) # Recchi 2009
        # 10 ** 15  # embedded cluster mass upper limit [solar mass]
        M_ecl_L = 5.  # embedded cluster mass lower limit [solar mass]
        I_str = 1.  # normalization factor in the Optimal Sampling condition equation
        M_str_L = 0.08  # star mass lower limit [solar mass]
        M_turn = 0.5  # IMF power-index breaking mass [solar mass]
        M_turn2 = 1.  # IMF power-index breaking mass [solar mass]
        M_str_U = 150  # star mass upper limit [solar mass]

        if printout == True:
            print("\n - GalIMF run in progress..")
        start_time = time.time()

        # --------------------------------------------------------------------------------------------------------------------------------
        # Construct IGIMF:
        # --------------------------------------------------------------------------------------------------------------------------------

        if printout == True:
            print("\nCalculating IGIMF......")

        galimf.function_galimf(
            "I",  # IorS ### "I" for IGIMF; "OS" for OSGIMF
            R14orNOT, # True or False
            SFR,  # Star Formation Rate [solar mass / yr]
            alpha3_model,  # IMF high-mass-end power-index model, see file 'alpha3.py'
            delta_t,  # star formation epoch [Myr]
            Z_over_X,  # M_over_H
            I_ecl,  # normalization factor in the Optimal Sampling condition equation
            M_ecl_U,  # embedded cluster mass upper limit [solar mass]
            M_ecl_L,  # embedded cluster mass lower limit [solar mass]
            beta_model,  ### ECMF power-index model, see file 'beta.py'
            I_str,  # normalization factor in the Optimal Sampling condition equation
            M_str_L,  # star mass lower limit [solar mass]
            alpha_1,  # IMF low-mass-end power-index
            alpha1_model,  # see file 'alpha1.py'
            M_turn,  # IMF power-index change point [solar mass]
            alpha_2,  # IMF middle-mass power-index
            alpha2_model,  # see file 'alpha2.py'
            M_turn2,  # IMF power-index change point [solar mass]
            M_str_U,  # star mass upper limit [solar mass]
            printout
        )

        if printout == True:
            ### Decorate the output file ###
            igimf_raw = np.loadtxt('GalIMF_IGIMF.txt')
            if M_str_U - igimf_raw[-1][0] > 0.01:
                file = open('GalIMF_IGIMF.txt', 'a')
                file.write("{} 0\n\n".format(igimf_raw[-1][0] + 0.01))
                file.write("{} 0".format(M_str_U))
                file.close()
            else:
                file = open('GalIMF_IGIMF.txt', 'a')
                file.write("{} 0".format(M_str_U))
                file.close()

        global masses, igimf

        masses = np.array(galimf.List_M_str_for_xi_str)
        igimf = np.array(galimf.List_xi)

        #######################################################
        # generated igimf is normalized by default to a total mass formed in 10 Myr given the SFR,
        # i.e., total stellar mass.
        # to change the normalization, uncomment the below commented part:
        #######################################################
        # Norm = simps(igimf*masses,masses) #- normalization to a total mass
        # Norm = simps(igimf,masses) #- normalization to number of stars
        # Mtot1Myr = SFR*10*1.e6 #total mass formed in 10 Myr
        # igimf = np.array(igimf)*Mtot1Myr/Norm
        #######################################################


        global length_of_igimf
        length_of_igimf = len(igimf)

        def write_imf_input2():
            global file, masses, igimf
            if SFR == 0:
                file = open('Generated_IGIMFs/igimf_SFR_Zero.py', 'w')
                file.write("def custom_imf(mass, time):  # there is no time dependence for IGIMF\n")
                file.write("    return 0\n")
                file.close()
            else:
                file = open('Generated_IGIMFs/igimf_SFR_{}_Fe_over_H_{}.py'.format(round(math.log(SFR, 10) * 100000),
                                                                                   round(Z_over_X * 100000)), 'w')
                file.write("# File to define a custom IMF\n"
                           "# The return value represents the chosen IMF value for the input mass\n\n\n")
                file.write("def custom_imf(mass, time):  # there is no time dependence for IGIMF\n")
                file.write("    if mass < 0.08:\n")
                file.write("        return 0\n")
                file.write("    elif mass < %s:\n" % masses[1])
                if masses[0] - masses[1] == 0:
                    k = 0
                    b = 0
                else:
                    k = (igimf[0] - igimf[1]) / (masses[0] - masses[1])
                    b = igimf[0] - k * masses[0]
                file.write("        return {} * mass + {}\n".format(k, b))
                write_imf_input_middle2(1)
                file.write("    else:\n")
                file.write("        return 0\n")
                file.close()
            return

        def write_imf_input_middle2(i):
            global file, length_of_igimf
            while i < length_of_igimf - 1:
                file.write("    elif mass < %s:\n" % masses[i + 1])
                if masses[i] - masses[i + 1] == 0:
                    k = 0
                    b = 0
                else:
                    k = (igimf[i] - igimf[i + 1]) / (masses[i] - masses[i + 1])
                    b = igimf[i] - k * masses[i]
                file.write("        return {} * mass + {}\n".format(k, b))
                (i) = (i + 3)
            return

        write_imf_input2()

        def write_imf_input3():
            global file, masses, igimf
            if SFR == 0:
                file = open('Generated_IGIMFs/igimf_epoch_{}.py'.format(sf_epoch), 'w')
                file.write("def custom_imf(mass, time):  # there is no time dependence for IGIMF\n")
                file.write("    return 0\n")
                file.close()
            else:
                file = open('Generated_IGIMFs/igimf_epoch_{}.py'.format(sf_epoch), 'w')
                file.write("# File to define a custom IMF\n"
                           "# The return value represents the chosen IMF value for the input mass\n\n\n")
                file.write("def custom_imf(mass, time):  # there is no time dependence for IGIMF\n")
                file.write("    if mass < 0.08:\n")
                file.write("        return 0\n")
                file.write("    elif mass < %s:\n" % masses[1])
                if masses[0] - masses[1] == 0:
                    k = 0
                    b = 0
                else:
                    k = (igimf[0] - igimf[1]) / (masses[0] - masses[1])
                    b = igimf[0] - k * masses[0]
                file.write("        return {} * mass + {}\n".format(k, b))
                write_imf_input_middle2(1)
                file.write("    else:\n")
                file.write("        return 0\n")
                file.close()
            return

        write_imf_input3()

        if printout == True:
            print("imf_input.py rewritten for SFR = {} and metallicity = {}\n".format(SFR, Z_over_X))

            file = open('../gimf_Fe_over_H.txt', 'w')
            file.write("{}".format(Z_over_X))
            file.close()

            file = open('../gimf_SFR.txt', 'w')
            file.write("{}".format(SFR))
            file.close()

            print(" - GalIMF run completed - Run time: %ss -\n\n" % round((time.time() - start_time), 2))
        return
    return


def function_element_abundunce(solar_abu_table, element_1_name, element_2_name, metal_1_mass, metal_2_mass, instant_ejection):
    # this function calculate the atom number ratio compare to solar value [metal/H]

    # The following warning might be due to too large a timestep.
    # Try applying the "high_time_resolution=True"
    # but the simulation will take longer time.
    if metal_2_mass == 0:
        if metal_1_mass == 0:
            metal_1_over_2 = 0
        elif metal_1_mass > 0:
            metal_1_over_2 = 6
        elif metal_1_mass < 0:
            if instant_ejection==False:
                print("Warning: current {} mass < 0. See galevo.py".format(element_1_name))
            metal_1_over_2 = -6
    elif metal_2_mass < 0:
        if instant_ejection == False:
            print("Warning: current {} mass < 0. See galevo.py".format(element_2_name))
        if metal_1_mass == 0:
            metal_1_over_2 = 6
        elif metal_1_mass > 0:
            metal_1_over_2 = 6
        elif metal_1_mass < 0:
            if instant_ejection == False:
                print("Warning: current {} mass < 0. See galevo.py".format(element_1_name))
            metal_1_over_2 = 0
    else:
        if metal_1_mass == 0:
            metal_1_over_2 = -6
        elif metal_1_mass < 0:
            if instant_ejection == False:
                print("Warning: current {} mass < 0. See galevo.py".format(element_1_name))
            metal_1_over_2 = -6
        else:
            solar_metal_1_logarithmic_abundances = element_abundances_solar.function_solar_element_abundances(
                solar_abu_table, element_1_name)
            solar_metal_2_logarithmic_abundances = element_abundances_solar.function_solar_element_abundances(
                solar_abu_table, element_2_name)
            metal_1_element_weight = element_weight_table.function_element_weight(element_1_name)
            metal_2_element_weight = element_weight_table.function_element_weight(element_2_name)

            metal_1_over_2 = math.log(metal_1_mass / metal_2_mass / metal_1_element_weight * metal_2_element_weight, 10) \
                             - (solar_metal_1_logarithmic_abundances - solar_metal_2_logarithmic_abundances)
    return metal_1_over_2


def function_get_avaliable_Z(str_yield_table):
    # extract avalible metallicity in the given grid table
    # stellar life-time table and metal production tables have different avalible metal grid.
    import os

    yield_path = 'yield_tables'
    if os.path.isdir(yield_path) == False:
        yield_path = '/galIMF/yield_tables'
        if os.path.isdir(yield_path) == False:
            cwd = os.getcwd()
            yield_path = cwd + '/galIMF/yield_tables'

    # list 1
    file_names_setllar_lifetime_from_str_yield_table = os.listdir(
        yield_path + '/rearranged___/setllar_lifetime_from_portinari98')
    Z_table_list = []
    for name in file_names_setllar_lifetime_from_str_yield_table:
        length_file_name = len(name)
        i = 0
        i_start = 0
        i_end = 0
        while i < length_file_name:
            if name[i] == '=':
                i_start = i
            if name[i] == '.':
                i_end = i
            (i) = (i + 1)
        i = i_start + 1
        Z = ''
        while i < i_end:
            Z += name[i]
            (i) = (i + 1)
        Z_table_list += [float(Z)]
    sorted_Z_table_list = sorted(Z_table_list)
    # list 2
    file_names_setllar_lifetime_from_str_yield_table = os.listdir(
        yield_path + '/rearranged___/setllar_Metal_eject_mass_from_{}'.format(str_yield_table))
    Z_table_list_2 = []
    for name in file_names_setllar_lifetime_from_str_yield_table:
        length_file_name = len(name)
        i = 0
        i_start = 0
        i_end = 0
        while i < length_file_name:
            if name[i] == '=':
                i_start = i
            if name[i] == '.':
                i_end = i
            (i) = (i + 1)
        i = i_start + 1
        Z = ''
        while i < i_end:
            Z += name[i]
            (i) = (i + 1)
        if Z != '':
            Z_table_list_2 += [float(Z)]
    sorted_Z_table_list_2 = sorted(Z_table_list_2)
    if str_yield_table != "portinari98":
        # list 3
        file_names_setllar_lifetime_from_str_yield_table = os.listdir(
            yield_path + '/rearranged___/setllar_Metal_eject_mass_from_marigo01')
        Z_table_list_3 = []
        for name in file_names_setllar_lifetime_from_str_yield_table:
            length_file_name = len(name)
            i = 0
            i_start = 0
            i_end = 0
            while i < length_file_name:
                if name[i] == '=':
                    i_start = i
                if name[i] == '.':
                    i_end = i
                (i) = (i + 1)
            i = i_start + 1
            Z = ''
            while i < i_end:
                Z += name[i]
                (i) = (i + 1)
            Z_table_list_3 += [float(Z)]
        sorted_Z_table_list_3 = sorted(Z_table_list_3)
    else:
        sorted_Z_table_list_3 = []
    return (sorted_Z_table_list, sorted_Z_table_list_2, sorted_Z_table_list_3)


def function_select_metal(Z, Z_table_list):
    # the list for stellar lifetime is
    # [0.0004, 0.0008, 0.0012, 0.0016, 0.002, 0.0024, 0.0028, 0.0032, 0.0036, 0.004, 0.008, 0.012]
    # the list for stellar metallicity is
    # [0.0004, 0.004, 0.008, 0.0127] or [0, 0.004, 0.02] for Kobayashi2006 massive star table
    if Z <= Z_table_list[0]:
        Z_select__ = Z_table_list[0]
        return ('out', Z_select__, Z_select__, Z_select__)
        # The 'out' flag means the current gas metallicity is outside the range of provided stellar yield table.
    elif Z >= Z_table_list[-1]:
        Z_select__ = Z_table_list[-1]
        return ('out', Z_select__, Z_select__, Z_select__)
    else:
        i = 1
        while i < len(Z_table_list):
            if Z < Z_table_list[i]:
                Z_select__low = Z_table_list[i - 1]
                Z_select__high = Z_table_list[i]
                return ('in', Z_select__low, Z, Z_select__high)
            (i) = (i + 1)


def fucntion_mass_boundary_SNIa_Greggio83(time, mass):
    # Accroding to Greggio 1983 eq. 5:
    logM = math.log(mass, 10)
    left = math.log(time, 10)
    right = 10 - 4.319*logM + 1.543*(logM)**2
    if left > right:
        return fucntion_mass_boundary_SNIa_Greggio83(time, mass*0.99)
    else:
        return mass



def fucntion_mass_boundary(time, mass_grid_for_lifetime, lifetime):
    # The adopted spline fit of portinari98 lifetime is not monotonic at the massive end. 
    # But this function ensures that lifetime is monotonically smaller for more massive stars.
    mass = mass_grid_for_lifetime
    length_list_lifetime = len(lifetime)
    x = round(length_list_lifetime / 2)
    loop_number_fucntion_mass_boundary = math.ceil(math.log(length_list_lifetime, 2))
    mass_boundary = 10000
    if lifetime[x] == time:
        mass_boundary = mass[x]
    else:
        i = 0
        low = 0
        high = length_list_lifetime
        while i < loop_number_fucntion_mass_boundary:
            if lifetime[x] > time:
                low = x
                x = x + round((high - x) / 2)
            else:
                high = x
                x = x - round((x - low) / 2)
            (i) = (i + 1)
        if x == length_list_lifetime - 1:
            mass_boundary = mass[x]
        else:
            if lifetime[x - 1] > time > lifetime[x]:
                x = x - 1
            mass_boundary = round(mass[x] + (mass[x + 1] - mass[x]) * (lifetime[x] - time) / (
                lifetime[x] - lifetime[x + 1]), 5)
    return mass_boundary


# def function_get_observed_mass(lower_limit, upper_limit, M_tot_for_one_epoch, SFR, integrated_igimf):
#     integrator = quad(function_get_xi_from_IGIMF, lower_limit, upper_limit, SFR, limit=40)[0]
#     observed_mass = M_tot_for_one_epoch * integrator / integrated_igimf
#     return observed_mass


# def function_xi_Kroupa_IMF(mass):
#     # integrate this function's output xi result in the number of stars in mass limits.
#     xi = Kroupa_IMF.custom_imf(mass, 0)
#     return xi


# def function_mass_Kroupa_IMF(mass):
#     # integrate this function's output m result in the total stellar mass for stars in mass limits.
#     m = mass * Kroupa_IMF.custom_imf(mass, 0)
#     return m


def text_output(imf, STF, SFR, SFEN, original_gas_mass, log_Z_0):
    print('Generating txt output files...')
    global time_axis
    # print("time:", time_axis)

    global all_sf_imf
    number_of_sf_epoch = len(all_sf_imf)

    # data = exec(open("simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/chemical_and_SN_evolution.txt".format(IMF, SFE[0], SFR[0], SFEN[0])).read())
    #
    # print(data)

    #################
    ### F05 study ###
    #################
    #
    # mass_range_1 = [0.3, 0.4]
    # mass_boundary_low = all_sf_imf[0][1]
    # mass_boundary_high = all_sf_imf[-1][1]
    # print("mass_range_2_boundary_low", all_sf_imf[0][1])
    # print("mass_range_2_boundary_high", all_sf_imf[-1][1])
    # mass_range_2 = [mass_boundary_low, mass_boundary_high]

    # mass_range_1 = [0.3, 0.4]
    # mass_range_2 = [0.08, 1]
    #
    # integrate_Kroupa_stellar_mass_range_1 = quad(function_mass_Kroupa_IMF, mass_range_1[0], mass_range_1[1], limit=40)[0]
    # integrate_Kroupa_stellar_mass_range_2 = quad(function_mass_Kroupa_IMF, mass_range_2[0], mass_range_2[1], limit=40)[0]
    # integrate_Kroupa_stellar_number_mass_range_1 = \
    #     quad(function_xi_Kroupa_IMF, mass_range_1[0], mass_range_1[1], limit=40)[0]
    # integrate_Kroupa_stellar_number_mass_range_2 = \
    #     quad(function_xi_Kroupa_IMF, mass_range_2[0], mass_range_2[1], limit=40)[0]
    #
    # F_mass_Kroupa_IMF = integrate_Kroupa_stellar_mass_range_1 / integrate_Kroupa_stellar_mass_range_2
    # F_number_Kroupa_IMF = integrate_Kroupa_stellar_number_mass_range_1 / integrate_Kroupa_stellar_number_mass_range_2
    #
    # integrate_IGIMF_stellar_mass_range_1 = 0
    # integrate_IGIMF_stellar_mass_range_2 = 0
    # integrate_IGIMF_stellar_number_mass_range_1 = 0
    # integrate_IGIMF_stellar_number_mass_range_2 = 0
    # i = 0
    # while i < number_of_sf_epoch:
    #     def function_xi_IGIMF(mass):
    #         xi = all_sf_imf[i][0].custom_imf(mass, 0)
    #         return xi
    #
    #     def function_mass_IGIMF(mass):
    #         m = mass * all_sf_imf[i][0].custom_imf(mass, 0)
    #         return m
    #
    #     integrate_IGIMF_stellar_mass_range_1 += quad(function_mass_IGIMF, mass_range_1[0], mass_range_1[1], limit=40)[0]
    #     integrate_IGIMF_stellar_mass_range_2 += quad(function_mass_IGIMF, mass_range_2[0], mass_range_2[1], limit=40)[0]
    #     integrate_IGIMF_stellar_number_mass_range_1 += \
    #         quad(function_xi_IGIMF, mass_range_1[0], mass_range_1[1], limit=40)[0]
    #     integrate_IGIMF_stellar_number_mass_range_2 += \
    #         quad(function_xi_IGIMF, mass_range_2[0], mass_range_2[1], limit=40)[0]
    #     (i) = (i + 1)
    #
    # F_mass_IGIMF = integrate_IGIMF_stellar_mass_range_1 / integrate_IGIMF_stellar_mass_range_2
    # F_number_IGIMF = integrate_IGIMF_stellar_number_mass_range_1 / integrate_IGIMF_stellar_number_mass_range_2
    #
    # print("F_mass_Kroupa_IMF", F_mass_Kroupa_IMF)
    # print("F_mass_IGIMF", F_mass_IGIMF)
    # print("F_number_Kroupa_IMF", F_number_Kroupa_IMF)
    # print("F_number_IGIMF", F_number_IGIMF)

    # print("Number of star formation event epoch (10^7 yr): ", number_of_sf_epoch)
    # print("modeled star formation duration:", number_of_sf_epoch/100, "Gyr")
    global total_energy_release_list
    # print("total number of SN: 10^", round(math.log(total_energy_release_list[-1], 10), 1))

    global BH_mass_list, NS_mass_list, WD_mass_list, remnant_mass_list, stellar_mass_list, ejected_gas_mass_list
    stellar_mass = round(math.log(stellar_mass_list[-1], 10), 4)
    # print("Mass of all alive stars at final time: 10 ^", stellar_mass)
    downsizing_relation__star_formation_duration = round(10 ** (2.38 - 0.24 * stellar_mass), 4)  # Recchi 2009
    # print("star formation duration (downsizing relation):", downsizing_relation__star_formation_duration, "Gyr")

    stellar_and_remnant_mass = round(math.log(stellar_mass_list[-1] + remnant_mass_list[-1], 10), 4)
    # print("Mass of stars and remnants at final time: 10 ^", stellar_and_remnant_mass)

    total_mas_in_box = original_gas_mass

    # # Dabringhausen 2008 eq.4
    # Dabringhausen_2008_a = 2.95
    # Dabringhausen_2008_b = 0.596
    # expansion_factor = 5 ################ the expansion_factor should be a funtion of galaxy mass and rise with the mass
    # log_binding_energy = round(
    #     math.log(4.3 * 6 / 5, 10) + 40 + (2 - Dabringhausen_2008_b) * math.log(total_mas_in_box, 10) - math.log(
    #         Dabringhausen_2008_a, 10) + 6 * Dabringhausen_2008_b + math.log(expansion_factor, 10), 1)
    # # print("the gravitational binding energy: 10^", log_binding_energy, "erg")

    global Fe_over_H_list, stellar_Fe_over_H_list, stellar_Fe_over_H_list_luminosity_weighted
    # print("Gas [Fe/H]:", round(Fe_over_H_list[-1], 3))
    # print("Stellar [Fe/H]:", round(stellar_Fe_over_H_list[-1], 3))

    global Mg_over_Fe_list, stellar_Mg_over_Fe_list, stellar_Mg_over_Fe_list_luminosity_weighted
    # print("Gas [Mg/Fe]:", round(Mg_over_Fe_list[-1], 3))
    # print("Stellar [Mg/Fe]:", round(stellar_Mg_over_Fe_list[-1], 3))

    global O_over_Fe_list, stellar_O_over_Fe_list, stellar_O_over_Fe_list_luminosity_weighted
    # print("Gas [O/Fe]:", round(O_over_Fe_list[-1], 3))
    # print("Stellar [O/Fe]:", round(stellar_O_over_Fe_list[-1], 3))

    global Mg_over_H_list, stellar_Mg_over_H_list, stellar_Mg_over_H_list_luminosity_weighted
    global C_over_H_list, stellar_C_over_H_list, stellar_C_over_H_list_luminosity_weighted
    global N_over_H_list, stellar_N_over_H_list, stellar_N_over_H_list_luminosity_weighted
    global Ca_over_H_list, stellar_Ca_over_H_list, stellar_Ca_over_H_list_luminosity_weighted
    global Si_over_H_list, stellar_Si_over_H_list, stellar_Si_over_H_list_luminosity_weighted
    global S_over_H_list, stellar_S_over_H_list, stellar_S_over_H_list_luminosity_weighted
    global Ne_over_H_list, stellar_Ne_over_H_list, stellar_Ne_over_H_list_luminosity_weighted
    # print("Gas [Mg/H]:", round(Mg_over_H_list[-1], 3))
    # print("Stellar [Mg/H]:", round(stellar_Mg_over_H_list[-1], 3))

    global O_over_H_list, stellar_O_over_H_list, stellar_O_over_H_list_luminosity_weighted
    # print("Gas [O/H]:", round(O_over_H_list[-1], 3))
    # print("Stellar [O/H]:", round(stellar_O_over_H_list[-1], 3))

    global gas_Z_over_X_list, stellar_Z_over_X_list, stellar_Z_over_X_list_luminosity_weighted, stellar_Z_over_H_list
    # print("Gas metallicity:", round(gas_Z_over_X_list[-1], 3))
    # print("Stellar metallicity:", round(stellar_Z_over_X_list[-1], 3))
    # print("Stellar [Z/H]:", round(stellar_Z_over_H_list[-1], 3))

    filename = "simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/chemical_and_SN_evolution.txt".format(imf, STF, SFR, SFEN, log_Z_0)
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    file = open(
        'simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/chemical_and_SN_evolution.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')

    print("simulation results saved in the file: "
          "simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/(plots/)...txt".format(imf, STF, SFR, SFEN,
                                                                                                   log_Z_0))

    file.write("# Number of star formation event epoch (10^7 yr):\n")
    file.write("%s\n" % number_of_sf_epoch)

    file.write("# Modeled star formation duration (Gyr):\n")
    file.write("{}\n".format(number_of_sf_epoch / 100))

    file.write("# Total number of SN (log_10):\n")
    file.write("%s\n" % round(math.log(total_energy_release_list[-1], 10), 1))

    file.write("# Mass of all alive stars at final time (log_10):\n")
    file.write("%s\n" % stellar_mass)

    file.write("# Star formation duration of this final stellar mass according to the downsizing relation, Gyr):\n")
    file.write("%s\n" % downsizing_relation__star_formation_duration)

    file.write("# Mass of stars and remnants at final time (log_10):\n")
    file.write("%s\n" % stellar_and_remnant_mass)

    file.write("# total mass in the closed-box model:\n")
    file.write("%s\n" % total_mas_in_box)

    length_of_time_axis = len(time_axis)

    file.write("# time step list:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % time_axis[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Number of SNIa + SNII (log_10):\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % total_energy_release_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [Fe/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % Fe_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted [Fe/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Fe_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [Mg/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % Mg_over_Fe_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted [Mg/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Mg_over_Fe_list[i])
        (i) = (i + 1)
    file.write("\n")

    Mass_weighted_stellar_Mg_over_Fe = stellar_Mg_over_Fe_list[-1]
    # print("Mass-weighted stellar [Mg/Fe] at final time:", Mass_weighted_stellar_Mg_over_Fe)

    file.write("# Gas [O/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % O_over_Fe_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted [O/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_O_over_Fe_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [Mg/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % Mg_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [C/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % C_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [N/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % N_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [Ca/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % Ca_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [Ne/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % Ne_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [Si/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % Si_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [S/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % S_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted [Mg/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Mg_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted [C/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_C_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted [N/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_N_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted [Ca/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Ca_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted [Ne/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Ne_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted [Si/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Si_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted [S/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_S_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [O/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % O_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted [O/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_O_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas metallicity, [Z/X]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % gas_Z_over_X_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar mass-weighted metallicity, [Z/X]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Z_over_X_list[i])
        (i) = (i + 1)
    file.write("\n")

    Mass_weighted_stellar_metallicity = stellar_Z_over_X_list[-1]
    # print("Mass-weighted stellar [Z/X] at final time:", Mass_weighted_stellar_metallicity)

    if SNIa_energy_release_list[-1] < 10 ** (-10):
        SNIa_energy_release_list[-1] = 10 ** (-10)
    file.write("# Total number of SNIa (log_10):\n")
    file.write("%s\n" % round(math.log(SNIa_energy_release_list[-1], 10), 1))

    if SNII_energy_release_list[-1] < 10 ** (-10):
        SNII_energy_release_list[-1] = 10 ** (-10)
    file.write("# Total number of SNII (log_10):\n")
    file.write("%s\n" % round(math.log(SNII_energy_release_list[-1], 10), 1))

    file.write("# Number of SNIa (log_10):\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % SNIa_energy_release_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Number of SNII (log_10):\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % SNII_energy_release_list[i])
        (i) = (i + 1)
    file.write("\n")

    global ejected_gas_Mg_over_Fe_list
    file.write("# [Mg/Fe] for total ejected gas till this time:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % ejected_gas_Mg_over_Fe_list[i])
        (i) = (i + 1)
    file.write("\n")

    global instant_ejected_gas_Mg_over_Fe_list
    file.write("# [Mg/Fe] for instant ejected gas at this time:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % instant_ejected_gas_Mg_over_Fe_list[i])
        (i) = (i + 1)
    file.write("\n")

    global ejected_metal_mass_list
    file.write("# total ejected metal mass till this time:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % ejected_metal_mass_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# initial [M/H]:\n")
    file.write("%s " % log_Z_0)
    file.write("\n")

    file.write("# Stellar [Z/H]. [Z/H] = [Fe/H] + A[Mg/Fe] where A=0.94:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Z_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted metallicity, [Z/X]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Z_over_X_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [Mg/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Mg_over_Fe_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [C/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_C_over_Fe_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [N/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_N_over_O_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [O/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_O_over_Fe_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [Ca/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Ca_over_Fe_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [Ne/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Ne_over_Fe_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [Si/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Si_over_Fe_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [S/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_S_over_Fe_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [C/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_C_over_H_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [N/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_N_over_H_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [O/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_O_over_H_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [Mg/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Mg_over_H_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [Fe/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Fe_over_H_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [Ca/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Ca_over_H_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [Ne/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Ne_over_H_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [Si/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Si_over_H_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar luminosity-weighted [S/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_S_over_H_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# X_list:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % X_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# stellar_X_list:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_X_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# stellar_X_list_luminosity_weighted:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_X_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Y_list:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % Y_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# stellar_Y_list:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Y_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# stellar_Y_list_luminosity_weighted:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Y_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Z_list:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % Z_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# stellar_Z_list:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Z_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# stellar_Z_list_luminosity_weighted:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Z_list_luminosity_weighted[i])
        (i) = (i + 1)
    file.write("\n")

    file.close()

    return


def plot_output(plot_show, plot_save, imf, igimf, SFR, SFEN, log_Z_0, STF):  # SFR is the maximum SFR.
    if plot_show is True:
        print('\nGenerating plot outputs...\n')
    # plot SFH
    global all_sfr
    SFR_list = []
    age_list = []
    age_list.append(0)
    SFR_list.append(-22)
    age_list.append(0.01)
    SFR_list.append(-22)
    for i in range(len(all_sfr)):
        age_list.append(all_sfr[i][1])
        SFR_list.append(math.log(all_sfr[i][0], 10))
        age_list.append(all_sfr[i][1] + 0.01)
        SFR_list.append(math.log(all_sfr[i][0], 10))
    age_list.append(all_sfr[i][1] + 0.01)
    SFR_list.append(-22)
    age_list.append(10)
    SFR_list.append(-22)

    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(0, figsize=(3, 2.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(age_list, SFR_list, color="k", lw=0.8)
        # with open('SFH_Lacchin.txt') as f:
        #     lines = f.readlines()
        #     time_Lacchin = [float(line.split()[0]) for line in lines]
        #     SFH_Lacchin = [float(line.split()[1]) for line in lines]
        # with open('SFH_Lacchin_igimf.txt') as f:
        #     lines = f.readlines()
        #     time_Lacchin_igimf = [float(line.split()[0]) for line in lines]
        #     SFH_Lacchin_igimf = [float(line.split()[1]) for line in lines]
        # # plt.plot(time_Lacchin, SFH_Lacchin, color="tab:red", label='Lacchin2019 Salpeter', ls='dashed', lw=0.7)
        # plt.plot(time_Lacchin_igimf, SFH_Lacchin_igimf, color="tab:blue", label='Lacchin2019 IGIMF', ls='dashed', lw=0.7)
        plt.xlabel('time [Gyr]')
        plt.ylabel(r'log$_{10}(SFR [M_\odot$/yr])')
        # plt.xlim(6.4, 1.01 * log_time_axis[-1])
        # plt.xlim(-0.1, 1.1)
        # plt.ylim(-7, -2.5)
        # plt.title('Star formation history', fontsize=10)
        plt.tight_layout()
        # plt.legend(prop={'size': 7})
        if plot_save is True:
            plt.savefig('galaxy_evolution_fig_SFH.pdf', dpi=250)

    filename = "simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/SFH.txt".format(imf, STF, SFR, SFEN, log_Z_0)
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    length_of_SFH_list = len(SFR_list)
    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/SFH.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# age_list\n")
    i = 0
    while i < length_of_SFH_list:
        file.write("{} ".format(age_list[i]))
        (i) = (i + 1)
    file.write("\n# SFR_list\n")
    i = 0
    while i < length_of_SFH_list:
        file.write("{} ".format(SFR_list[i]))
        (i) = (i + 1)
    file.write("\n")
    file.close()

    # # plot IMF
    global all_sf_imf
    number_of_sf_epoch = len(all_sf_imf)
    mass_list = []
    xi_last_time = []
    xi_Kroupa = []
    xi_Salpeter = []
    xi_Kroupa_not_log = []
    xi_observe = []
    xi_each_epoch = []
    xi_each_time_log = []
    xi_each_time = []
    i = 0
    while i < number_of_sf_epoch:
        xi_each_epoch.append([])
        xi_each_time_log.append([])
        xi_each_time.append([])
        mass = 200
        while mass > 0.05:
            xi_each_epoch__ = all_sf_imf[i][0].custom_imf(mass, 0)
            if xi_each_epoch__ == 0:
                xi_each_epoch[i] += [-10]
            else:
                xi_each_epoch[i] += [math.log(xi_each_epoch__, 10)]
            j = 0
            xi_each_time__ = 0
            while j < i + 1:
                xi_each_time__ += all_sf_imf[j][0].custom_imf(mass, 0)
                (j) = (j + 1)
            if xi_each_time__ == 0:
                xi_each_time_log[i] += [-10]
                xi_each_time[i] += [0]
            else:
                xi_each_time_log[i] += [math.log(xi_each_time__, 10)]
                xi_each_time[i] += [xi_each_time__]
            (mass) = (mass * 0.99)
        (i) = (i + 1)

    j = 0
    xi_1_last_time = 0
    while j < number_of_sf_epoch:
        xi_1_last_time += all_sf_imf[j][0].custom_imf(1, 0)
        (j) = (j + 1)
    from IMFs import Salpeter_IMF
    normal_Kroupa = xi_1_last_time / Kroupa_IMF.custom_imf(1, 0)
    normal_Salpeter = xi_1_last_time / Salpeter_IMF.custom_imf(1, 0)

    mass = 200
    while mass > 0.05:
        mass_list += [mass]
        xi_last_time += [all_sf_imf[-1][0].custom_imf(mass, 0)]
        # xi_last_time += [igimf.custom_imf(mass, 0)]
        xi_observe__ = 0
        for i in range(number_of_sf_epoch):
            xi_observe__ += all_sf_imf[i][0].custom_imf(mass, 0)
            # if mass < all_sf_imf[i][1]:
            #     xi_observe__ += all_sf_imf[i][0].custom_imf(mass, 0)
        xi_observe += [xi_observe__]
        xi_Kroupa__ = Kroupa_IMF.custom_imf(mass, 0) * normal_Kroupa
        if xi_Kroupa__ == 0:
            xi_Kroupa += [-10]
        else:
            xi_Kroupa += [math.log(xi_Kroupa__, 10)]
        xi_Kroupa_not_log += [xi_Kroupa__]
        xi_Salpeter__ = Salpeter_IMF.custom_imf(mass, 0) * normal_Salpeter
        if xi_Salpeter__ == 0:
            xi_Salpeter += [-10]
        else:
            xi_Salpeter += [math.log(xi_Salpeter__, 10)]
        (mass) = (mass * 0.99)

    i = 0
    while i < number_of_sf_epoch:
        time = round(all_sf_imf[i][2] / 10 ** 6)
        file = open('Generated_IGIMFs/imf_at_time_{}_Myr.txt'.format(time), 'w')
        file.write("# This file gives the total number of stars in a unit mass interval for a given stellar mass "
                   "for the entire galaxy, i.e., galaxy-wide Initial Mass Function (gwIMF), at {} Myr.\n".format(time))
        file.write("# Below showing the given stellar mass on the left and the corresponding xi on the right, "
                   "where xi = d number / d mass.\n")
        mass_list_length = len(mass_list)
        j = mass_list_length - 1
        while j > 0:
            file.write("{} {}\n".format(mass_list[j], xi_each_time[i][j]))
            (j) = (j - 1)
        file.close()
        (i) = (i + 1)

    i = 0
    while i < number_of_sf_epoch:
        time = round(all_sf_imf[i][2] / 10 ** 6)
        length_of_xi = len(mass_list)
        file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/imf_at_time_{}_Myr.txt'.format(imf, STF, SFR, SFEN, log_Z_0, time), 'w')
        file.write("# mass_list\n")
        j = 0
        while j < length_of_xi:
            file.write("{} ".format(mass_list[j]))
            (j) = (j + 1)
        file.write("\n")
        file.write("# xi_each_time\n")
        j = 0
        while j < length_of_xi:
            file.write("{} ".format(xi_each_time[i][j]))
            (j) = (j + 1)
        file.write("\n")
        file.write("# xi_Kroupa\n")
        j = 0
        while j < length_of_xi:
            file.write("{} ".format(xi_Kroupa_not_log[j]))
            (j) = (j + 1)
        file.write("\n")
        file.close()
        (i) = (i + 1)

    for i in range(len(mass_list)):
        mass_list[i] = math.log(mass_list[i], 10)
        if xi_last_time[i] == 0:
            xi_last_time[i] = -10
        else:
            xi_last_time[i] = math.log(xi_last_time[i], 10)
        if xi_observe[i] == 0:
            xi_observe[i] = -10
        else:
            xi_observe[i] = math.log(xi_observe[i], 10)
            # if xi_Kroupa[i] == 0:
            #     xi_Kroupa[i] = -10
            # else:
            #     xi_Kroupa[i] = math.log(xi_Kroupa[i], 10)

    if plot_show is True or plot_save is True:
        # plt.rc('font', family='serif')
        # plt.rc('xtick', labelsize='x-small')
        # plt.rc('ytick', labelsize='x-small')
        # fig = plt.figure(1, figsize=(3, 2.5))
        # fig.add_subplot(1, 1, 1)
        # # i = 0
        # # while i < number_of_sf_epoch:
        # #     time = round(all_sf_imf[i][2] / 10**6)
        # #     if i < 3:
        # #         plt.plot(mass_list, xi_each_time_log[i], label='TIgwIMF at {} Myr'.format(time))
        # #     else:
        # #         plt.plot(mass_list, xi_each_time_log[i])
        # #     (i) = (i + 1)
        # # plt.plot(mass_list, xi_observe, label='final TIgwIMF')
        # plt.plot(mass_list, xi_observe, color='k', label='IGIMF', lw=0.9)
        # plt.plot(mass_list, xi_Salpeter, linestyle='dashed', color='r', label='Salpeter IMF')
        # plt.plot(mass_list, xi_Kroupa, linestyle='dotted', color='b', label='Kroupa IMF', lw=1.5)
        # plt.xlabel(r'log$_{10}(M_\star$ [$M_\odot$])')
        # plt.ylabel(r'log$_{10}(\xi_\star)$')
        # plt.ylim(-4, 8)
        # # plt.title('Time Integrated galaxy-wide IMF', fontsize=10)
        # plt.legend(prop={'size': 7})
        # plt.tight_layout()
        #
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(2, figsize=(3, 2.5))
        fig.add_subplot(1, 1, 1)
        i = 0
        while i < number_of_sf_epoch:
            plt.plot(mass_list, xi_each_epoch[i], lw=0.3, color='0.5')
            (i) = (i + 1)
        plt.plot(mass_list, xi_Salpeter, linestyle='dashed', color='tab:blue', label='Salpeter IMF', lw=0.7)
        plt.plot(mass_list, xi_observe, label='TIgwIMF', color='k')
        plt.plot(mass_list, xi_Kroupa, linestyle='dotted', color='r', label='Canonical IMF', lw=1.5)
        plt.xlabel(r'log$_{10}(M_\star)$ [$M_{\odot}$]')
        plt.ylabel(r'log$_{10}(\xi_\star)$')
        # plt.title('galaxy-wide IMF of each star formation epoch', fontsize=10)
        # plt.ylim(-4, 7)
        plt.legend(prop={'size': 7})
        plt.tight_layout()
        if plot_save is True:
            plt.savefig('galaxy_evolution_fig_TIgwIMF.pdf', dpi=250)


    ################# evolution plots #################
    global time_axis
    length_of_time_axis = len(time_axis)
    log_time_axis = []
    for i in range(length_of_time_axis):
        if time_axis[i] != 0:
            log_time_axis += [math.log((time_axis[i]), 10)]
        else:
            log_time_axis += [0]

    global X_list, stellar_X_list, stellar_X_list_luminosity_weighted
    global Y_list, stellar_Y_list, stellar_Y_list_luminosity_weighted
    global Z_list, stellar_Z_list, stellar_Z_list_luminosity_weighted
    global primary_He_mass_fraction
    DY_over_Z_list = [(yyy-primary_He_mass_fraction+1e-8)/zzz for zzz, yyy in zip(Z_list, Y_list)]
    stellar_DY_over_Z_list = [(yyy-primary_He_mass_fraction+1e-8)/zzz for zzz, yyy in zip(stellar_Z_list, stellar_Y_list)]
    stellar_DY_over_Z_list_luminosity_weighted = [(yyy-primary_He_mass_fraction+1e-8)/zzz for zzz, yyy in zip(stellar_Z_list_luminosity_weighted, stellar_Y_list_luminosity_weighted)]

    # dY_over_dZ = []
    # x_axis_for_dY_over_dZ = []
    # length_Z_list = len(Z_list)
    # i = 2
    # while i < length_Z_list:
    #     __dY_over_dZ__ = (Z_list[i] - Z_list[i-2])/(Y_list[i] - Y_list[i-2])
    #     __x_axis_for_dY_over_dZ__ = Z_list[i-1]
    #     dY_over_dZ.append(__dY_over_dZ__)
    #     x_axis_for_dY_over_dZ.append(__x_axis_for_dY_over_dZ__)
    #     (i) = (i+1)

    # print(stellar_Y_list)
    # print(stellar_Z_list)
    # print(stellar_DY_over_Z_list)

    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(61, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.stackplot(log_time_axis, Y_list, X_list, Z_list, labels=["Y", "X", "Z"])
    #     if plot_save is not True:
    #         plt.title('gas-phase H, He, and metal mass fraction', fontsize=10)
    #     plt.xlim(7, log_time_axis[-1])
    #     plt.ylim(0, 1)
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('stacked mass fraction')
    #     plt.legend(loc='lower left', prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('XYZ_gas_phase.pdf', dpi=250)
    #     fig = plt.figure(62, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.stackplot(log_time_axis, stellar_Y_list, stellar_X_list, stellar_Z_list, labels=["Y", "X", "Z"])
    #     if plot_save is not True:
    #         plt.title('stellar H, He, and metal mass fraction', fontsize=10)
    #     plt.xlim(7, log_time_axis[-1])
    #     plt.ylim(0, 1)
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('stacked mass fraction')
    #     plt.legend(loc='lower left', prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('XYZ_star_MW.pdf', dpi=250)
    #     fig = plt.figure(63, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.stackplot(log_time_axis, stellar_Y_list_luminosity_weighted, stellar_X_list_luminosity_weighted, stellar_Z_list_luminosity_weighted, labels=["Y", "X", "Z"])
    #     if plot_save is not True:
    #         plt.title('stellar luminosity-weighted H, He, and metal mass fraction', fontsize=10)
    #     plt.xlim(7, log_time_axis[-1])
    #     plt.ylim(0, 1)
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('stacked mass fraction')
    #     plt.legend(loc='lower left', prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('XYZ_star_LW.pdf', dpi=250)

    # global mm, zz
    # fig = plt.figure(0, figsize=(3, 2.5))
    # plt.plot(mm, zz)

    global Fe_over_H_list, stellar_Fe_over_H_list, stellar_Fe_over_H_list_luminosity_weighted

    print('plot stellar_Fe_over_H final', stellar_Fe_over_H_list[-1])

    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(3, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, Fe_over_H_list, label='gas')
    #     plt.plot(log_time_axis, stellar_Fe_over_H_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_Fe_over_H_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[Fe/H]')
    #     plt.title('Element abundance evolution', fontsize=10)
    #     # if imf == 'igimf':
    #     #     plt.title('IGIMF')
    #     # elif imf == 'Kroupa':
    #     #     plt.title('Kroupa IMF')
    #     # plt.legend(scatterpoints=1, numpoints=1, loc=0, prop={'size': 7}, ncol=2)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-5, 1)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_FeH_{}.pdf'.format(imf), dpi=250)

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/Fe_over_H_time.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# log_time_axis\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(log_time_axis[i]))
        (i) = (i + 1)
    file.write("\n# gas_Fe_over_H_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(Fe_over_H_list[i]))
        (i) = (i + 1)
    file.write("\n# stellar_Fe_over_H_list\n")
    i = 0
    while i < length_of_time_axis:
        if stellar_Fe_over_H_list[i] is None:
            file.write("{} ".format(-6))
        else:
            file.write("{} ".format(stellar_Fe_over_H_list[i]))
        (i) = (i + 1)
    file.write("\n# stellar_Fe_over_H_list_luminosity_weighted\n")
    i = 0
    while i < length_of_time_axis:
        if stellar_Fe_over_H_list_luminosity_weighted[i] is None:
            file.write("{} ".format(-6))
        else:
            file.write("{} ".format(stellar_Fe_over_H_list_luminosity_weighted[i]))
        (i) = (i + 1)
    file.write("\n")
    file.close()

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/Fe_over_H_mass.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# gas_Fe_over_H\n")
    file.write("{} ".format(Fe_over_H_list[-1]))
    file.write("\n# stellar_Fe_over_H\n")
    file.write("{} ".format(stellar_Fe_over_H_list[-1]))
    file.write("\n# stellar_Fe_over_H_list_luminosity_weighted\n")
    file.write("{} ".format(stellar_Fe_over_H_list_luminosity_weighted[-1]))
    file.write("\n")
    file.close()
    #
    global O_over_H_list, stellar_O_over_H_list, stellar_O_over_H_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(4, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, O_over_H_list, label='gas')
    #     plt.plot(log_time_axis, stellar_O_over_H_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_O_over_H_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[O/H]')
    #     plt.title('Element abundance evolution', fontsize=10)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-5, 1)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_OH_{}.pdf'.format(imf), dpi=250)
    #
    global Mg_over_H_list, stellar_Mg_over_H_list, stellar_Mg_over_H_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(5, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, Mg_over_H_list, label='gas')
    #     # print(stellar_Mg_over_H_list)
    #     plt.plot(log_time_axis, stellar_Mg_over_H_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_Mg_over_H_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[Mg/H]')
    #     plt.title('Element abundance evolution', fontsize=10)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-5, 1)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_MgH_{}.pdf'.format(imf), dpi=250)
    ###
    global C_over_H_list, stellar_C_over_H_list, stellar_C_over_H_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(31, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, C_over_H_list, label='gas')
    #     plt.plot(log_time_axis, stellar_C_over_H_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_C_over_H_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[C/H]')
    #     plt.title('Element abundance evolution', fontsize=10)
    #     # if imf == 'igimf':
    #     #     plt.title('IGIMF')
    #     # elif imf == 'Kroupa':
    #     #     plt.title('Kroupa IMF')
    #     # plt.legend(scatterpoints=1, numpoints=1, loc=0, prop={'size': 7}, ncol=2)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-5, 1)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_CH_{}.pdf'.format(imf), dpi=250)
    #
    global N_over_H_list, stellar_N_over_H_list, stellar_N_over_H_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(32, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, N_over_H_list, label='gas')
    #     plt.plot(log_time_axis, stellar_N_over_H_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_N_over_H_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[N/H]')
    #     plt.title('Element abundance evolution', fontsize=10)
    #     # if imf == 'igimf':
    #     #     plt.title('IGIMF')
    #     # elif imf == 'Kroupa':
    #     #     plt.title('Kroupa IMF')
    #     # plt.legend(scatterpoints=1, numpoints=1, loc=0, prop={'size': 7}, ncol=2)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-5, 1)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_NH_{}.pdf'.format(imf), dpi=250)
    #
    global Ca_over_H_list, stellar_Ca_over_H_list, stellar_Ca_over_H_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(33, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, Ca_over_H_list, label='gas')
    #     plt.plot(log_time_axis, stellar_Ca_over_H_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_Ca_over_H_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[Ca/H]')
    #     plt.title('Element abundance evolution', fontsize=10)
    #     # if imf == 'igimf':
    #     #     plt.title('IGIMF')
    #     # elif imf == 'Kroupa':
    #     #     plt.title('Kroupa IMF')
    #     # plt.legend(scatterpoints=1, numpoints=1, loc=0, prop={'size': 7}, ncol=2)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-5, 1)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_CaH_{}.pdf'.format(imf), dpi=250)
    #
    global Ne_over_H_list, stellar_Ne_over_H_list, stellar_Ne_over_H_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(331, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, Ne_over_H_list, label='gas')
    #     plt.plot(log_time_axis, stellar_Ne_over_H_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_Ne_over_H_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[Ne/H]')
    #     plt.title('Element abundance evolution', fontsize=10)
    #     # if imf == 'igimf':
    #     #     plt.title('IGIMF')
    #     # elif imf == 'Kroupa':
    #     #     plt.title('Kroupa IMF')
    #     # plt.legend(scatterpoints=1, numpoints=1, loc=0, prop={'size': 7}, ncol=2)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-5, 1)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_NeH_{}.pdf'.format(imf), dpi=250)
    #
    global Si_over_H_list, stellar_Si_over_H_list, stellar_Si_over_H_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(332, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, Si_over_H_list, label='gas')
    #     plt.plot(log_time_axis, stellar_Si_over_H_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_Si_over_H_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[Si/H]')
    #     plt.title('Element abundance evolution', fontsize=10)
    #     # if imf == 'igimf':
    #     #     plt.title('IGIMF')
    #     # elif imf == 'Kroupa':
    #     #     plt.title('Kroupa IMF')
    #     # plt.legend(scatterpoints=1, numpoints=1, loc=0, prop={'size': 7}, ncol=2)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-5, 1)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_SiH_{}.pdf'.format(imf), dpi=250)
    #
    global S_over_H_list, stellar_S_over_H_list, stellar_S_over_H_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(333, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, S_over_H_list, label='gas')
    #     plt.plot(log_time_axis, stellar_S_over_H_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_S_over_H_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[S/H]')
    #     plt.title('Element abundance evolution', fontsize=10)
    #     # if imf == 'igimf':
    #     #     plt.title('IGIMF')
    #     # elif imf == 'Kroupa':
    #     #     plt.title('Kroupa IMF')
    #     # plt.legend(scatterpoints=1, numpoints=1, loc=0, prop={'size': 7}, ncol=2)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-5, 1)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_SH_{}.pdf'.format(imf), dpi=250)
    ###
    # Reference: Serenelli & Basu 2010, Determining the Initial Helium Abundance of the Sun, DOI: 10.1088/0004-637X/719/1/865
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(64, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, Y_list, label='gas')
    #     plt.plot(log_time_axis, stellar_Y_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_Y_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [Y_solar, Y_solar], color='red',
    #              ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('Y')
    #     plt.title('Helium mass fraction evolution', fontsize=10)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_Y_{}.pdf'.format(imf), dpi=250)

    # if True: # plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(65, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, DY_over_Z_list, label='gas')
    #     plt.plot(log_time_axis, stellar_DY_over_Z_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_DY_over_Z_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [(Y_solar-primary_He_mass_fraction+1e-8)/Z_solar, (Y_solar-primary_He_mass_fraction+1e-8)/Z_solar], color='red',
    #              ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel(r'$\Delta$Y/Z')
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_DY_over_Z_{}.pdf'.format(imf), dpi=250)

    # global gas_Z_over_X_list, stellar_Z_over_X_list, stellar_Z_over_X_list_luminosity_weighted
    # if True: # plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(66, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(gas_Z_over_X_list, DY_over_Z_list, label='gas')
    #     plt.plot(stellar_Z_over_X_list, stellar_DY_over_Z_list, label='stellar MW')
    #     plt.plot(stellar_Z_over_X_list_luminosity_weighted, stellar_DY_over_Z_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([gas_Z_over_X_list[0], gas_Z_over_X_list[-1]], [(Y_solar-primary_He_mass_fraction+1e-8)/Z_solar, (Y_solar-primary_He_mass_fraction+1e-8)/Z_solar], color='red',
    #              ls='dashed', label='solar')
    #     plt.xlabel('[Z/X]')
    #     plt.ylabel('$\Delta$Y/Z')
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_DY_over_Z_over_Z_over_X_{}.pdf'.format(imf), dpi=250)
        
    # if True: # plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(67, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(Z_list, DY_over_Z_list, label='gas')
    #     plt.plot(stellar_Z_list, stellar_DY_over_Z_list, label='stellar MW')
    #     plt.plot(stellar_Z_list_luminosity_weighted, stellar_DY_over_Z_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([Z_list[0], Z_list[-1]], [(Y_solar-primary_He_mass_fraction+1e-8)/Z_solar, (Y_solar-primary_He_mass_fraction+1e-8)/Z_solar], color='red',
    #              ls='dashed', label='solar')
    #     plt.xlabel('Z')
    #     plt.ylabel('$\Delta$Y/Z')
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_DY_over_Z_over_Z_{}.pdf'.format(imf), dpi=250)

    # if True: # plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(68, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(Z_list, Y_list, label='gas')
    #     plt.plot(stellar_Z_list, stellar_Y_list, label='stellar MW')
    #     plt.plot(stellar_Z_list_luminosity_weighted, stellar_Y_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([Z_list[0], Z_list[-1]], [Y_solar, Y_solar], color='red',
    #              ls='dashed', label='solar')
    #     plt.xlabel('Z')
    #     plt.ylabel('Y')
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_Y_over_Z_{}.pdf'.format(imf), dpi=250)

    # if True: # plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(69, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(x_axis_for_dY_over_dZ, dY_over_dZ, label='gas')
    #     # plt.plot(stellar_Z_list, stellar_Y_list, label='stellar MW')
    #     # plt.plot(stellar_Z_list_luminosity_weighted, stellar_Y_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([x_axis_for_dY_over_dZ[0], x_axis_for_dY_over_dZ[-1]], [(Y_solar-primary_He_mass_fraction+1e-8)/Z_solar, (Y_solar-primary_He_mass_fraction+1e-8)/Z_solar], color='red',
    #              ls='dashed', label='solar')
    #     plt.xlabel('Z')
    #     plt.ylabel('dY/dZ')
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_dY_over_dZ_over_Z_{}.pdf'.format(imf), dpi=250)

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/Y_time.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# log_time_axis\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(log_time_axis[i]))
        (i) = (i + 1)
    file.write("\n# gas_Y_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(Y_list[i]))
        (i) = (i + 1)
    file.write("\n# stellar_Y_list\n")
    i = 0
    while i < length_of_time_axis:
        if stellar_Y_list[i] is None:
            file.write("0.001")
        else:
            file.write("{} ".format(stellar_Y_list[i]))
        (i) = (i + 1)
    file.write("\n# stellar_Y_list_luminosity_weighted\n")
    i = 0
    while i < length_of_time_axis:
        if stellar_Y_list_luminosity_weighted[i] is None:
            file.write("0.001")
        else:
            file.write("{} ".format(stellar_Y_list_luminosity_weighted[i]))
        (i) = (i + 1)
    file.write("\n")
    file.close()

    #
    global Mg_over_Fe_list, stellar_Mg_over_Fe_list, stellar_Mg_over_Fe_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(7, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, Mg_over_Fe_list, label='gas')
    #     plt.plot(log_time_axis, stellar_Mg_over_Fe_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_Mg_over_Fe_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[Mg/Fe]')
    #     plt.title('Element number ratio evolution', fontsize=10)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-1, 3.5)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_MgFe_{}.pdf'.format(imf), dpi=250)

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/Mg_over_Fe_time.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# log_time_axis\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(log_time_axis[i]))
        (i) = (i + 1)
    file.write("\n# gas_Mg_over_Fe_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(Mg_over_Fe_list[i]))
        (i) = (i + 1)
    file.write("\n# stellar_Mg_over_Fe_list\n")
    i = 0
    while i < length_of_time_axis:
        if stellar_Mg_over_Fe_list[i] is None:
            file.write("-6 ")
        else:
            file.write("{} ".format(stellar_Mg_over_Fe_list[i]))
        (i) = (i + 1)
    file.write("\n# stellar_Mg_over_Fe_list_luminosity_weighted\n")
    i = 0
    while i < length_of_time_axis:
        if stellar_Mg_over_Fe_list_luminosity_weighted[i] is None:
            file.write("-6 ")
        else:
            file.write("{} ".format(stellar_Mg_over_Fe_list_luminosity_weighted[i]))
        (i) = (i + 1)
    file.write("\n")
    file.close()

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/Mg_over_Fe_mass.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# gas_Mg_over_Fe\n")
    file.write("{} ".format(Mg_over_Fe_list[-1]))
    file.write("\n# stellar_Mg_over_Fe\n")
    file.write("{} ".format(stellar_Mg_over_Fe_list[-1]))
    file.write("\n# stellar_Mg_over_Fe_list_luminosity_weighted\n")
    file.write("{} ".format(stellar_Mg_over_Fe_list_luminosity_weighted[-1]))
    file.write("\n")
    file.close()

    #
    global O_over_Fe_list, stellar_O_over_Fe_list, stellar_O_over_Fe_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(8, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, O_over_Fe_list, label='gas')
    #     plt.plot(log_time_axis, stellar_O_over_Fe_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_O_over_Fe_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[O/Fe]')
    #     plt.title('Element number ratio evolution', fontsize=10)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-1, 3.5)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_MgFe_{}.pdf'.format(imf), dpi=250)
    #
    global Ca_over_Fe_list, stellar_Ca_over_Fe_list, stellar_Ca_over_Fe_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(81, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, Ca_over_Fe_list, label='gas')
    #     plt.plot(log_time_axis, stellar_Ca_over_Fe_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_Ca_over_Fe_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[Ca/Fe]')
    #     plt.title('Element number ratio evolution', fontsize=10)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-1, 3.5)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_CaFe_{}.pdf'.format(imf), dpi=250)
    #
    global Ne_over_Fe_list, stellar_Ne_over_Fe_list, stellar_Ne_over_Fe_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(811, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, Ne_over_Fe_list, label='gas')
    #     plt.plot(log_time_axis, stellar_Ne_over_Fe_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_Ne_over_Fe_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[Ne/Fe]')
    #     plt.title('Element number ratio evolution', fontsize=10)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-1, 3.5)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_NeFe_{}.pdf'.format(imf), dpi=250)
    #
    global Si_over_Fe_list, stellar_Si_over_Fe_list, stellar_Si_over_Fe_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(812, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, Si_over_Fe_list, label='gas')
    #     plt.plot(log_time_axis, stellar_Si_over_Fe_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_Si_over_Fe_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[Si/Fe]')
    #     plt.title('Element number ratio evolution', fontsize=10)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-1, 3.5)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_SiFe_{}.pdf'.format(imf), dpi=250)
    #
    global S_over_Fe_list, stellar_S_over_Fe_list, stellar_S_over_Fe_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(813, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, S_over_Fe_list, label='gas')
    #     plt.plot(log_time_axis, stellar_S_over_Fe_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_S_over_Fe_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[S/Fe]')
    #     plt.title('Element number ratio evolution', fontsize=10)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-1, 3.5)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_SFe_{}.pdf'.format(imf), dpi=250)
    #
    global C_over_Fe_list, stellar_C_over_Fe_list, stellar_C_over_Fe_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(82, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, C_over_Fe_list, label='gas')
    #     plt.plot(log_time_axis, stellar_C_over_Fe_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_C_over_Fe_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[C/Fe]')
    #     plt.title('Element number ratio evolution', fontsize=10)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-1, 3.5)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_CFe_{}.pdf'.format(imf), dpi=250)
    # #
    # global N_over_O_list, stellar_N_over_O_list, stellar_N_over_O_list_luminosity_weighted
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(83, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(log_time_axis, N_over_O_list, label='gas')
    #     plt.plot(log_time_axis, stellar_N_over_O_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_N_over_O_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel('[N/Fe]')
    #     plt.title('Element number ratio evolution', fontsize=10)
    #     # plt.xlim(6.4, 1.01 * log_time_axis[-1])
    #     # plt.ylim(-1, 3.5)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_NFe_{}.pdf'.format(imf), dpi=250)
    #
    del stellar_Fe_over_H_list[0]
    del stellar_Fe_over_H_list_luminosity_weighted[0]
    del stellar_Mg_over_Fe_list[0]
    del stellar_Mg_over_Fe_list_luminosity_weighted[0]
    del stellar_O_over_Fe_list[0]
    del stellar_O_over_Fe_list_luminosity_weighted[0]
    del stellar_Si_over_Fe_list[0]
    del stellar_Si_over_Fe_list_luminosity_weighted[0]
    del stellar_Ca_over_Fe_list[0]
    del stellar_Ca_over_Fe_list_luminosity_weighted[0]

    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(9, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     Fe_over_H_list[0] = -99
    #     stellar_Fe_over_H_list[0] = -99
    #     stellar_Fe_over_H_list_luminosity_weighted[0] = -99
    #     plt.plot(Fe_over_H_list, Mg_over_Fe_list, label='gas')
    #     # plt.scatter(Fe_over_H_list, Mg_over_Fe_list, alpha=0.5, s=10)
    #     plt.plot(stellar_Fe_over_H_list, stellar_Mg_over_Fe_list, label='stellar MW')
    #     plt.plot(stellar_Fe_over_H_list_luminosity_weighted, stellar_Mg_over_Fe_list_luminosity_weighted,
    #              label='stellar LW')
    #     plt.plot([-4, -0.6], [-0.8, -0.8], color='red', label='Lacchin2019', lw=0.5)
    #     plt.plot([-4, -0.6], [1, 1], color='red', lw=0.5)
    #     plt.plot([-4, -4], [-0.8, 1], color='red', lw=0.5)
    #     plt.plot([-0.6, -0.6], [-0.8, 1], color='red', lw=0.5)
    #     plt.plot([-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.6], [0.44, 0.445, 0.45, 0.44, 0.37, 0.17, -0.14, -0.44], color='red', ls='dashed')
    #     plt.scatter([-3.7, -3.2, -2.92, -2.7, -2.5, -2.3, -2.3, -2.2, -2.05, -1.9, -1.8],
    #                 [0.47, 0.4, 0.65, 0.15, 0.2, 0.32, 0.35, 0.35, 0.3, 0.2, 0.48], color='k', alpha=0.5)
    #     with open('Mg_Lacchin.txt') as f:
    #         lines = f.readlines()
    #         time_Mg_Lacchin_igimf = [float(line.split()[0]) for line in lines]
    #         Mg_Lacchin_igimf = [float(line.split()[1]) for line in lines]
    #     plt.plot(time_Mg_Lacchin_igimf, Mg_Lacchin_igimf, color="tab:blue", label='Lacchin2019 IGIMF', ls='dashed')
    #     plt.xlabel('[Fe/H]')
    #     plt.ylabel('[Mg/Fe]')
    #     plt.xlim(-4.5, 0)
    #     plt.ylim(-1, 1.5)
    #     plt.legend(loc='lower left', prop={'size': 6})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_MgFe-FeH_{}.pdf'.format(imf), dpi=250)
    #
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(91, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(O_over_H_list, N_over_O_list, label='gas')
    #     plt.plot(stellar_O_over_H_list, stellar_N_over_O_list, label='stellar MW')
    #     plt.plot(stellar_O_over_H_list_luminosity_weighted, stellar_N_over_O_list_luminosity_weighted,
    #              label='stellar LW')
    #     plt.plot([-5, 1], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.plot([0, 0], [-1, 3.5], color='red', ls='dashed')
    #     plt.xlabel('[O/H]')
    #     plt.ylabel('[N/O]')
    #     # plt.xlim(-5, 1)
    #     # plt.ylim(-1, 3.5)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_NO-OH_{}.pdf'.format(imf), dpi=250)
    # #
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(92, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(O_over_H_list, C_over_H_list, label='gas')
    #     plt.plot(stellar_O_over_H_list, stellar_C_over_H_list, label='stellar MW')
    #     plt.plot(stellar_O_over_H_list_luminosity_weighted, stellar_C_over_H_list_luminosity_weighted,
    #              label='stellar LW')
    #     plt.plot([-5, 1], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.plot([0, 0], [-1, 3.5], color='red', ls='dashed')
    #     plt.xlabel('[O/H]')
    #     plt.ylabel('[C/H]')
    #     # plt.xlim(-5, 1)
    #     # plt.ylim(-1, 3.5)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_CH-OH_{}.pdf'.format(imf), dpi=250)
    # #
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(93, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(Fe_over_H_list, Si_over_Fe_list, label='gas')
    #     # plt.scatter(Fe_over_H_list, Si_over_Fe_list, alpha=0.5, s=10)
    #     plt.plot(stellar_Fe_over_H_list, stellar_Si_over_Fe_list, label='stellar MW')
    #     plt.plot(stellar_Fe_over_H_list_luminosity_weighted, stellar_Si_over_Fe_list_luminosity_weighted,
    #              label='stellar LW')
    #     plt.plot([-4, -0.6], [-0.5, -0.5], color='red', label='Lacchin2019', lw=0.5)
    #     plt.plot([-4, -0.6], [1, 1], color='red', lw=0.5)
    #     plt.plot([-4, -4], [-0.5, 1], color='red', lw=0.5)
    #     plt.plot([-0.6, -0.6], [-0.5, 1], color='red', lw=0.5)
    #     plt.plot([-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.6], [0.68, 0.64, 0.6, 0.56, 0.5, 0.33, 0.1, -0.07], color='red', ls='dashed')
    #     plt.scatter([-3.68, -2.05, -1.93],
    #                 [0.77, 0.1, 0.13], color='k', alpha=0.5)
    #     with open('Si_Lacchin.txt') as f:
    #         lines = f.readlines()
    #         time_Si_Lacchin_igimf = [float(line.split()[0]) for line in lines]
    #         Si_Lacchin_igimf = [float(line.split()[1]) for line in lines]
    #     plt.plot(time_Si_Lacchin_igimf, Si_Lacchin_igimf, color="tab:blue", label='Lacchin2019 IGIMF', ls='dashed')
    #     plt.xlabel('[Fe/H]')
    #     plt.ylabel('[Si/Fe]')
    #     plt.xlim(-4.5, 0)
    #     plt.ylim(-1, 1.5)
    #     plt.legend(loc='lower left', prop={'size': 6})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_SiFe-FeH_{}.pdf'.format(imf), dpi=250)
    # #
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(94, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(Fe_over_H_list, Ca_over_Fe_list, label='gas')
    #     # plt.scatter(Fe_over_H_list, Ca_over_Fe_list, alpha=0.5, s=10)
    #     plt.plot(stellar_Fe_over_H_list, stellar_Ca_over_Fe_list, label='stellar MW')
    #     plt.plot(stellar_Fe_over_H_list_luminosity_weighted, stellar_Ca_over_Fe_list_luminosity_weighted,
    #              label='stellar LW')
    #     plt.plot([-4, -0.6], [-0.5, -0.5], color='red', label='Lacchin2019', lw=0.5)
    #     plt.plot([-4, -0.6], [0.7, 0.7], color='red', lw=0.5)
    #     plt.plot([-4, -4], [-0.5, 0.7], color='red', lw=0.5)
    #     plt.plot([-0.6, -0.6], [-0.5, 0.7], color='red', lw=0.5)
    #     plt.plot([-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.6], [0.3, 0.27, 0.22, 0.19, 0.14, 0.01, -0.15, -0.26], color='red', ls='dashed')
    #     plt.xlabel('[Fe/H]')
    #     plt.ylabel('[Ca/Fe]')
    #     plt.scatter([-3.68, -3.2, -2.91, -2.7, -2.5, -2.28, -2.28, -2.2, -2.06, -1.94, -1.8],
    #                 [0.52, 0.46, 0.32, 0.18, 0.24, 0.28, 0.32, -0.01, 0.15, 0.1, 0.22], color='k', alpha=0.5)
    #     with open('Ca_Lacchin.txt') as f:
    #         lines = f.readlines()
    #         time_Ca_Lacchin_igimf = [float(line.split()[0]) for line in lines]
    #         Ca_Lacchin_igimf = [float(line.split()[1]) for line in lines]
    #     plt.plot(time_Ca_Lacchin_igimf, Ca_Lacchin_igimf, color="tab:blue", label='Lacchin2019 IGIMF', ls='dashed')
    #     plt.xlim(-4.5, 0)
    #     plt.ylim(-1, 1.5)
    #     plt.legend(loc='lower left', prop={'size': 6})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_CaFe-FeH_{}.pdf'.format(imf), dpi=250)
    #
    #
    #
    # if plot_show is True or plot_save is True:
    #     fig, axs = plt.subplots(3, 1, sharex=True, figsize=(3, 5))
    #     axs[0].plot(Fe_over_H_list, Mg_over_Fe_list, color='k')
    #     with open('Mg_Lacchin_best.txt') as f:
    #         lines = f.readlines()
    #         time_Mg_Lacchin_best = [float(line.split()[0]) for line in lines]
    #         Mg_Lacchin_best = [float(line.split()[1]) for line in lines]
    #     Lai2011_Fe_H = [-2.34 - 0.05, -2.37 - 0.05, -3.09 - 0.05, -2.39 - 0.05, -2.89 - 0.05, -2.2 - 0.05, -2.49 - 0.05,
    #                     -2.48 - 0.05, -2.65 - 0.05, -2.43 - 0.05, -2.48 - 0.05, -2.49 - 0.05, -2.57 - 0.05, -2.89 - 0.05,
    #                     -2.51 - 0.05, -3.29 - 0.05, -2.42 - 0.05, -3.79 - 0.05, -2.87 - 0.05, -2.9 - 0.05, -1.65 - 0.05,
    #                     -2.76 - 0.05, -2.31 - 0.05, -1.86 - 0.05]
    #     Lai2011_Mg_Fe = [0.12 + 0.05, 0.28 + 0.05, 0.43 + 0.05, 0.35 + 0.05, 0.05 + 0.05, 0.32 + 0.05, 0.16 + 0.05,
    #                      0.18 + 0.05, 0.37 + 0.05, 0.05 + 0.05, 0.04 + 0.05, 0.23 + 0.05, 0.34 + 0.05, 0.3 + 0.05, 0.13 + 0.05,
    #                      0.28 + 0.05, 0.13 + 0.05, 0.27 + 0.05, 0.12 + 0.05, 0.16 + 0.05, 0.46 + 0.05, 0.25 + 0.05,
    #                      0.06 + 0.05, 0.51 + 0.05]
    #     Lai2011_Fe_H_error = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
    #                           0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
    #     Lai2011_Mg_Fe_error = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
    #                            0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    #     Lacchin2019_Fe_H = [-3.6613, -3.2003, -2.9075, -2.7063, -2.5025, -2.2810, -2.2810, -2.1969, -2.0492, -1.9219, -1.8048]
    #     Lacchin2019_Mg_Fe = [0.47409, 0.39382, 0.63630, 0.13418, 0.20614, 0.35587, 0.32298, 0.35595, 0.30228, 0.20672, 0.47594]
    #     Lacchin2019_Fe_H_error = [0.1, 0.15, 0.15, 0.15, 0.15, 0.3, 0.15, 0.15, 0.15, 0.15, 0.15]
    #     Lacchin2019_Mg_Fe_error = [0.1, 0.05, 0.12, 0.135, 0.2, 0.23, 0.235, 0.13, 0.23, 0.14, 0.24]
    #     # axs[0].errorbar(Lai2011_Fe_H, Lai2011_Mg_Fe, xerr=Lai2011_Fe_H_error, yerr=Lai2011_Mg_Fe_error, markersize=0,
    #     #                 ecolor='k', ls='none', elinewidth=0.5, alpha=0.4)
    #     axs[0].errorbar(Lacchin2019_Fe_H, Lacchin2019_Mg_Fe, xerr=Lacchin2019_Fe_H_error, yerr=Lacchin2019_Mg_Fe_error, markersize=0,
    #                     ecolor='tab:blue', ls='none', elinewidth=0.5, alpha=0.4)
    #     if imf == 'igimf':
    #         axs[0].plot(time_Mg_Lacchin_igimf, Mg_Lacchin_igimf, color='tab:blue', ls='dashed', lw=0.7)
    #         axs[0].plot(time_Mg_Lacchin_best, Mg_Lacchin_best, color='tab:blue', ls='-.', lw=0.7)
    #     elif imf == 'Salpeter':
    #         axs[0].plot([-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.6], [0.44, 0.445, 0.45, 0.44, 0.37, 0.17, -0.14, -0.44],
    #                  color='red', ls='dashed', lw=0.7)
    #     axs[0].scatter(Lacchin2019_Fe_H, Lacchin2019_Mg_Fe, color='tab:blue', alpha=0.5, label='Lacchin2019', s=20)
    #     # axs[0].scatter(Lai2011_Fe_H, Lai2011_Mg_Fe, color='k', alpha=0.5, marker='*', s=20, label='Lai2011')
    #     axs[0].set_xlabel('[Fe/H]')
    #     axs[0].set_ylabel('[Mg/Fe]')
    #     axs[0].set_xlim(-4, -0.6)
    #     axs[0].set_ylim(-0.8, 1)
    #     # axs[0].legend(prop={'size': 6}, loc='lower left')
    #
    #     with open('Si_Lacchin_best.txt') as f:
    #         lines = f.readlines()
    #         time_Si_Lacchin_best = [float(line.split()[0]) for line in lines]
    #         Si_Lacchin_best = [float(line.split()[1]) for line in lines]
    #     Lacchin2019_Fe_H = [-3.6624, -2.0549, -1.9241]
    #     Lacchin2019_Si_Fe = [0.77012, 0.10533, 0.13018]
    #     Lacchin2019_Fe_H_error = [0.11, 0.15, 0.15]
    #     Lacchin2019_Si_Fe_error = [0.15, 0.26, 0.23]
    #     axs[1].plot(Fe_over_H_list, Si_over_Fe_list, color='k')
    #     axs[1].errorbar(Lacchin2019_Fe_H, Lacchin2019_Si_Fe, xerr=Lacchin2019_Fe_H_error, yerr=Lacchin2019_Si_Fe_error,
    #                     markersize=0, ecolor='tab:blue', ls='none', elinewidth=0.5, alpha=0.4)
    #     if imf == 'igimf':
    #         axs[1].plot(time_Si_Lacchin_best, Si_Lacchin_best, color='tab:blue', ls='-.', lw=0.7)
    #         axs[1].plot(time_Si_Lacchin_igimf, Si_Lacchin_igimf, color='tab:blue', ls='dashed', lw=0.7)
    #     elif imf == 'Salpeter':
    #         axs[1].plot([-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.6], [0.68, 0.64, 0.6, 0.56, 0.5, 0.33, 0.1, -0.07],
    #                     color='red', ls='dashed', lw=0.7)
    #     axs[1].scatter(Lacchin2019_Fe_H, Lacchin2019_Si_Fe, color='tab:blue', alpha=0.5, label='Lacchin2019', s=20)
    #     axs[1].set_xlabel('[Fe/H]')
    #     axs[1].set_ylabel('[Si/Fe]')
    #     axs[1].set_xlim(-4, -0.6)
    #     axs[1].set_ylim(-0.5, 1)
    #
    #     with open('Ca_Lacchin_best.txt') as f:
    #         lines = f.readlines()
    #         time_Ca_Lacchin_best = [float(line.split()[0]) for line in lines]
    #         Ca_Lacchin_best = [float(line.split()[1]) for line in lines]
    #     Lacchin2019_Fe_H = [-3.6584, -3.1993, -2.9072, -2.7045, -2.5011, -2.2812, -2.2813, -2.2820, -2.1992, -2.0522, -1.9244, -1.8019, ]
    #     Lacchin2019_Ca_Fe = [0.52167, 0.45923, 0.32191, 0.18145, 0.23876, 0.31201, 0.28968, 0.13178, -0.010205, 0.14921, 0.098101, 0.22084]
    #     Lacchin2019_Fe_H_error = [0.11, 0.15, 0.16, 0.16, 0.16, 0.11, 0.3, 0.16, 0.16, 0.16, 0.16, 0.16]
    #     Lacchin2019_Ca_Fe_error = [0.1, 0.05, 0.05, 0.05, 0.08, 0.05, 0.28, 0.05, 0.11, 0.07, 0.05, 0.06]
    #     axs[2].plot(Fe_over_H_list, Ca_over_Fe_list, color='k')
    #     axs[2].errorbar(Lacchin2019_Fe_H, Lacchin2019_Ca_Fe, xerr=Lacchin2019_Fe_H_error, yerr=Lacchin2019_Ca_Fe_error,
    #                     markersize=0, ecolor='tab:blue', ls='none', elinewidth=0.5, alpha=0.4)
    #     if imf == 'igimf':
    #         axs[2].plot(time_Ca_Lacchin_best, Ca_Lacchin_best, color='tab:blue', ls='-.', lw=0.7)
    #         axs[2].plot(time_Ca_Lacchin_igimf, Ca_Lacchin_igimf, color='tab:blue', ls='dashed', lw=0.7)
    #     elif imf == 'Salpeter':
    #         axs[2].plot([-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.6], [0.3, 0.27, 0.22, 0.19, 0.14, 0.01, -0.15, -0.26], color='red', ls='dashed', lw=0.7)
    #     axs[2].scatter(Lacchin2019_Fe_H, Lacchin2019_Ca_Fe, color='tab:blue', alpha=0.5, label='Lacchin2019', s=20)
    #     axs[2].set_xlabel('[Fe/H]')
    #     axs[2].set_ylabel('[Ca/Fe]')
    #     axs[2].set_xlim(-4, -0.6)
    #     axs[2].set_ylim(-0.5, 0.7)
    #
    #     plt.tight_layout()
    #     fig.subplots_adjust(hspace=0)  # Remove horizontal space between axes
    #     plt.savefig('MgSiCa.pdf', dpi=250)
    #
    #
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(10, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(Fe_over_H_list, O_over_Fe_list, label='gas')
    #     plt.plot(stellar_Fe_over_H_list, stellar_O_over_Fe_list, label='stellar MW')
    #     plt.plot(stellar_Fe_over_H_list_luminosity_weighted, stellar_O_over_Fe_list_luminosity_weighted,
    #              label='stellar LW')
    #     plt.plot([-5, 1], [0, 0], color='red', ls='dashed', label='solar')
    #     plt.plot([0, 0], [-1, 3.5], color='red', ls='dashed')
    #     plt.xlabel('[Fe/H]')
    #     plt.ylabel('[O/Fe]')
    #     # plt.xlim(-5, 1)
    #     # plt.ylim(-1, 3.5)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_OFe-FeH_{}.pdf'.format(imf), dpi=250)
    #
    global SNIa_number_per_century, SNII_number_per_century
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(11, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.loglog(time_axis, SNIa_number_per_century, label='SNIa', color="tab:orange",
    #                ls='dotted')  # Number per century
    #     plt.loglog(time_axis, SNII_number_per_century, label='SNII', color="tab:orange")  # Number per century
    #     # plt.loglog(time_axis, SN_number_per_century, ls="dotted", label='total')
    #     plt.xlabel(r'time [yr]')
    #     plt.ylabel(r'# of SN per century')
    #     plt.title('Supernova rate evolution', fontsize=10)
    #     plt.xlim(10 ** 7, 14 * 10 ** 9)
    #     # plt.ylim(1e-2, 1e6)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_SN_number.pdf'.format(imf), dpi=250)
    #
    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(111, figsize=(3, 2.5))
        fig.add_subplot(1, 1, 1)
        SNIa_number_per_yr = [x / 100 for x in SNIa_number_per_century]
        time_in_Gyr = [x / 10**9 for x in time_axis]
        plt.plot(time_in_Gyr, SNIa_number_per_yr, color="k", lw=0.8)
        # with open('SNIa_Lacchin.txt') as f:
        #     lines = f.readlines()
        #     time_Lacchin = [float(line.split()[0]) for line in lines]
        #     SNIa_Lacchin = [10 ** float(line.split()[1]) for line in lines]
        # with open('SNIa_Lacchin_igimf.txt') as f:
        #     lines = f.readlines()
        #     time_Lacchin_igimf = [float(line.split()[0]) for line in lines]
        #     SNIa_Lacchin_igimf = [10 ** float(line.split()[1]) for line in lines]
        # # plt.plot(time_Lacchin, SNIa_Lacchin, color="tab:red", label='Lacchin2019 Salpeter', ls='dashed', lw=0.7)
        # plt.plot(time_Lacchin_igimf, SNIa_Lacchin_igimf, color="tab:blue", label='Lacchin2019 IGIMF', ls='dashed', lw=0.7)
        plt.yscale('log')
        plt.xlabel(r'time [Gyr]')
        plt.ylabel(r'# of SNIa per yr')
        # plt.xlim(0, 14)
        # plt.ylim(4e-10, 5e-7)
        # plt.ylim(1e-9, 1e-6)
        # plt.legend(prop={'size': 7})
        plt.tight_layout()
        if plot_save is True:
            plt.savefig('galaxy_evolution_SNIa_number_loglinear.pdf'.format(imf), dpi=250)
    #
    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(112, figsize=(3, 2.5))
        fig.add_subplot(1, 1, 1)
        SNII_number_per_yr = [x / 100 for x in SNII_number_per_century]
        time_in_Gyr = [x / 10**9 for x in time_axis]
        plt.plot(time_in_Gyr, SNII_number_per_yr, color="k", lw=0.8)
        # with open('SNII_Lacchin.txt') as f:
        #     lines = f.readlines()
        #     time_Lacchin = [float(line.split()[0]) for line in lines]
        #     SNII_Lacchin = [10**float(line.split()[1]) for line in lines]
        # with open('SNII_Lacchin_igimf.txt') as f:
        #     lines = f.readlines()
        #     time_Lacchin_igimf = [float(line.split()[0]) for line in lines]
        #     SNII_Lacchin_igimf = [10**float(line.split()[1]) for line in lines]
        # # plt.plot(time_Lacchin, SNII_Lacchin, color="tab:red", label='Lacchin2019 Salpeter', ls='dashed', lw=0.7)
        # plt.plot(time_Lacchin_igimf, SNII_Lacchin_igimf, color="tab:blue", label='Lacchin2019 IGIMF', ls='dashed', lw=0.7)
        plt.yscale('log')
        plt.xlabel(r'time [Gyr]')
        plt.ylabel(r'# of SNII per yr')
        # plt.xlim(0, 1)
        # plt.ylim(1e-9, 1e-5)
        # plt.legend(prop={'size': 7})
        plt.tight_layout()
        if plot_save is True:
            plt.savefig('galaxy_evolution_SNII_number_loglinear.pdf'.format(imf), dpi=250)

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/SN_number_evolution.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# time_axis\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(time_axis[i]))
        (i) = (i + 1)
    file.write("\n# SNIa_number_per_century\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(SNIa_number_per_century[i]))
        (i) = (i + 1)
    file.write("\n# SNII_number_per_century\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(SNII_number_per_century[i]))
        (i) = (i + 1)
    file.write("\n# SN_number_per_century\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(SN_number_per_century[i]))
        (i) = (i + 1)
    file.write("\n")
    file.close()

    # calculate the gravitational binding energy:
    global expansion_factor_list, original_gas_mass
    # consider the following to calculate the energy:
    # galaxy_mass_without_gas_at_this_time,
    # original_gas_mass,
    # total_gas_mass_at_this_time,
    # ejected_gas_mass_at_this_time,
    # gas_mass = max(ejected_gas_mass_at_this_time, 1)

    # galaxy mass--radii relation adopted from Dabringhausen 2008 eq.4
    Dabringhausen_2008_a = 2.95
    Dabringhausen_2008_b = 0.596
    final_expansion_factor = expansion_factor_list[-1]
    binding_energy_list = []
    SN_energy_per_current_crossing_time_list = []
    SN_energy_per_final_crossing_time_list = []
    length_expansion_factor_list = len(expansion_factor_list)

    global remnant_mass_list, stellar_mass_list
    gravitational_constant = 6.674
    finial_galaxy_inner_mass = remnant_mass_list[-1] + stellar_mass_list[-1]
    final_galaxy_radii = 3 * (finial_galaxy_inner_mass / 10 ** 6) ** 0.6  # in parsec
    final_crossing_time = 1 / (gravitational_constant * 10 ** (-11) * finial_galaxy_inner_mass * 2 * 10 ** 30 / (
    final_galaxy_radii * 3 * 10 ** 16) ** 3) ** (0.5) / (60 * 60 * 24 * 365 * 10 ** 6)  # in Myr

    i = 0
    while i < length_expansion_factor_list:
        ### binding energy ###
        current_shrink_factor = final_expansion_factor / expansion_factor_list[i]
        log_binding_energy = round(
            math.log(3 / 5 * gravitational_constant * 1.989 ** 2 / 3.086, 10) + 40 + (2 - Dabringhausen_2008_b) *
            math.log(original_gas_mass, 10) - math.log(Dabringhausen_2008_a, 10) +
            6 * Dabringhausen_2008_b + math.log(current_shrink_factor, 10), 3)
        # 40 = 30 (solar mass) * 2 - 11 (Gravitational constant) - 16 (pc to meter) + 7 (J to erg)
        binding_energy = 10 ** log_binding_energy  # [erg]
        binding_energy_list.append(binding_energy)
        ### crossing time ###
        current_galaxy_inner_mass = remnant_mass_list[i] + stellar_mass_list[i]
        current_galaxy_radii = final_galaxy_radii / current_shrink_factor
        crossing_time = 1 / (gravitational_constant * 10 ** (-11) * current_galaxy_inner_mass * 2 * 10 ** 30 / (
        current_galaxy_radii * 3 * 10 ** 16) ** 3) ** (0.5) / (60 * 60 * 24 * 365 * 10 ** 6)  # in Myr
        SN_energy_per_current_crossing_time = SN_number_per_century[i] * crossing_time * 10 ** 4 * 10 ** 51
        SN_energy_per_final_crossing_time = SN_number_per_century[i] * final_crossing_time * 10 ** 4 * 10 ** 51
        SN_energy_per_current_crossing_time_list.append(SN_energy_per_current_crossing_time)
        SN_energy_per_final_crossing_time_list.append(SN_energy_per_final_crossing_time)
        (i) = (i + 1)

    global log_binding_energy_initial
    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(12, figsize=(3, 2.5))
        fig.add_subplot(1, 1, 1)
        # plt.loglog(time_axis, SN_energy_per_current_crossing_time_list, label='in a instant crossing time')
        # plt.loglog(time_axis, SN_energy_per_final_crossing_time_list, label='in a final crossing time')
        plt.loglog(time_axis, SNIa_energy_release_list, label='SNIa', ls='dashed', c='k')
        plt.loglog(time_axis, SNII_energy_release_list, label='SNII', c='k', lw=1.2)
        plt.loglog(time_axis, total_energy_release_list, ls="dotted", label='SNIa+SNII', c='k')
        # plt.loglog(time_axis, binding_energy_list, label='binding')
        plt.loglog([time_axis[0], time_axis[-1]], [10**log_binding_energy_initial, 10**log_binding_energy_initial], label='binding', lw=0.5, c='k')
        # plt.loglog(time_axis, total_gas_kinetic_energy_list, label='gas kinetic')
        plt.xlabel(r'time [yr]')
        plt.ylabel(r'Energy [erg]')
        plt.xlim(1e7, 1.1*1e10)
        # plt.ylim(5e49, 2e53)
        # plt.title('Energy produced by supernovae (within a crossing time)', fontsize=10)
        # plt.xlim(6, 1.01 * log_time_axis[-1])
        # plt.ylim(8.5, 11.6)
        plt.legend(prop={'size': 7})
        plt.tight_layout()
        if plot_save is True:
            plt.savefig('energy_evolution.pdf'.format(imf), dpi=250)

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/energy_evolution.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# time_axis\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(time_axis[i]))
        (i) = (i + 1)
    file.write("\n# SNIa_energy_release_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(SNIa_energy_release_list[i]))
        (i) = (i + 1)
    file.write("\n# SNII_energy_release_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(SNII_energy_release_list[i]))
        (i) = (i + 1)
    file.write("\n# SN_energy_per_current_crossing_time_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(SN_energy_per_current_crossing_time_list[i]))
        (i) = (i + 1)
    file.write("\n# SN_energy_per_final_crossing_time_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(SN_energy_per_final_crossing_time_list[i]))
        (i) = (i + 1)
    file.write("\n# total_energy_release_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(total_energy_release_list[i]))
        (i) = (i + 1)
    file.write("\n# binding_energy_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(binding_energy_list[i]))
        (i) = (i + 1)
    file.write("\n# total_gas_kinetic_energy_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(total_gas_kinetic_energy_list[i]))
        (i) = (i + 1)
    file.write("\n")
    file.close()

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/energy_mass.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# final SN_energy_release\n")
    file.write("{}\n".format(total_energy_release_list[-1]))
    file.write("# final binding_energy\n")
    file.write("{}\n".format(binding_energy_list[-1]))
    file.write("# final gas_kinetic_energy\n")
    file.write("{}\n".format(total_gas_kinetic_energy_list[-1]))
    file.close()

    # #
    # if plot_show is True or plot_save is True:
    #     plt.rc('font', family='serif')
    #     plt.rc('xtick', labelsize='x-small')
    #     plt.rc('ytick', labelsize='x-small')
    #     fig = plt.figure(13, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     # time_axis[0] = 1
    #     # time_axis_G = [0]*length_of_time_axis
    #     # for i in range(length_of_time_axis):
    #     #     time_axis_G[i] = time_axis[i]/10**9
    #     # gas_Z_over_X_list[i]=math.log(gas_Z_over_X_list[i], 10)
    #     plt.plot(log_time_axis, gas_Z_over_X_list, label='gas')
    #     plt.plot(log_time_axis, stellar_Z_over_X_list, label='stellar MW')
    #     plt.plot(log_time_axis, stellar_Z_over_X_list_luminosity_weighted, label='stellar LW')
    #     plt.plot([6, 11], [0, 0], color='red', ls='dashed', label='solar')
    #     # The [Z/X]s where the applied portinari98 stellar yield table will be changed for Z=0.0127, 0.008, 0.004, 0.0004.
    #     plt.plot([6, 11], [-1.173, -1.173], color='red', ls='dotted', lw=0.5)
    #     plt.plot([6, 11], [-0.523, -0.523], color='red', ls='dotted', lw=0.5)
    #     plt.plot([6, 11], [-0.272, -0.272], color='red', ls='dotted', lw=0.5)
    #     plt.xlabel(r'log(time [Gyr])')
    #     plt.ylabel('[Z/X]')
    #     plt.title('Metallicity evolution', fontsize=10)
    #     # plt.ylim(-2, 1)
    #     # if imf == 'igimf':
    #     #     plt.title('IGIMF')
    #     # elif imf == 'Kroupa':
    #     #     plt.title('Kroupa IMF')
    #     # plt.legend(scatterpoints=1, numpoints=1, loc=0, prop={'size': 7.5}, ncol=2)
    #     # plt.xlim(6.4, 1.01*time_axis[-1])
    #     # plt.ylim(-0.4, 0.2)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()
    #     if plot_save is True:
    #         plt.savefig('galaxy_evolution_fig_Z_{}.pdf'.format(imf), dpi=250)

    for i in range(length_of_time_axis):
        time_axis[i] = math.log(time_axis[i], 10)

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/Z_over_X_time.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# log_time_axis\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(log_time_axis[i]))
        (i) = (i + 1)
    file.write("\n# gas_Z_over_X_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(gas_Z_over_X_list[i]))
        (i) = (i + 1)
    file.write("\n# stellar_Z_over_X_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(stellar_Z_over_X_list[i]))
        (i) = (i + 1)
    file.write("\n# stellar_Z_over_X_list_luminosity_weighted\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(stellar_Z_over_X_list_luminosity_weighted[i]))
        (i) = (i + 1)
    file.write("\n")
    file.close()

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/Z_over_X_mass.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# gas_Z_over_X\n")
    file.write("{} ".format(gas_Z_over_X_list[-1]))
    file.write("\n# stellar_Z_over_X_list\n")
    file.write("{} ".format(stellar_Z_over_X_list[-1]))
    file.write("\n# stellar_Z_over_X_list_luminosity_weighted\n")
    file.write("{} ".format(stellar_Z_over_X_list_luminosity_weighted[-1]))
    file.write("\n")
    file.close()

    #
    global BH_mass_list, NS_mass_list, WD_mass_list, total_gas_mass_list, ejected_gas_mass_list
    for i in range(length_of_time_axis):
        if remnant_mass_list[i] < 10 ** (-10):
            remnant_mass_list[i] = 10 ** (-10)
        remnant_mass_list[i] = math.log(remnant_mass_list[i], 10)
        if total_gas_mass_list[i] < 10 ** (-10):
            total_gas_mass_list[i] = 10 ** (-10)
        total_gas_mass_list[i] = math.log(total_gas_mass_list[i], 10)
        if stellar_mass_list[i] < 10 ** (-10):
            stellar_mass_list[i] = 10 ** (-10)
        stellar_mass_list[i] = math.log(stellar_mass_list[i], 10)
        ejected_gas_mass_list[i] = math.log(ejected_gas_mass_list[i], 10)
        if WD_mass_list[i] < 10 ** (-10):
            WD_mass_list[i] = 10 ** (-10)
        WD_mass_list[i] = math.log(WD_mass_list[i], 10)
        if NS_mass_list[i] < 10 ** (-10):
            NS_mass_list[i] = 10 ** (-10)
        NS_mass_list[i] = math.log(NS_mass_list[i], 10)
        if BH_mass_list[i] < 10 ** (-10):
            BH_mass_list[i] = 10 ** (-10)
        BH_mass_list[i] = math.log(BH_mass_list[i], 10)
    # time_axis[0] = time_axis[1]
    if plot_show is True or plot_save is True:
        fig = plt.figure(14, figsize=(3, 2.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(time_axis, total_gas_mass_list, lw=1.5, label='gas', ls='dotted', c='k')
        # plt.plot(time_axis, ejected_gas_mass_list, lw=2, label='ejected gas')
        plt.plot(time_axis, stellar_mass_list, lw=1.5, label='living stars', c='k')
        print('plot stellar_mass final', stellar_mass_list[-1])
        plt.plot(time_axis, remnant_mass_list, lw=1.5, label='stellar remnants', ls='dashed', c='k')
        # plt.plot(time_axis, BH_mass_list, lw=2, label='black holes')
        # plt.plot(time_axis, NS_mass_list, lw=2, label='neutron stars')
        # plt.plot(time_axis, WD_mass_list, lw=2, label='white dwarfs')
        plt.xlabel(r'log$_{10}$(time [yr])')
        plt.ylabel(r'log$_{10}$(Mass [$M_\odot$])')
        # plt.title('Mass evolution', fontsize=10)
        # if imf == 'igimf':
        #     plt.title('IGIMF')
        # elif imf == 'Kroupa':
        #     plt.title('Kroupa IMF')
        plt.legend(prop={'size': 7})
        # plt.xlim(6.4, 1.01 * time_axis[-1])
        plt.xlim(6.4, 10.1)
        # plt.ylim(0, 7)
        plt.tight_layout()
        if plot_save is True:
            plt.savefig('mass_evolution.pdf'.format(imf), dpi=250)

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/mass_evolution.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# time_axis\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(time_axis[i]))
        (i) = (i + 1)
    file.write("\n# total_gas_mass_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(total_gas_mass_list[i]))
        (i) = (i + 1)
    file.write("\n# ejected_gas_mass_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(ejected_gas_mass_list[i]))
        (i) = (i + 1)
    file.write("\n# stellar_mass_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(stellar_mass_list[i]))
        (i) = (i + 1)
    file.write("\n# remnant_mass_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(remnant_mass_list[i]))
        (i) = (i + 1)
    file.write("\n# BH_mass_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(BH_mass_list[i]))
        (i) = (i + 1)
    file.write("\n# NS_mass_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(NS_mass_list[i]))
        (i) = (i + 1)
    file.write("\n# WD_mass_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(WD_mass_list[i]))
        (i) = (i + 1)
    file.write("\n")
    file.close()

    final_alive_stellar_mass = stellar_mass_list[-1]
    final_remnant_stellar_mass = remnant_mass_list[-1]
    final_alive_and_remnant_stellar_mass = math.log((10 ** final_alive_stellar_mass + 10 ** final_remnant_stellar_mass),
                                                    10)
    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/mass_ratio.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# final alive stellar mass in log\n")
    file.write("{}\n".format(final_alive_stellar_mass))
    file.write("# final alive + remnant mass in log\n")
    file.write("{}\n".format(final_alive_and_remnant_stellar_mass))
    file.write("# (alive + remnant) / alive IN LOG\n")
    file.write("{}\n".format(final_alive_and_remnant_stellar_mass - final_alive_stellar_mass))
    file.close()

    global total_star_formed
    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/SN_number_mass.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# final SNIa_number per stellar mass formed\n")
    file.write("{}\n".format(SNIa_number_list[-1] / total_star_formed))
    # print("total SNIa number per solar mass of star formed:", SNIa_number_list[-1]/total_star_formed)
    file.write("# final SNII_number per stellar mass formed\n")
    file.write("{}\n".format(SNII_number_list[-1] / total_star_formed))
    file.close()


    # global expansion_factor_instantaneous_list, expansion_factor_adiabat_list
    # if plot_show is True or plot_save is True:
    #     fig = plt.figure(15, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(time_axis, expansion_factor_instantaneous_list, label='instantaneous')
    #     plt.plot(time_axis, expansion_factor_adiabat_list, label='slow')
    #     plt.plot(time_axis, expansion_factor_list, label='average')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel(r'expansion factor')
    #     plt.legend(prop={'size': 7})
    #     plt.title('Galaxy size evolution', fontsize=10)
    #     # plt.xlim(6.4, 1.01 * time_axis[-1])
    #     # plt.ylim(7.3, 12.2)
    #     # plt.ylim(6, 12)
    #     plt.tight_layout()

    file = open('simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/plots/expansion_factor.txt'.format(imf, STF, SFR, SFEN, log_Z_0), 'w')
    file.write("# time_axis\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(time_axis[i]))
        (i) = (i + 1)
    file.write("\n# expansion_factor_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(expansion_factor_list[i]))
        (i) = (i + 1)
    file.write("\n# expansion_factor_instantaneous_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(expansion_factor_instantaneous_list[i]))
        (i) = (i + 1)
    file.write("\n# expansion_factor_adiabatic_list\n")
    i = 0
    while i < length_of_time_axis:
        file.write("{} ".format(expansion_factor_adiabat_list[i]))
        (i) = (i + 1)
    file.write("\n")
    file.close()

    # global ejected_gas_Mg_over_Fe_list, instant_ejected_gas_Mg_over_Fe_list
    # if plot_show is True or plot_save is True:
    #     fig = plt.figure(16, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(time_axis, ejected_gas_Mg_over_Fe_list, label='total')
    #     plt.plot(time_axis, instant_ejected_gas_Mg_over_Fe_list, label='instant')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel(r'[Mg/Fe]')
    #     # plt.xlim(6.4, 1.01 * time_axis[-1])
    #     # plt.ylim(7.3, 12.2)
    #     # plt.ylim(6, 12)
    #     plt.legend(prop={'size': 7})
    #     plt.title('[Mg/Fe] of the ejected gas at different time', fontsize=10)
    #     plt.tight_layout()
    #
    # global ejected_metal_mass_list
    # if plot_show is True or plot_save is True:
    #     fig = plt.figure(17, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    #     plt.plot(time_axis, ejected_metal_mass_list, label='total')
    #     plt.xlabel(r'log$_{10}$(time [yr])')
    #     plt.ylabel(r'ejected metal mass')
    #     # plt.xlim(6.4, 1.01 * time_axis[-1])
    #     # plt.ylim(7.3, 12.2)
    #     # plt.ylim(6, 12)
    #     plt.title('IMF averaged yield at different time', fontsize=10)
    #     plt.legend(prop={'size': 7})
    #     plt.tight_layout()

        #
        # global ejected_O_mass_till_this_time_tot_list, ejected_O_mass_till_this_time_SNIa_list, ejected_O_mass_till_this_time_SNII_list
        # if plot_show is True or plot_save is True:
        #     plt.rc('font', family='serif')
        #     plt.rc('xtick', labelsize='x-small')
        #     plt.rc('ytick', labelsize='x-small')
        #     fig = plt.figure(21, figsize=(3, 2.5))
        #     fig.add_subplot(1, 1, 1)
        #     plt.plot(log_time_axis, ejected_O_mass_till_this_time_tot_list, label='tot')
        #     plt.plot(log_time_axis, ejected_O_mass_till_this_time_SNIa_list, label='from SNIa')
        #     plt.plot(log_time_axis, ejected_O_mass_till_this_time_SNII_list, label='from SNII')
        #     plt.xlabel(r'log$_{10}$(time [yr])')
        #     plt.ylabel(r'ejected O [$M_\odot$]')
        #     plt.title('IMF averaged yield', fontsize=10)
        #     plt.legend(prop={'size': 7})
        #     plt.tight_layout()
        #
        # global ejected_Mg_mass_till_this_time_tot_list, ejected_Mg_mass_till_this_time_SNIa_list, ejected_Mg_mass_till_this_time_SNII_list
        # if plot_show is True or plot_save is True:
        #     plt.rc('font', family='serif')
        #     plt.rc('xtick', labelsize='x-small')
        #     plt.rc('ytick', labelsize='x-small')
        #     fig = plt.figure(22, figsize=(3, 2.5))
        #     fig.add_subplot(1, 1, 1)
        #     plt.plot(log_time_axis, ejected_Mg_mass_till_this_time_tot_list, label='tot')
        #     plt.plot(log_time_axis, ejected_Mg_mass_till_this_time_SNIa_list, label='from SNIa')
        #     plt.plot(log_time_axis, ejected_Mg_mass_till_this_time_SNII_list, label='from SNII')
        #     plt.xlabel(r'log$_{10}$(time [yr])')
        #     plt.ylabel(r'ejected Mg [$M_\odot$]')
        #     plt.title('IMF averaged yield from different type of SN', fontsize=10)
        #     plt.legend(prop={'size': 7})
        #     plt.tight_layout()
        #
        # global ejected_Fe_mass_till_this_time_tot_list, ejected_Fe_mass_till_this_time_SNIa_list, ejected_Fe_mass_till_this_time_SNII_list
        # if plot_show is True or plot_save is True:
        #     plt.rc('font', family='serif')
        #     plt.rc('xtick', labelsize='x-small')
        #     plt.rc('ytick', labelsize='x-small')
        #     fig = plt.figure(23, figsize=(3, 2.5))
        #     fig.add_subplot(1, 1, 1)
        #     plt.plot(log_time_axis, ejected_Fe_mass_till_this_time_tot_list, label='tot')
        #     plt.plot(log_time_axis, ejected_Fe_mass_till_this_time_SNIa_list, label='from SNIa')
        #     plt.plot(log_time_axis, ejected_Fe_mass_till_this_time_SNII_list, label='from SNII')
        #     plt.xlabel(r'log$_{10}$(time [yr])')
        #     plt.ylabel(r'ejected Fe [$M_\odot$]')
        #     plt.title('IMF averaged yield from different type of SN', fontsize=10)
        #     plt.legend(prop={'size': 7})
        #     plt.tight_layout()
        #     if plot_save is True:
        #         plt.savefig('Fe_production.pdf', dpi=250)
        #
        # global ejected_Ca_mass_till_this_time_tot_list, ejected_Ca_mass_till_this_time_SNIa_list, ejected_Ca_mass_till_this_time_SNII_list
        # if plot_show is True or plot_save is True:
        #     plt.rc('font', family='serif')
        #     plt.rc('xtick', labelsize='x-small')
        #     plt.rc('ytick', labelsize='x-small')
        #     fig = plt.figure(24, figsize=(3, 2.5))
        #     fig.add_subplot(1, 1, 1)
        #     plt.plot(log_time_axis, ejected_Ca_mass_till_this_time_tot_list, label='tot')
        #     plt.plot(log_time_axis, ejected_Ca_mass_till_this_time_SNIa_list, label='from SNIa')
        #     plt.plot(log_time_axis, ejected_Ca_mass_till_this_time_SNII_list, label='from SNII')
        #     plt.xlabel(r'log$_{10}$(time [yr])')
        #     plt.ylabel(r'ejected Ca [$M_\odot$]')
        #     plt.title('IMF averaged yield from different type of SN', fontsize=10)
        #     plt.legend(prop={'size': 7})
        #     plt.tight_layout()
        #     if plot_save is True:
        #         plt.savefig('Ca_production.pdf', dpi=250)
        #
        # global ejected_S_mass_till_this_time_tot_list, ejected_S_mass_till_this_time_SNIa_list, ejected_S_mass_till_this_time_SNII_list
        # if plot_show is True or plot_save is True:
        #     plt.rc('font', family='serif')
        #     plt.rc('xtick', labelsize='x-small')
        #     plt.rc('ytick', labelsize='x-small')
        #     fig = plt.figure(25, figsize=(3, 2.5))
        #     fig.add_subplot(1, 1, 1)
        #     plt.plot(log_time_axis, ejected_S_mass_till_this_time_tot_list, label='tot')
        #     plt.plot(log_time_axis, ejected_S_mass_till_this_time_SNIa_list, label='from SNIa')
        #     plt.plot(log_time_axis, ejected_S_mass_till_this_time_SNII_list, label='from SNII')
        #     plt.xlabel(r'log$_{10}$(time [yr])')
        #     plt.ylabel(r'ejected S [$M_\odot$]')
        #     plt.title('IMF averaged yield from different type of SN', fontsize=10)
        #     plt.legend(prop={'size': 7})
        #     plt.tight_layout()
        #     if plot_save is True:
        #         plt.savefig('S_production.pdf', dpi=250)
        #
        # global ejected_Si_mass_till_this_time_tot_list, ejected_Si_mass_till_this_time_SNIa_list, ejected_Si_mass_till_this_time_SNII_list
        # if plot_show is True or plot_save is True:
        #     plt.rc('font', family='serif')
        #     plt.rc('xtick', labelsize='x-small')
        #     plt.rc('ytick', labelsize='x-small')
        #     fig = plt.figure(26, figsize=(3, 2.5))
        #     fig.add_subplot(1, 1, 1)
        #     plt.plot(log_time_axis, ejected_Si_mass_till_this_time_tot_list, label='tot')
        #     plt.plot(log_time_axis, ejected_Si_mass_till_this_time_SNIa_list, label='from SNIa')
        #     plt.plot(log_time_axis, ejected_Si_mass_till_this_time_SNII_list, label='from SNII')
        #     plt.xlabel(r'log$_{10}$(time [yr])')
        #     plt.ylabel(r'ejected Si [$M_\odot$]')
        #     plt.title('IMF averaged yield from different type of SN', fontsize=10)
        #     plt.legend(prop={'size': 7})
        #     plt.tight_layout()
        #     if plot_save is True:
        #         plt.savefig('Si_production.pdf', dpi=250)
        #
        # global ejected_Ne_mass_till_this_time_tot_list, ejected_Ne_mass_till_this_time_SNIa_list, ejected_Ne_mass_till_this_time_SNII_list
        # if plot_show is True or plot_save is True:
        #     plt.rc('font', family='serif')
        #     plt.rc('xtick', labelsize='x-small')
        #     plt.rc('ytick', labelsize='x-small')
        #     fig = plt.figure(27, figsize=(3, 2.5))
        #     fig.add_subplot(1, 1, 1)
        #     plt.plot(log_time_axis, ejected_Ne_mass_till_this_time_tot_list, label='tot')
        #     plt.plot(log_time_axis, ejected_Ne_mass_till_this_time_SNIa_list, label='from SNIa')
        #     plt.plot(log_time_axis, ejected_Ne_mass_till_this_time_SNII_list, label='from SNII')
        #     plt.xlabel(r'log$_{10}$(time [yr])')
        #     plt.ylabel(r'ejected Ne [$M_\odot$]')
        #     plt.title('IMF averaged yield from different type of SN', fontsize=10)
        #     plt.legend(prop={'size': 7})
        #     plt.tight_layout()
        #     if plot_save is True:
        #         plt.savefig('Ne_production.pdf', dpi=250)

    if True: # plot_show is True:
        plt.show()
    return


def generate_SFH(distribution, Log_SFR, SFEN, sfr_tail, skewness, location):
    if distribution == "skewnorm":
        generate_sfh_skewnorm(Log_SFR, SFEN, sfr_tail, skewness, location)
    elif distribution == "flat":
        generate_sfh_flat(Log_SFR, SFEN)
    elif distribution == "flat_tail":
        generate_sfh_flat_tail(Log_SFR, SFEN)
    elif distribution == "lognorm":
        generate_sfh_lognorm(Log_SFR, SFEN)
    elif distribution == "given":
        generate_sfh_given()
    else:
        print('Warning: input unrecognized distribution name for galevo.generate_SFH')
    return


# def generate_sfh_given():
#     file = open('SFH.txt', 'w')
#     file.write("0\n")
#     j = 0
#     while j < 1:
#         file.write("0.0025\n")
#         (j) = (j + 1)
#     j = 1
#     while j < 10:
#         file.write("{}\n".format(0.0025 * 0.6 ** j))
#         (j) = (j + 1)
#     j = 0
#     while j < 1301 - SFEN:
#         file.write("0\n")
#         (j) = (j + 1)
#     file.write("# The value in each line stand for the SFR [solar mass / yr]\n")
#     file.write("# in a star formation epoch (10 Myr)\n")
#     file.write("# start from time 0 for the first line.\n")
#     file.write("# Warning! Effective line number must be larger than 1.\n")
#     file.write("# Add a '0' in the next line if there is only one line.\n")
#
#     file.close()
#     return

def generate_sfh_given():
    file = open('SFH.txt', 'w')
    file.write("0\n")
    file.write("0.0007\n")
    j = 1
    while j < 10:
        file.write("0.0009\n")
        (j) = (j + 1)
    j = 1
    while j < 90:
        file.write("{}\n".format(0.0009 * 0.9 ** j))
        (j) = (j + 1)
    j = 0
    while j < 1301 - SFEN:
        file.write("0\n")
        (j) = (j + 1)
    file.write("# The value in each line stand for the SFR [solar mass / yr]\n")
    file.write("# in a star formation epoch (10 Myr)\n")
    file.write("# start from time 0 for the first line.\n")
    file.write("# Warning! Effective line number must be larger than 1.\n")
    file.write("# Add a '0' in the next line if there is only one line.\n")

    file.close()
    return

# def generate_sfh_given():
#     file = open('SFH.txt', 'w')
#     file.write("0\n")
#     j = 0
#     while j < 67:
#         file.write("0.00005\n")
#         (j) = (j + 1)
#     j = 0
#     while j < 20:
#         file.write("{}\n".format(0.00005 * 0.9 ** j))
#         (j) = (j + 1)
#     j = 0
#     while j < 1301 - SFEN:
#         file.write("0\n")
#         (j) = (j + 1)
#     file.write("# The value in each line stand for the SFR [solar mass / yr]\n")
#     file.write("# in a star formation epoch (10 Myr)\n")
#     file.write("# start from time 0 for the first line.\n")
#     file.write("# Warning! Effective line number must be larger than 1.\n")
#     file.write("# Add a '0' in the next line if there is only one line.\n")
#
#     file.close()
#     return

def generate_sfh_flat(Log_SFR, SFEN):
    # Flat distribution for star formation history
    # took input: star formation rate, star formation event number
    file = open('SFH.txt', 'w')
    file.write("0\n")
    j = 0
    while j < SFEN:
        file.write("{}\n".format(10 ** Log_SFR))
        (j) = (j + 1)
    j = 0
    while j < 1301 - SFEN:
        file.write("0\n")
        (j) = (j + 1)
    file.write("# The value in each line stand for the SFR [solar mass / yr]\n")
    file.write("# in a star formation epoch (10 Myr)\n")
    file.write("# start from time 0 for the first line.\n")
    file.write("# Warning! Effective line number must be larger than 1.\n")
    file.write("# Add a '0' in the next line if there is only one line.\n")

    file.close()
    return


def generate_sfh_flat_tail(Log_SFR, SFEN):
    # Flat distribution for star formation history
    # took input: star formation rate, star formation event number
    file = open('SFH.txt', 'w')
    file.write("0\n")
    j = 0
    while j < SFEN:
        file.write("{}\n".format(10 ** Log_SFR))
        (j) = (j + 1)
    j = 0
    while j < 1389 - SFEN:
        file.write("0\n")
        (j) = (j + 1)
    j = 0
    while j < 10:
        file.write("{}\n".format(10 ** (Log_SFR-2)))
        (j) = (j + 1)
    file.write("# The value in each line stand for the SFR [solar mass / yr]\n")
    file.write("# in a star formation epoch (10 Myr)\n")
    file.write("# start from time 0 for the first line.\n")
    file.write("# Warning! Effective line number must be larger than 1.\n")
    file.write("# Add a '0' in the next line if there is only one line.\n")

    file.close()
    return


def generate_sfh_skewnorm(Log_SFR, SFEN, sfr_tail, skewness, location):
    tot_sf_set = 10 ** Log_SFR * SFEN
    tot_sf = 0
    while tot_sf < tot_sf_set:
        SFEN += 1
        result_cal_tot_sf = cal_tot_sf(Log_SFR, SFEN, skewness, location)
        (tot_sf) = (result_cal_tot_sf[0])

    file = open('SFH.txt', 'w')
    file.write("0\n")
    sfr_for_this_epoch = 0
    result_starburst_sf = 0
    result_tail_sf = 0
    j = 0
    while j < SFEN:
        sfr_for_this_epoch = result_cal_tot_sf[1] * result_cal_tot_sf[2][j]
        file.write("{}\n".format(sfr_for_this_epoch))
        if sfr_for_this_epoch > 10 ** Log_SFR / 2:
            result_starburst_sf += sfr_for_this_epoch
        else:
            result_tail_sf += sfr_for_this_epoch
        (j) = (j + 1)
    sfr_for_the_tail_epoch = sfr_for_this_epoch / 2
    if sfr_tail == 0:
        j = 0
        while j < 1301 - SFEN:
            file.write("0\n")
            (j) = (j + 1)
    elif sfr_tail == 1:
        j = 0
        while j < 101:
            file.write("{}\n".format(sfr_for_the_tail_epoch))
            result_tail_sf += sfr_for_the_tail_epoch
            (j) = (j + 1)
        while j < 1301 - SFEN:
            file.write("0\n")
            (j) = (j + 1)
    file.write("# The value in each line stand for the SFR [solar mass / yr]\n")
    file.write("# in a star formation epoch (10 Myr)\n")
    file.write("# start from time 0 for the first line.\n")
    file.write("# Warning! Effective line number must be larger than 1.\n")
    file.write("# Add a '0' in the next line if there is only one line.\n")

    if sfr_tail == 1:
        print("star formation tail (after the SFR is lower than half of the maximum value) contributes {}% "
              "of the total star formation.".format(
            round(result_tail_sf / (result_starburst_sf + result_tail_sf) * 100, 2)))

    file.close()
    return


def generate_sfh_lognorm(Log_SFR, SFEN):
    tot_sf_set = 10 ** Log_SFR * SFEN
    time_length_in_Gyr = 13
    time_step_number = time_length_in_Gyr * 100

    from scipy.stats import lognorm
    s = 1
    sc = SFEN / 2
    time_list = np.linspace(0, time_step_number, time_step_number)
    star_formation_rate = tot_sf_set * lognorm.pdf(time_list, s, scale=sc)

    i = 20
    while i < time_step_number:
        # if star_formation_rate[i] < star_formation_rate[1]/10:
        #     star_formation_rate[i] = 0
        star_formation_rate[i] = 0
        (i) = (i+1)

    file = open('SFH.txt', 'w')
    for i in range(time_step_number):
        file.write("{}\n".format(star_formation_rate[i]))

    file.write("# The value in each line stand for the SFR [solar mass / yr]\n")
    file.write("# in a star formation epoch (10 Myr)\n")
    file.write("# start from time 0 for the first line.\n")
    file.write("# Warning! Effective line number must be larger than 1.\n")
    file.write("# Add a '0' in the next line if there is only one line.\n")

    file.close()

    # plt.plot(time_list, star_formation_rate, label='lognorm SFH')
    # plt.xlabel('time step')
    # plt.ylabel(r'SFR [solar mass/year]')
    # plt.show()

    return


def cal_tot_sf(SFR, SFEN, skewness, location):
    # Skew normal distribution for star formation history
    # took input: maximum star formation rate, star formation event number
    # from scipy.stats import f
    from scipy.stats import skewnorm
    x = np.linspace(skewnorm.ppf(0.01, skewness, location, 1), skewnorm.ppf(0.999999999, skewness, location, 1), SFEN)
    y = skewnorm.pdf(x, skewness, location, 1)
    # skewnorm.pdf(x, a, loc, scale) is the location and scale parameters,
    #   [identically equivalent to skewnorm.pdf(y, a) / scale with y = (x - loc) / scale]
    # The scale is not used as the SFEN & SFR setup the scale through parameter tot_sf_set & mult.
    mult = 10 ** SFR / max(y)
    j = 0
    tot_sf = 0
    while j < SFEN:
        sf = mult * y[j]
        tot_sf += sf
        (j) = (j + 1)
    return tot_sf, mult, y


if __name__ == '__main__':
    SFEN = "None"

    ### Generate a new SFH.txt file according to the following given parameters ###

    SFEN = 3  # the number of the 10 Myr star formation epoch (thus 10 stand for a star formation timescale of 100 Myr)
    Log_SFR = -2.47  # logarithmic characteristic star formation rate
    location = 0  # SFH shape parameter
    skewness = 10  # SFH shape parameter
    sfr_tail = 0  # SFH shape parameter
    # generate_SFH("given", Log_SFR, SFEN, sfr_tail, skewness, location)
    ## input "flat", "lognorm", "skewnorm", or "given" to generate a boxy, lognormal, or skewnorm SFH, respectively.

    ####################################

    # str_yield_table='WW95' or 'portinari98' or 'Kobayashi06', specify the stellar evolution table
    # imf='igimf' or 'Kroupa' or 'Salpeter' or 'diet_Salpeter'...(see line 652 above), specify galaxy IMF model.
    # SFH_model='provided' or 'gas_mass_dependent' specify the star formation history.
    # The 'provided' SFH is given in SFH.txt;
    # The 'gas_mass_dependent' use SFH.txt to setup the initial condition
    # then recalculate SFR at each timestep, resulting a SFH similar to SFH.txt but gas mass dependent.
    # SNIa_yield_table='Thielemann1993' or 'Seitenzahl2013' or 'Iwamoto1999'
    # solar_abu_table='Anders1989' or 'Asplund2009'

    galaxy_evol(imf='igimf', STF=0.05, SFEN=SFEN, Z_0=0.015*1e-16, solar_mass_component="Asplund2009_mass",
                str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
                time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
                SFH_model='provided', SFE=0.004, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
                solar_abu_table='Asplund2009',
                high_time_resolution=None, plot_show=None, plot_save=None, outflow=100, check_igimf=None)
    
    # Use plot_show=True on persenal computer to view the simualtion result immidiately after the computation
    # Use plot_show=None if running on a computer cluster to avoid possible issues.
    # In both cases, the simulation results are saved as txt files.

    # galaxy_evol(imf='Salpeter', STF=0.039, SFEN=SFEN, Z_0=0.02 * 1e-1, solar_mass_component="Asplund2009_mass",
    #             str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
    #             time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
    #             SFH_model='provided', SFE=0.0015, SNIa_ON='SD', SNIa_yield_table='Iwamoto1999',
    #             solar_abu_table='Asplund2009',
    #             high_time_resolution=True, plot_show=None, plot_save=None, outflow=100, check_igimf=None)




    # # Model Reproduce L19-Salpeter
    # # Change line 1128 53.7 to 53.75
    # # # Reproducing the Salpeter IMF model:
    # #
    # galaxy_evol(imf='Salpeter', STF=0.01688, SFEN=SFEN, Z_0=0.02*1e-16, solar_mass_component="Asplund2009_mass",
    #             str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
    #             time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
    #             SFH_model='gas_mass_dependent', SFE=0.001, SNIa_ON='SD', SNIa_yield_table='Iwamoto1999',
    #             solar_abu_table='Asplund2009',
    #             high_time_resolution=None, plot_show=True, plot_save=True, outflow=100, check_igimf=None)


    # # # Model Reproduce L19-R14
    # # # Change line 1128 53.68 STF=0.0194
    # # # # Reproducing the IGIMF model:
    # # #
    # galaxy_evol(imf='igimf', STF=0.01688, SFEN=SFEN, Z_0=0.02*1e-7, solar_mass_component="Asplund2009_mass",
    #             str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
    #             time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
    #             SFH_model='gas_mass_dependent', SFE=0.00082, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
    #             solar_abu_table='Asplund2009',
    #             high_time_resolution=None, plot_show=True, plot_save=True, outflow=105, check_igimf=None)


    # # Model IGIMF-R14
    # # R14 change galimf.py line 1050, 1053, 1063;   1166, 1169, 1173
    # # SFH with 0.7, 0.9... to line 48
    # # line 1801
    # galaxy_evol(imf='igimf', STF=0.042, SFEN=SFEN, Z_0=0.02 * 1e-7, solar_mass_component="Asplund2009_mass",
    #             str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
    #             time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
    #             SFH_model='gas_mass_dependent', SFE=0.0017, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
    #             solar_abu_table='Asplund2009',
    #             high_time_resolution=None, plot_show=True, plot_save=True, outflow=100, check_igimf=None)


    # # #
    # # Model IGIMF-R14-SD with the Greggio 1983 SNIa rate
    # galaxy_evol(imf='igimf', STF=0.039, SFEN=SFEN, Z_0=0.02 * 1e-7, solar_mass_component="Asplund2009_mass",
    #             str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
    #             time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
    #             SFH_model='gas_mass_dependent', SFE=0.0015, SNIa_ON='SD', SNIa_yield_table='Iwamoto1999',
    #             solar_abu_table='Asplund2009',
    #             high_time_resolution=None, plot_show=True, plot_save=None, outflow=100, check_igimf=None)


    # # Model ?
    # # # with IGIMF3: from R14 to IGIMF2 change galimf.py line 1050, 1053, 1063; 1166, 1168, 1169, 1172, 1173
    # galaxy_evol(imf='igimf', STF=0.08, SFEN=SFEN, Z_0=0.02*1e-7, solar_mass_component="Asplund2009_mass",
    #             str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
    #             time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
    #             SFH_model='gas_mass_dependent', SFE=0.004, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
    #             solar_abu_table='Asplund2009',
    #             high_time_resolution=None, plot_show=True, plot_save=True, outflow=100, check_igimf=None)



    # # Model with different SFH?
    # galaxy_evol(imf='igimf', STF=0.0218, SFEN=SFEN, Z_0=0.02*1e-16, solar_mass_component="Asplund2009_mass",
    #             str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
    #             time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
    #             SFH_model='gas_mass_dependent', SFE=0.004, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
    #             solar_abu_table='Asplund2009',
    #             high_time_resolution=None, plot_show=True, plot_save=True, outflow=100, check_igimf=None)


    # # Model IGIMF2
    # galaxy_evol(imf='igimf', STF=0.073, SFEN=SFEN, Z_0=0.02 * 1e-7, solar_mass_component="Asplund2009_mass",
    #             str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
    #             time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
    #             SFH_model='gas_mass_dependent', SFE=0.0045, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
    #             solar_abu_table='Asplund2009',
    #             high_time_resolution=None, plot_show=True, plot_save=True, outflow=100, check_igimf=None)


    # # Model IGIMF3
    # galaxy_evol(imf='igimf', STF=0.03, SFEN=SFEN, Z_0=0.02 * 1e-7, solar_mass_component="Asplund2009_mass",
    #             str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
    #             time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
    #             SFH_model='gas_mass_dependent', SFE=0.002, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
    #             solar_abu_table='Asplund2009',
    #             high_time_resolution=None, plot_show=True, plot_save=True, outflow=100, check_igimf=None)


    # # Model IGIMF4
    # galaxy_evol(imf='igimf', STF=0.052, SFEN=SFEN, Z_0=0.02 * 1e-7, solar_mass_component="Asplund2009_mass",
    #             str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
    #             time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
    #             SFH_model='gas_mass_dependent', SFE=0.0035, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
    #             solar_abu_table='Asplund2009',
    #             high_time_resolution=None, plot_show=True, plot_save=True, outflow=100, check_igimf=None)

    # # Model IGIMF-Z
    # galaxy_evol(imf='igimf', STF=0.043, SFEN=SFEN, Z_0=0.02 * 1e-7, solar_mass_component="Asplund2009_mass",
    #             str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
    #             time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
    #             SFH_model='gas_mass_dependent', SFE=0.0035, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
    #             solar_abu_table='Asplund2009',
    #             high_time_resolution=None, plot_show=True, plot_save=True, outflow=100, check_igimf=None)
