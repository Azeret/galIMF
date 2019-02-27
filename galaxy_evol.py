import numpy as np
import time
import math
import importlib
from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import gc
import sys
sys.path.insert(0, 'Generated_IGIMFs')
sys.path.insert(0, 'IMFs')
import element_weight_table, solar_element_abundances
from IMFs import Kroupa_IMF, diet_Salpeter_IMF


def galaxy_evol(imf='igimf', unit_SFR=1, SFE=0.3, SFEN=1, Z_0=0.000000134, Z_solar=0.01886, str_evo_table='portinari98',
                IMF_name='Kroupa', steller_mass_upper_bound=150,
                time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
                SNIa_ON=True, high_time_resolution=True, plot_show=True, plot_save=None, outflow=None, check_igimf=False):

    start_time = time.time()

    ######################
    # If imf='igimf', the model will use variable IMF, imf='Kroupa' will use Kroupa IMF
    # unit_SFR correspond to SFH.txt. A 1 in SFH.txt stand for SFR = 1 * unit_SFR [solar mass/year] in a 10 Myr epoch.
    # SFE is the total star formation efficency (total stellar mass/total gas mass in 13Gyr), which determines the initial gas mass.
    # Z_0 is the initial metallicity
    ######################
    global igimf_mass_function, mass_grid_table, mass_grid_table2, Mfinal_table, Mmetal_table, M_element_table
    global time_axis, gas_Z_over_X_list, total_star_formed, \
        O_over_H_list, Mg_over_H_list, Fe_over_H_list, Y_list, \
        Mg_over_Fe_list, C_over_Fe_list, N_over_Fe_list, Ca_over_Fe_list, O_over_Fe_list, \
        stellar_O_over_H_list, stellar_Mg_over_H_list, stellar_Y_list, stellar_Fe_over_H_list, \
        stellar_Mg_over_Fe_list, stellar_C_over_Fe_list, stellar_N_over_Fe_list, stellar_Ca_over_Fe_list, \
        stellar_O_over_Fe_list, stellar_Z_over_X_list, stellar_Z_over_H_list, \
        stellar_O_over_H_list_luminosity_weighted, stellar_Mg_over_H_list_luminosity_weighted, stellar_Y_list_luminosity_weighted, stellar_Fe_over_H_list_luminosity_weighted, \
        stellar_Mg_over_Fe_list_luminosity_weighted, stellar_C_over_Fe_list_luminosity_weighted, stellar_N_over_Fe_list_luminosity_weighted, stellar_Ca_over_Fe_list_luminosity_weighted, \
        stellar_O_over_Fe_list_luminosity_weighted, stellar_Z_over_X_list_luminosity_weighted, stellar_Z_over_H_list_luminosity_weighted, \
        remnant_mass_list, total_gas_mass_list, stellar_mass_list, ejected_gas_mass_list, ejected_metal_mass_list, \
        ejected_gas_Mg_over_Fe_list, instant_ejected_gas_Mg_over_Fe_list, expansion_factor_list, \
        expansion_factor_instantaneous_list, expansion_factor_adiabat_list
    global SNIa_energy_release_list, SNIa_number_list, SNIa_number_per_century, \
        SNII_energy_release_list, SNII_number_list, SNII_number_per_century, \
        total_energy_release_list, SN_number_per_century, total_gas_kinetic_energy_list, original_gas_mass#, binding_energy_list
    global BH_mass_list, NS_mass_list, WD_mass_list
    global all_sf_imf, all_sfr
    global times_calculate_igimf
    times_calculate_igimf = 0
    ###################
    ### preparation ###
    ###################

    # get all avaliable metallicity from stellar evolution table
    (Z_list, Z_list_2, Z_list_3) = function_get_avaliable_Z(str_evo_table)

    # read in SFH
    SFH_input = np.loadtxt('SFH.txt')
    length_list_SFH_input = len(SFH_input)

    i = 0
    total_SF = 0
    while i < length_list_SFH_input:
        total_SF += SFH_input[i]
        (i) = (i + 1)

    # Star_formation_efficiency = SFE
    total_star_formed = unit_SFR * 10**7 * total_SF
    original_gas_mass = total_star_formed / SFE  # in solar mass unit

    # Create the time steps (x axis) for final output
    time_axis = []
    time_resolution = time_resolution_in_Myr * 10 ** 5 * 10
    if True:  # high_time_resolution==True:
        # for i in range(10 ** 7, 10 ** 8, time_resolution * 10):
        #     time_axis += [i]
        # for i in range(10 ** 8, 10 ** 9, time_resolution * 100):
        #     time_axis += [i]
        for i in range(10 ** 9, 15 * 10 ** 9, time_resolution * 1000):
            time_axis += [i]
    else:
        plot_at_age = [5 * 10 ** 7, 1 * 10 ** 8, 5 * 10 ** 8, 1 * 10 ** 9, 9 * 10 ** 9, 10 * 10 ** 9, 11 * 10 ** 9]
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
                time_axis_for_SFH_input += [i * 10 ** 7]
                time_axis_for_SFH_input += [i * 10 ** 7 + 5 * 10 ** 5]
                time_axis_for_SFH_input += [i * 10 ** 7 + 1 * 10 ** 6]
                time_axis_for_SFH_input += [i * 10 ** 7 + 2 * 10 ** 6]
                time_axis_for_SFH_input += [i * 10 ** 7 + 4 * 10 ** 6]
                time_axis_for_SFH_input += [i * 10 ** 7 + 6 * 10 ** 6]
                time_axis_for_SFH_input += [i * 10 ** 7 + 8 * 10 ** 6]
                time_axis_for_SFH_input += [i * 10 ** 7 + 10 * 10 ** 6]
            else:
                time_axis_for_SFH_input += [i * 10 ** 7]
                time_axis_for_SFH_input += [i * 10 ** 7 + 9 * 10 ** 6]
        (i)=(i+1)

    # the final time axis is the sorted combination of the two
    time_axis = sorted(list(set(time_axis + time_axis_for_SFH_input)))
    print("calculate at time:", time_axis)
    length_list_time_step = len(time_axis)

    ###################
    ###  main loop  ###
    ###################
    S_F_R_of_this_epoch_list = []
    S_F_F = 1
    # define an array save SF event informations that will be used in every latter time steps
    all_sf_imf = []
    all_sfr = []
    epoch_info = [] # save [S_F_R_of_this_epoch, M_tot_of_this_epoch, igimf_mass_function, igimf_normalization]
    BH_mass_list = []
    NS_mass_list = []
    WD_mass_list = []
    gas_Z_over_X_list = []
    O_over_H_list = []
    Y_list = []
    Mg_over_H_list = []
    Fe_over_H_list = []
    Mg_over_Fe_list = []
    C_over_Fe_list = []
    N_over_Fe_list = []
    Ca_over_Fe_list = []
    O_over_Fe_list = []
    stellar_O_over_H_list = []
    stellar_Mg_over_H_list = []
    stellar_Y_list = []
    stellar_Fe_over_H_list = []
    stellar_Mg_over_Fe_list = []
    stellar_C_over_Fe_list = []
    stellar_N_over_Fe_list = []
    stellar_Ca_over_Fe_list = []
    stellar_O_over_Fe_list = []
    stellar_Z_over_X_list = []
    stellar_Z_over_H_list = []
    stellar_O_over_H_list_luminosity_weighted = []
    stellar_Mg_over_H_list_luminosity_weighted = []
    stellar_Y_list_luminosity_weighted = []
    stellar_Fe_over_H_list_luminosity_weighted = []
    stellar_Mg_over_Fe_list_luminosity_weighted = []
    stellar_C_over_Fe_list_luminosity_weighted = []
    stellar_N_over_Fe_list_luminosity_weighted = []
    stellar_Ca_over_Fe_list_luminosity_weighted = []
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
    Z_over_X = math.log(Z_0 / Z_solar, 10)
    # do calculation for each time start from time 0
    time_step = 0
    gc_collect_check = 1
    # do calculation for each time to the end time
    while time_step < length_list_time_step:
        # get time
        this_time = time_axis[time_step]
        # calculated the array index (line number in SFH.txt) this_time has reached
        epoch_index_limit = (this_time+1) / 10**7
        if epoch_index_limit > length_list_SFH_input:
            epoch_index_limit = length_list_SFH_input
        last_time_age = 0
        age_of_this_epoch = 0
        number_in_SNIa_boundary = 0
        # get masses and metallicity at this time (values are calculated by the end of last time step)
        # initialize values
        total_energy_release = 0
        SNIa_energy_release = 0
        SNIa_number = 0
        SNII_energy_release = 0
        SNII_number = 0
        if time_step == 0:
            eject_H_mass = 0
            eject_C_mass = 0
            eject_N_mass = 0
            eject_O_mass = 0
            eject_Mg_mass = 0
            eject_Ca_mass = 0
            eject_Fe_mass = 0
            eject_metal_mass = 0

            total_gas_mass_at_this_time = 0
            ejected_gas_mass_at_this_time = 0
            ejected_metal_mass_at_last_time = 0

            M_tot_of_last_time = 0
            M_tot_of_this_time = 0
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
            total_Fe_mass_at_this_time = 0

            total_H_mass_at_last_time = original_gas_mass * (0.7381-Z_0)
            H_weight = element_weight_table.function_element_weight("H")
            total_He_mass_at_last_time = original_gas_mass * 0.2485
            total_C_mass_at_last_time = original_gas_mass * 0.7381 / H_weight\
                                        * 10**(solar_element_abundances.function_solar_element_abundances("C")-12) \
                                        * element_weight_table.function_element_weight("C") * Z_0 / Z_solar
            total_N_mass_at_last_time = original_gas_mass * 0.7381/H_weight\
                                        * 10**(solar_element_abundances.function_solar_element_abundances("N")-12) \
                                        * element_weight_table.function_element_weight("N") * Z_0 / Z_solar
            total_O_mass_at_last_time = original_gas_mass * 0.7381/H_weight\
                                        * 10**(solar_element_abundances.function_solar_element_abundances("O")-12) \
                                        * element_weight_table.function_element_weight("O") * Z_0 / Z_solar
            total_Mg_mass_at_last_time = original_gas_mass * 0.7381/H_weight\
                                        * 10**(solar_element_abundances.function_solar_element_abundances("Mg")-12) \
                                        * element_weight_table.function_element_weight("Mg") * Z_0 / Z_solar
            total_Ca_mass_at_last_time = original_gas_mass * 0.7381/H_weight\
                                        * 10**(solar_element_abundances.function_solar_element_abundances("Ca")-12) \
                                        * element_weight_table.function_element_weight("Ca") * Z_0 / Z_solar
            total_Fe_mass_at_last_time = original_gas_mass * 0.7381/H_weight\
                                        * 10**(solar_element_abundances.function_solar_element_abundances("Fe")-12) \
                                        * element_weight_table.function_element_weight("Fe") * Z_0 / Z_solar
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

            stellar_metal_luminosity_at_this_time = 0
            stellar_H_luminosity_at_this_time = 0
            stellar_He_luminosity_at_this_time = 0
            stellar_C_luminosity_at_this_time = 0
            stellar_N_luminosity_at_this_time = 0
            stellar_O_luminosity_at_this_time = 0
            stellar_Mg_luminosity_at_this_time = 0
            stellar_Ca_luminosity_at_this_time = 0
            stellar_Fe_luminosity_at_this_time = 0

            metal_mass_fraction_in_gas = [Z_gas_this_time_step, 0.7381-Z_0, 0.2485,
                                          total_C_mass_at_last_time/original_gas_mass,
                            total_N_mass_at_last_time/original_gas_mass, total_O_mass_at_last_time/original_gas_mass,
                            total_Mg_mass_at_last_time/original_gas_mass, total_Ca_mass_at_last_time/original_gas_mass,
                            total_Fe_mass_at_last_time/original_gas_mass] # H He C N O Mg Ca Fe
        else:
            total_gas_mass_at_last_time = total_gas_mass_at_this_time
            total_gas_mass_at_this_time = 0
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
            total_Fe_mass_at_last_time = total_Fe_mass_at_this_time
            total_Fe_mass_at_this_time = 0
            M_tot_of_last_time = M_tot_of_this_time
            M_tot_of_this_time = 0
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
            ejected_Fe_mass_till_this_time = 0
            ejected_metal_mass_at_last_time = ejected_metal_mass_at_this_time
            # Fe_H_mass_ratio_at_last_time = Fe_H_mass_ratio_at_this_time
            Z_gas_this_time_step = total_metal_mass_in_gas_at_last_time/total_gas_mass_at_last_time
            metal_mass_fraction_in_gas = [Z_gas_this_time_step,
                                          total_H_mass_at_last_time/total_gas_mass_at_last_time,
                                          total_He_mass_at_last_time/total_gas_mass_at_last_time,
                                          total_C_mass_at_last_time/total_gas_mass_at_last_time,
                                          total_N_mass_at_last_time/total_gas_mass_at_last_time,
                                          total_O_mass_at_last_time/total_gas_mass_at_last_time,
                                          total_Mg_mass_at_last_time/total_gas_mass_at_last_time,
                                          total_Ca_mass_at_last_time/total_gas_mass_at_last_time,
                                          total_Fe_mass_at_last_time/total_gas_mass_at_last_time]
            stellar_metal_mass_at_this_time = 0
            stellar_H_mass_at_this_time = 0
            stellar_He_mass_at_this_time = 0
            stellar_C_mass_at_this_time = 0
            stellar_N_mass_at_this_time = 0
            stellar_O_mass_at_this_time = 0
            stellar_Mg_mass_at_this_time = 0
            stellar_Ca_mass_at_this_time = 0
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
        # add up metals contributed by SSP from each SF epoch
        # consider only the SF event (epoch) that had happend
        Fe_production_SNII = 0
        Fe_production_SNIa = 0
        Mg_production_SNII = 0
        Mg_production_SNIa = 0
        O_production_SNII = 0
        O_production_SNIa = 0

        epoch_index = 0
        while epoch_index < epoch_index_limit:
            # get age
            age_of_this_epoch = this_time - epoch_index * 10 ** 7
            # get SFR, M_tot, igimf, integrated igimf, stellar lifetime and stellar remnant mass for this metallicity
            if epoch_index == len(epoch_info):
                # SFR
                S_F_R_of_this_epoch = SFH_input[epoch_index] * unit_SFR

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
                        igimf_of_this_epoch = function_get_igimf_for_this_epoch(S_F_R_of_this_epoch, Z_over_X, this_time, epoch_index, check_igimf)  # Fe_over_H_number_ratio)
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
                    age_of_this_epoch_at_end = (length_list_SFH_input - epoch_index -1 ) * 10 ** 7
                    #
                    def igimf_xi_function(mass):
                        return igimf_of_this_epoch.custom_imf(mass, this_time)

                    def igimf_mass_function(mass):
                        return igimf_of_this_epoch.custom_imf(mass, this_time) * mass

                    def igimf_luminous_function(mass):
                        if mass < 0.23 ** (1 / 1.7):
                            return igimf_of_this_epoch.custom_imf(mass, this_time) * 0.23 * mass ** 2.3
                        elif mass < 1.96:
                            return igimf_of_this_epoch.custom_imf(mass, this_time) * mass ** 4
                        elif mass < (32000 / 1.4) ** (1 / 2.5):
                            return igimf_of_this_epoch.custom_imf(mass, this_time) * 1.4 * mass ** 3.5
                        else:
                            return igimf_of_this_epoch.custom_imf(mass, this_time) * 32000 * mass

                    # integrated igimf_mass_function from 0.08 to steller_mass_upper_bound
                    integrate_igimf_mass = quad(igimf_mass_function, 0.08, steller_mass_upper_bound, limit=40)[0]
                    # as the integration of the IGIMF always has a small (at least for low SFRs) computational error,
                    # it need to be fixed by mutiplying a calibration factor which is close to 1:
                    mass_calibration_factor = M_tot_of_this_epoch / integrate_igimf_mass
                    # print("mass_calibration_factor:", mass_calibration_factor)
                    
                    # integrate_igimf_mass_l = quad(igimf_mass_function, 0.08, 3, limit=40)[0]
                    # integrate_igimf_mass_h = quad(igimf_mass_function, 8, steller_mass_upper_bound, limit=40)[0]
                    # integrate_igimf_mass_m = quad(igimf_mass_function, 1.5, 8, limit=40)[0]
                    # print("high mass star mass ratio:", integrate_igimf_mass_h/integrate_igimf_mass)
                    # print("middle mass star mass ratio:", integrate_igimf_mass_m/integrate_igimf_mass)
                    # print("Low mass star mass ratio:", integrate_igimf_mass_l/integrate_igimf_mass)
                    # integrate_igimf_number = quad(igimf_xi_function, 0.08, steller_mass_upper_bound, limit=40)[0]
                    # integrate_igimf_number_l = quad(igimf_xi_function, 0.08, 3, limit=40)[0]
                    # integrate_igimf_number_h = quad(igimf_xi_function, 8, steller_mass_upper_bound, limit=40)[0]
                    # integrate_igimf_number_m = quad(igimf_xi_function, 1.5, 8, limit=40)[0]
                    # print("high mass star number ratio:", integrate_igimf_number_h/integrate_igimf_number)
                    # print("middle mass star number ratio:", integrate_igimf_number_m/integrate_igimf_number)
                    # print("Low mass star number ratio:", integrate_igimf_number_l/integrate_igimf_number)

                    # Choose the closest metallicity
                    Z_select_in_table = function_select_metal(Z_gas_this_time_step, Z_list)
                    Z_select_in_table_2 = function_select_metal(Z_gas_this_time_step, Z_list_2)
                    if str_evo_table != "portinari98":
                        Z_select_in_table_3 = function_select_metal(Z_gas_this_time_step, Z_list_3)
                    else:
                        Z_select_in_table_3 = None
                    # read in interpolated stellar lifetime table
                    (mass_1, mass, lifetime_table) = function_read_lifetime(str_evo_table, Z_select_in_table)
                    # read in interpolated stellar final mass
                    (mass_12, Mfinal_table) = function_read_Mfinal(str_evo_table, Z_select_in_table)
                    # read in interpolated stellar ejected metal mass
                    (mass_2, mass2, Mmetal_table) = function_read_Mmetal(str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    # read in interpolated stellar ejected elements mass
                    MH_table = function_read_M_element("H", str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    MHe_table = function_read_M_element("He", str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    MC_table = function_read_M_element("C", str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    MN_table = function_read_M_element("N", str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    MO_table = function_read_M_element("O", str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    MMg_table = function_read_M_element("Mg", str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    MNe_table = function_read_M_element("Ne", str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    MSi_table = function_read_M_element("Si", str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    MS_table = function_read_M_element("S", str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    MCa_table = function_read_M_element("Ca", str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    MFe_table = function_read_M_element("Fe", str_evo_table, Z_select_in_table_2, Z_select_in_table_3)
                    M_element_table = [MH_table, MHe_table, MC_table, MN_table, MO_table, MMg_table, MNe_table, MSi_table, MS_table, MCa_table, MFe_table]

                    # check if the in put lifetime and final mass table used the same mass grid
                    if mass_1 != mass_12:
                        print('Error! Stellar lifetime and final mass input data do not match.\n'
                              'Check the table file: yield_tables/rearranged/setllar_final_mass_from_portinari98/portinari98_Z={}.txt\n'
                              'and table file: yield_tables/rearranged/setllar_lifetime_from_portinari98/portinari98_Z={}.txt'.format(
                                                                                               Z_select_in_table,
                                                                                               Z_select_in_table))
                    else:
                        mass_grid_table = mass
                        mass_grid_table2 = mass2
                    last_time_age = age_of_this_epoch
                    number_in_SNIa_boundary = mass_calibration_factor * quad(igimf_xi_function, 1.5, 8, limit=40)[0]  # see function_number_SNIa below
                    number_all = quad(igimf_xi_function, 0.08, steller_mass_upper_bound, limit=40)[0]  # see function_number_SNIa below
                    # number_low = quad(igimf_xi_function, 0.08, 2, limit=40)[0]  # see function_number_SNIa below
                    # number_up = quad(igimf_xi_function, 8, steller_mass_upper_bound, limit=40)[0]  # see function_number_SNIa below
                    # print("up", number_up/number_all)
                    SNIa_number_prob = number_in_SNIa_boundary**2 / number_all * 10**2 * 0.61
                    # number_in_SNIa_boundary = SNIa_number_prob
                    # SNIa_number_prob = number_in_SNIa_boundary / integrate_igimf_mass
                    # print("SNIa SNIa_number_prob:", SNIa_number_prob)
                    # print("total star number", number_all)
                    # print("low", number_low/number_all)

                    mass_boundary_at_end = fucntion_mass_boundary(age_of_this_epoch_at_end, mass_grid_table,lifetime_table)
                    all_sf_imf.append([igimf, mass_boundary_at_end, this_time])
                    all_sfr.append([S_F_R_of_this_epoch, age_of_this_epoch_at_end])
                    epoch_info.append(
                        [S_F_R_of_this_epoch, M_tot_of_this_epoch, igimf_of_this_epoch, integrate_igimf_mass,
                         mass_grid_table, lifetime_table, Mfinal_table, mass_grid_table2, Mmetal_table, M_element_table,
                         last_time_age, number_in_SNIa_boundary, metal_mass_fraction_in_gas, mass_calibration_factor])
                    metal_in_gas = metal_mass_fraction_in_gas
                else:  # if SFR == 0
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
                number_in_SNIa_boundary = epoch_info[epoch_index][11]
                metal_in_gas = epoch_info[epoch_index][12]
                mass_calibration_factor = epoch_info[epoch_index][13]
                def igimf_xi_function(mass):
                    return igimf_of_this_epoch.custom_imf(mass, this_time)
                def igimf_mass_function(mass):
                    return igimf_of_this_epoch.custom_imf(mass, this_time) * mass
                def igimf_luminous_function(mass):
                    if mass < 0.23**(1/1.7):
                        return igimf_of_this_epoch.custom_imf(mass, this_time) * 0.23 * mass**2.3
                    elif mass < 1.96:
                        return igimf_of_this_epoch.custom_imf(mass, this_time) * mass**4
                    elif mass < (32000/1.4)**(1/2.5):
                        return igimf_of_this_epoch.custom_imf(mass, this_time) * 1.4 * mass ** 3.5
                    else:
                        return igimf_of_this_epoch.custom_imf(mass, this_time) * 32000 * mass

            if S_F_R_of_this_epoch > 0:
                # get M_tot (total initial mass of all star ever formed)
                M_tot_of_this_time += M_tot_of_this_epoch
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

                    integrate_star_mass = quad(igimf_mass_function, 0.08, mass_boundary, limit=40)[0] # normalized mass
                    integrate_star_luminosity = quad(igimf_luminous_function, 0.08, mass_boundary, limit=40)[0] # normalized mass
                    current_stellar_mass_of_this_epoch = mass_calibration_factor * integrate_star_mass # real mass
                    current_stellar_luminosity_of_this_epoch = mass_calibration_factor * integrate_star_luminosity

                    # apprent metal mass (neglect stellar evolution, only account for the initial metal mass when SF):
                    stellar_metal_mass_of_this_epoch = current_stellar_mass_of_this_epoch * metal_in_gas[0]
                    stellar_H_mass_of_this_epoch = current_stellar_mass_of_this_epoch * metal_in_gas[1]
                    stellar_He_mass_of_this_epoch = current_stellar_mass_of_this_epoch * metal_in_gas[2]
                    # stellar_C_mass_of_this_epoch = current_stellar_mass_of_this_epoch * metal_in_gas[3]
                    # stellar_N_mass_of_this_epoch = current_stellar_mass_of_this_epoch * metal_in_gas[4]
                    stellar_O_mass_of_this_epoch = current_stellar_mass_of_this_epoch * metal_in_gas[5]
                    stellar_Mg_mass_of_this_epoch = current_stellar_mass_of_this_epoch * metal_in_gas[6]
                    # stellar_Ca_mass_of_this_epoch = current_stellar_mass_of_this_epoch * metal_in_gas[7]
                    stellar_Fe_mass_of_this_epoch = current_stellar_mass_of_this_epoch * metal_in_gas[8]

                    # apprent metal luminosity (neglect stellar evolution, only account for main-sequence stars):
                    stellar_metal_luminosity_of_this_epoch = current_stellar_luminosity_of_this_epoch * metal_in_gas[0]
                    stellar_H_luminosity_of_this_epoch = current_stellar_luminosity_of_this_epoch * metal_in_gas[1]
                    stellar_He_luminosity_of_this_epoch = current_stellar_luminosity_of_this_epoch * metal_in_gas[2]
                    # stellar_C_luminosity_of_this_epoch = current_stellar_luminosity_of_this_epoch * metal_in_gas[3]
                    # stellar_N_luminosity_of_this_epoch = current_stellar_luminosity_of_this_epoch * metal_in_gas[4]
                    stellar_O_luminosity_of_this_epoch = current_stellar_luminosity_of_this_epoch * metal_in_gas[5]
                    stellar_Mg_luminosity_of_this_epoch = current_stellar_luminosity_of_this_epoch * metal_in_gas[6]
                    # stellar_Ca_luminosity_of_this_epoch = current_stellar_luminosity_of_this_epoch * metal_in_gas[7]
                    stellar_Fe_luminosity_of_this_epoch = current_stellar_luminosity_of_this_epoch * metal_in_gas[8]
                    #
                    BH_mass_of_this_epoch = get_BH_mass(mass_boundary, 1, 1, mass_calibration_factor, steller_mass_upper_bound)
                    NS_mass_of_this_epoch = get_NS_mass(mass_boundary, 1, 1, mass_calibration_factor)
                    WD_mass_of_this_epoch = get_WD_mass(mass_boundary, 1, 1, mass_calibration_factor)
                    remnant_mass_of_this_epoch = WD_mass_of_this_epoch + NS_mass_of_this_epoch + BH_mass_of_this_epoch
                    ejected_gas_mass_of_this_epoch = M_tot_of_this_epoch - current_stellar_mass_of_this_epoch - remnant_mass_of_this_epoch
                    #
                    # # consider direct black hole as in Heger et al. (2003) (maybe not self-consistant with the stellar evolution table)
                    # if mass_boundary > 100:
                    #     SNII_number_of_this_epoch_1 = quad(igimf_mass_function, mass_boundary, steller_mass_upper_bound, limit=40)[0]
                    #     SNII_number_of_this_epoch_2 = 0
                    # elif mass_boundary > 40:
                    #     SNII_number_of_this_epoch_1 = quad(igimf_mass_function, 100, steller_mass_upper_bound, limit=40)[0]
                    #     SNII_number_of_this_epoch_2 = 0
                    # elif mass_boundary > 8:
                    #     SNII_number_of_this_epoch_1 = quad(igimf_mass_function, 100, steller_mass_upper_bound, limit=40)[0]
                    #     SNII_number_of_this_epoch_2 = quad(igimf_mass_function, mass_boundary, 40, limit=40)[0]
                    # else:
                    #     SNII_number_of_this_epoch_1 = quad(igimf_mass_function, 100, steller_mass_upper_bound, limit=40)[0]
                    #     SNII_number_of_this_epoch_2 = quad(igimf_mass_function, 8, 40, limit=40)[0]
                    # SNII_number_of_this_epoch = (SNII_number_of_this_epoch_1 + SNII_number_of_this_epoch_2) * mass_calibration_factor
                    if mass_boundary > 8:
                        SNII_number_of_this_epoch = quad(igimf_mass_function, mass_boundary, steller_mass_upper_bound, limit=40)[0]
                        SNII_ejected_mass_of_this_epoch = quad(igimf_mass_function, mass_boundary, steller_mass_upper_bound, limit=40)[0]
                    else:
                        SNII_number_of_this_epoch = quad(igimf_mass_function, 8, steller_mass_upper_bound, limit=40)[0]
                    SNII_number_of_this_epoch = SNII_number_of_this_epoch * mass_calibration_factor
                    SNII_energy_release_per_event = 10**51
                    SNII_number += SNII_number_of_this_epoch
                    SNII_energy_release += SNII_energy_release_per_event * SNII_number_of_this_epoch
                    # ejected_ :
                    metal_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound, 2, 2, mass_calibration_factor)
                    H_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound, 2, "H", mass_calibration_factor)
                    He_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound, 2, "He", mass_calibration_factor)
                    # C_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound, 2, "C", mass_calibration_factor)
                    # N_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound, 2, "N", mass_calibration_factor)
                    O_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound, 2, "O", mass_calibration_factor)
                    Mg_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound, 2, "Mg", mass_calibration_factor)
                    # Ca_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound, 2, "Ca", mass_calibration_factor)
                    Fe_mass_of_this_epoch = function_get_target_mass_in_range(mass_boundary, steller_mass_upper_bound, 2, "Fe", mass_calibration_factor)
                    Fe_production_SNII += Fe_mass_of_this_epoch
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
                    current_stellar_mass_of_this_epoch = 0
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
                    Fe_mass_of_this_epoch = 0
                # if consider SNIa
                if SNIa_ON == True:
                    # read in SNIa yield table
                    # TNH93: see Gibson, B. K., Loewenstein, M., & Mushotzky, R. F. 1997, MNRAS, 290, 623
                    # based on the work of Thielemann et al. (1993)
                    Fe_mass_eject = 0.744  # Nomoto 1984 0.613  # TNH93 0.744 i99CDD1 0.56, i99CDD2 0.76, i99W7 0.63, ivo12/13 0.62-0.67, t03 0.74, t86 0.63
                    Si_mass_eject = 0.158
                    O_mass_eject = 0.148 # Nomoto 1984 0.140 TNH93 0.148 i99CDD1 0.09, i99CDD2 0.06, i99W7 0.14, ivo12/13 0.09-0.1, t03 0.14, t86 0.13
                    S_mass_eject = 0.086
                    Mg_mass_eject = 0.009  # Nomoto 1984 0.023 TNH93 0.009 i99CDD1 0.0077, i99CDD2 0.0042, i99W7 0.0085, ivo12/13 0.015-0.029, t03 0.013, t86 0.016
                    Ne_mass_eject = 0.005
                    Chandrasekhar_mass = 1.44
                    pre_SNIa_NS_mass = 1
                    SNIa_energy_release_per_event = 10**51 # in the unit of 10^51 erg
                    # integrate SNIa number from last_delay_time to this_delay_time contributed by this SF epoch
                    SNIa_number_from_this_epoch_till_this_time = function_number_SNIa(0, age_of_this_epoch,
                                                                                      number_in_SNIa_boundary,
                                                                                      S_F_R_of_this_epoch)

                    # the following should result in 0.0022+-50% for a SSP,
                    # but now calibrate to a different value to fit with galaxy [Fe/H] observation
                    if age_of_this_epoch == 10*10**9-1*10**7:
                        # print(function_number_SNIa(0, 10 * 10 ** 9, 1, 0))
                        # print("SN number per star in range:", SNIa_number_from_this_epoch_till_this_time/number_in_SNIa_boundary)
                        print("SNIa number within 10Gyr per solar mass of star:", SNIa_number_from_this_epoch_till_this_time/M_tot_of_this_epoch)
                    # update the element masses
                    ejected_gas_mass_of_this_epoch += pre_SNIa_NS_mass * SNIa_number_from_this_epoch_till_this_time
                    metal_mass_of_this_epoch += (Chandrasekhar_mass - (Chandrasekhar_mass - pre_SNIa_NS_mass) *
                                                 Z_gas_this_time_step) * SNIa_number_from_this_epoch_till_this_time
                    O_mass_of_SNIa = O_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    Mg_mass_of_SNIa = Mg_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    Fe_mass_of_SNIa = (Fe_mass_eject
                                              #- (Chandrasekhar_mass - pre_SNIa_NS_mass) * Fe_H_mass_ratio_at_last_time * 0.7057 # this term is small and can be neglected
                                              ) * SNIa_number_from_this_epoch_till_this_time
                    O_mass_of_this_epoch += O_mass_of_SNIa
                    Mg_mass_of_this_epoch += Mg_mass_of_SNIa
                    Fe_mass_of_this_epoch += Fe_mass_of_SNIa
                    Fe_production_SNIa += Fe_mass_of_this_epoch
                    Mg_production_SNIa += Mg_mass_of_this_epoch
                    O_production_SNIa += O_mass_of_this_epoch
                    # if age_of_this_epoch == 1 * 10 ** 9:
                    #     print("Fe_production_SNIa", Fe_mass_of_SNIa)
                    #     print("O_production_SNIa", O_mass_of_SNIa)
                    #     print("Mg_production_SNIa", Mg_mass_of_SNIa)

                    remnant_mass_of_this_epoch -= pre_SNIa_NS_mass * SNIa_number_from_this_epoch_till_this_time
                    WD_mass_of_this_epoch -= pre_SNIa_NS_mass * SNIa_number_from_this_epoch_till_this_time
                    SNIa_number += SNIa_number_from_this_epoch_till_this_time
                    SNIa_energy_release += SNIa_energy_release_per_event * SNIa_number_from_this_epoch_till_this_time
                #
                stellar_mass_at_this_time += current_stellar_mass_of_this_epoch
                stellar_metal_mass_at_this_time += stellar_metal_mass_of_this_epoch
                stellar_H_mass_at_this_time += stellar_H_mass_of_this_epoch
                stellar_He_mass_at_this_time += stellar_He_mass_of_this_epoch
                stellar_O_mass_at_this_time += stellar_O_mass_of_this_epoch
                stellar_Mg_mass_at_this_time += stellar_Mg_mass_of_this_epoch
                stellar_Fe_mass_at_this_time += stellar_Fe_mass_of_this_epoch
                #
                stellar_luminosity_at_this_time += current_stellar_luminosity_of_this_epoch
                stellar_metal_luminosity_at_this_time += stellar_metal_luminosity_of_this_epoch
                stellar_H_luminosity_at_this_time += stellar_H_luminosity_of_this_epoch
                stellar_He_luminosity_at_this_time += stellar_He_luminosity_of_this_epoch
                stellar_O_luminosity_at_this_time += stellar_O_luminosity_of_this_epoch
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
                ejected_Mg_mass_till_this_time += Mg_mass_of_this_epoch
                ejected_Fe_mass_till_this_time += Fe_mass_of_this_epoch
            # Goes to the next SF epoch until all SF event before this time step is accounted:
            (epoch_index) = (epoch_index + 1)
        # output of this time step
        total_energy_release = SNIa_energy_release + SNII_energy_release

        ### yeilds at this time step:
        ejected_gas_mass_at_this_time = ejected_gas_mass_till_this_time - ejected_gas_mass_till_last_time
        ejected_metal_mass_at_this_time = ejected_metal_mass_till_this_time - ejected_metal_mass_till_last_time
        ejected_H_mass_at_this_time = ejected_H_mass_till_this_time - ejected_H_mass_till_last_time
        ejected_He_mass_at_this_time = ejected_He_mass_till_this_time - ejected_He_mass_till_last_time
        ejected_C_mass_at_this_time = ejected_C_mass_till_this_time - ejected_C_mass_till_last_time
        ejected_N_mass_at_this_time = ejected_N_mass_till_this_time - ejected_N_mass_till_last_time
        ejected_O_mass_at_this_time = ejected_O_mass_till_this_time - ejected_O_mass_till_last_time
        ejected_Mg_mass_at_this_time = ejected_Mg_mass_till_this_time - ejected_Mg_mass_till_last_time
        ejected_Ca_mass_at_this_time = ejected_Ca_mass_till_this_time - ejected_Ca_mass_till_last_time
        ejected_Fe_mass_at_this_time = ejected_Fe_mass_till_this_time - ejected_Fe_mass_till_last_time
        ejected_gas_Mg_over_Fe_till_this_time = function_element_abundunce("Mg", "Fe", ejected_Mg_mass_till_this_time, ejected_Fe_mass_till_this_time)
        ejected_gas_Mg_over_Fe_at_this_time = function_element_abundunce("Mg", "Fe", ejected_Mg_mass_at_this_time, ejected_Fe_mass_at_this_time)
        M_tot_at_this_time = M_tot_of_this_time - M_tot_of_last_time # new SF mass added at this time step
        #
        galaxy_mass_without_gas_at_this_time = stellar_mass_at_this_time + remnant_mass_at_this_time
        if galaxy_mass_without_gas_at_this_time == 0 or ejected_gas_mass_at_this_time == 0:
            expansion_factor_instantaneous = 1
            expansion_factor_adiabat = 1
        elif galaxy_mass_without_gas_at_this_time < ejected_gas_mass_at_this_time:
            print('Warning: galaxy_mass < ejected_gas_mass. '
                  'This is due to too large a timestep. '
                  'It is easy to aviod this issue by applying the "high_time_resolution=True". '
                  'But the simulation will take much longer time.')
            expansion_factor_instantaneous = 10
            expansion_factor_adiabat = (galaxy_mass_without_gas_at_this_time + ejected_gas_mass_at_this_time) / galaxy_mass_without_gas_at_this_time
        else:
            expansion_factor_instantaneous = galaxy_mass_without_gas_at_this_time / (galaxy_mass_without_gas_at_this_time - ejected_gas_mass_at_this_time)
            expansion_factor_adiabat = (galaxy_mass_without_gas_at_this_time + ejected_gas_mass_at_this_time) / galaxy_mass_without_gas_at_this_time

        expansion_factor = 10 ** (
            (math.log(expansion_factor_instantaneous, 10) + math.log(expansion_factor_adiabat, 10)) / 2)

        ### Element abundances in the gas phase (in solar unit):
        if outflow is True:
            lockup_and_outflow_mass = M_tot_at_this_time * 2  # lockup gas mass in BDs is about 4% thus neglected while the uniform outflow is often assumed to be the same value as the formed stellar mass.
        else:
            lockup_and_outflow_mass = M_tot_at_this_time
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
        total_Fe_mass_at_this_time = total_Fe_mass_at_last_time - lockup_and_outflow_mass * (
            total_Fe_mass_at_last_time / total_gas_mass_at_last_time) + ejected_Fe_mass_at_this_time
        if total_Fe_mass_at_this_time < 0.0001:
            total_Fe_mass_at_this_time = 0.0001

        # # calculate the gravitational binding engergy:
        #
        # # galaxy_mass_without_gas_at_this_time, original_gas_mass, total_gas_mass_at_this_time, ejected_gas_mass_at_this_time
        # # gas_mass = max(ejected_gas_mass_at_this_time, 1)
        #
        # # galaxy mass--radii relation adopted from Dabringhausen 2008 eq.4
        # Dabringhausen_2008_a = 2.95
        # Dabringhausen_2008_b = 0.596
        # initial_expansion_factor = 10000  # need change for every simulation. use the expansion_factor at final time
        # # initial_expansion_factor = expansion_factor_list[-1]
        # if expansion_factor_list == []:
        #     current_expansion_factor = initial_expansion_factor
        # else:
        #     current_expansion_factor = initial_expansion_factor - expansion_factor_list[-1]
        #
        # log_binding_energy = round(
        #     math.log(3 / 5 * gravitational_constant * 1.989**2 / 3.086, 10) + 40 + (2 - Dabringhausen_2008_b) *
        #     math.log(original_gas_mass, 10) - math.log(Dabringhausen_2008_a, 10) +
        #     6 * Dabringhausen_2008_b + math.log(current_expansion_factor, 10), 3)
        # # 40 = 30 (solar mass) * 2 - 11 (Gravitational constant) - 16 (pc to meter) + 7 (J to erg)
        # binding_energy = 10 ** log_binding_energy  # [erg]

        # calculate the kinetic energy of the gas if they have a uniform temperature of 2 keV:
        all_masses = total_H_mass_at_this_time + total_He_mass_at_this_time + total_metal_mass_at_this_time
        X_for_H = total_H_mass_at_this_time / all_masses
        Y_for_He = total_He_mass_at_this_time / all_masses
        Z_for_metal = total_metal_mass_at_this_time / all_masses
        mean_molecular_weight = 1 / (2*X_for_H + 3/4*Y_for_He + Z_for_metal/2) * \
                                element_weight_table.function_element_weight("H") / 6.022140857 / 1.9891
        # / 10**23 / 10**33 (i.e., 10**56) mean_molecular_weight in solar mass unit.
        log_mean_molecular_weight = math.log(mean_molecular_weight, 10) - 56  # log [M_sun]
        log_total_number_of_molecule = math.log(total_gas_mass_at_this_time, 10) - log_mean_molecular_weight  # log [Number]
        # 1 [keV] = 1.60217662 * 10**(-9) [erg]
        log_energy_per_molecule = math.log(2 * 1.60217662, 10)-9  # [erg]
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
        Z_over_X = math.log(total_metal_mass_at_this_time/total_H_mass_at_this_time, 10) - math.log(Z_solar/0.7381, 10)
        # Fe_H_mass_ratio_at_this_time = total_Fe_mass_at_this_time / total_H_mass_at_this_time
        gas_Y_at_this_time = total_He_mass_at_this_time / total_gas_mass_at_this_time
        O_over_H_number_ratio = function_element_abundunce("O", "H", total_O_mass_at_this_time, total_H_mass_at_this_time)
        Mg_over_H_number_ratio = function_element_abundunce("Mg", "H", total_Mg_mass_at_this_time, total_H_mass_at_this_time)
        Fe_over_H_number_ratio = function_element_abundunce("Fe", "H", total_Fe_mass_at_this_time, total_H_mass_at_this_time)
        C_over_Fe_number_ratio = function_element_abundunce("C", "Fe", total_C_mass_at_this_time, total_Fe_mass_at_this_time)
        N_over_Fe_number_ratio = function_element_abundunce("N", "Fe", total_N_mass_at_this_time, total_Fe_mass_at_this_time)
        O_over_Fe_number_ratio = function_element_abundunce("O", "Fe", total_O_mass_at_this_time, total_Fe_mass_at_this_time)
        Mg_over_Fe_number_ratio = function_element_abundunce("Mg", "Fe", total_Mg_mass_at_this_time, total_Fe_mass_at_this_time)
        Ca_over_Fe_number_ratio = function_element_abundunce("Ca", "Fe", total_Ca_mass_at_this_time, total_Fe_mass_at_this_time)

        ### Element abundances in of stars (consider only the metal of stellar surface, i.e., neglect stellar evolution
        # This raises errors from very low mass stars which are fully convective but may not be observationally important):
        ##### mass weighted abundances
        # (total metal in stars / total H in stars):
        mass_weighted_stellar_Y = stellar_He_mass_at_this_time / stellar_mass_at_this_time
        mass_weighted_stellar_O_over_H = function_element_abundunce("O", "H", stellar_O_mass_at_this_time, stellar_H_mass_at_this_time)
        mass_weighted_stellar_Mg_over_H = function_element_abundunce("Mg", "H", stellar_Mg_mass_at_this_time, stellar_H_mass_at_this_time)
        mass_weighted_stellar_Fe_over_H = function_element_abundunce("Fe", "H", stellar_Fe_mass_at_this_time, stellar_H_mass_at_this_time)
        # mass_weighted_stellar_C_over_Fe = function_element_abundunce("C", "Fe", stellar_C_mass_at_this_time, stellar_Fe_mass_at_this_time)
        # mass_weighted_stellar_N_over_Fe = function_element_abundunce("N", "Fe", stellar_N_mass_at_this_time, stellar_Fe_mass_at_this_time)
        mass_weighted_stellar_O_over_Fe = function_element_abundunce("O", "Fe", stellar_O_mass_at_this_time, stellar_Fe_mass_at_this_time)
        mass_weighted_stellar_Mg_over_Fe = function_element_abundunce("Mg", "Fe", stellar_Mg_mass_at_this_time, stellar_Fe_mass_at_this_time)
        # mass_weighted_stellar_Ca_over_Fe = function_element_abundunce("Ca", "Fe", stellar_Ca_mass_at_this_time, stellar_Fe_mass_at_this_time)

        ##### luminosity weighted abundances
        # (total metal in stars / total H in stars):
        luminosity_weighted_stellar_Y = stellar_He_luminosity_at_this_time / stellar_luminosity_at_this_time
        luminosity_weighted_stellar_O_over_H = function_element_abundunce("O", "H", stellar_O_luminosity_at_this_time, stellar_H_luminosity_at_this_time)
        luminosity_weighted_stellar_Mg_over_H = function_element_abundunce("Mg", "H", stellar_Mg_luminosity_at_this_time, stellar_H_luminosity_at_this_time)
        luminosity_weighted_stellar_Fe_over_H = function_element_abundunce("Fe", "H", stellar_Fe_luminosity_at_this_time, stellar_H_luminosity_at_this_time)
        # luminosity_weighted_stellar_C_over_Fe = function_element_abundunce("C", "Fe", stellar_C_luminosity_at_this_time, stellar_Fe_luminosity_at_this_time)
        # luminosity_weighted_stellar_N_over_Fe = function_element_abundunce("N", "Fe", stellar_N_luminosity_at_this_time, stellar_Fe_luminosity_at_this_time)
        luminosity_weighted_stellar_O_over_Fe = function_element_abundunce("O", "Fe", stellar_O_luminosity_at_this_time, stellar_Fe_luminosity_at_this_time)
        luminosity_weighted_stellar_Mg_over_Fe = function_element_abundunce("Mg", "Fe", stellar_Mg_luminosity_at_this_time, stellar_Fe_luminosity_at_this_time)
        # luminosity_weighted_stellar_Ca_over_Fe = function_element_abundunce("Ca", "Fe", stellar_Ca_luminosity_at_this_time, stellar_Fe_luminosity_at_this_time)


        if stellar_H_mass_at_this_time == 0:
            mass_weighted_stellar_Z_over_X = math.log(Z_0 / Z_solar, 10)
            luminosity_weighted_stellar_Z_over_X = mass_weighted_stellar_Z_over_X
        else:
            mass_weighted_stellar_Z_over_X = math.log(stellar_metal_mass_at_this_time / stellar_H_mass_at_this_time, 10) \
                                         - math.log(Z_solar / 0.7381, 10)
            luminosity_weighted_stellar_Z_over_X = math.log(stellar_metal_luminosity_at_this_time / stellar_H_luminosity_at_this_time, 10) \
                                         - math.log(Z_solar / 0.7381, 10)

        # the so called [Z/H] as determined by observation with equation: [Z/H] = [Fe/H] + A[Mg/Fe] where A=0.94 (Thomas 2003)
        if mass_weighted_stellar_Fe_over_H is None or mass_weighted_stellar_Mg_over_Fe is None:
            mass_weighted_stellar_Z_over_H = None
        else:
            mass_weighted_stellar_Z_over_H = mass_weighted_stellar_Fe_over_H + 0.94 * mass_weighted_stellar_Mg_over_Fe
            luminosity_weighted_stellar_Z_over_H = luminosity_weighted_stellar_Fe_over_H + 0.94 * luminosity_weighted_stellar_Mg_over_Fe

        ##########################################
        ##### luminosity weighted abundances #####
        # (......):
        ##########################################

        if BH_mass_till_this_time == 0:
            BH_mass_list += [0.1]
        else:
            BH_mass_list += [BH_mass_till_this_time]

        if NS_mass_till_this_time == 0:
            NS_mass_list += [0.1]
        else:
            NS_mass_list += [NS_mass_till_this_time]

        if WD_mass_till_this_time == 0:
            WD_mass_list += [0.1]
        else:
            WD_mass_list += [WD_mass_till_this_time]

        if Z_over_X == 0:
            print("Warning: Z_over_X == 0")
            gas_Z_over_X_list += [math.log(Z_0 / Z_solar, 10)]
        else:
            gas_Z_over_X_list += [Z_over_X]

        Y_list += [gas_Y_at_this_time]
        O_over_H_list += [O_over_H_number_ratio]
        Mg_over_H_list += [Mg_over_H_number_ratio]
        Fe_over_H_list += [Fe_over_H_number_ratio]
        C_over_Fe_list += [C_over_Fe_number_ratio]
        N_over_Fe_list += [N_over_Fe_number_ratio]
        O_over_Fe_list += [O_over_Fe_number_ratio]
        Mg_over_Fe_list += [Mg_over_Fe_number_ratio]
        Ca_over_Fe_list += [Ca_over_Fe_number_ratio]


        stellar_O_over_H_list += [mass_weighted_stellar_O_over_H]
        stellar_Mg_over_H_list += [mass_weighted_stellar_Mg_over_H]
        stellar_Y_list += [mass_weighted_stellar_Y]
        stellar_Fe_over_H_list += [mass_weighted_stellar_Fe_over_H]
        # stellar_C_over_Fe_list += [mass_weighted_stellar_C_over_Fe]
        # stellar_N_over_Fe_list += [mass_weighted_stellar_N_over_Fe]
        stellar_O_over_Fe_list += [mass_weighted_stellar_O_over_Fe]
        stellar_Mg_over_Fe_list += [mass_weighted_stellar_Mg_over_Fe]
        # stellar_Ca_over_Fe_list += [mass_weighted_stellar_Ca_over_Fe]
        stellar_Z_over_X_list += [mass_weighted_stellar_Z_over_X]
        stellar_Z_over_H_list += [mass_weighted_stellar_Z_over_H]

        stellar_O_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_O_over_H]
        stellar_Mg_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_Mg_over_H]
        stellar_Y_list_luminosity_weighted += [luminosity_weighted_stellar_Y]
        stellar_Fe_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_Fe_over_H]
        # stellar_C_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_C_over_Fe]
        # stellar_N_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_N_over_Fe]
        stellar_O_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_O_over_Fe]
        stellar_Mg_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_Mg_over_Fe]
        # stellar_Ca_over_Fe_list_luminosity_weighted += [luminosity_weighted_stellar_Ca_over_Fe]
        stellar_Z_over_X_list_luminosity_weighted += [luminosity_weighted_stellar_Z_over_X]
        stellar_Z_over_H_list_luminosity_weighted += [luminosity_weighted_stellar_Z_over_H]

        if remnant_mass_at_this_time == 0:
            remnant_mass_list += [0.1]
        else:
            remnant_mass_list += [remnant_mass_at_this_time]

        if total_gas_mass_at_this_time == 0:
            total_gas_mass_list += [0.1]
        else:
            total_gas_mass_list += [total_gas_mass_at_this_time]

        if ejected_gas_mass_till_this_time == 0 or ejected_gas_mass_till_this_time < 0:
            ejected_gas_mass_list += [0.1]
        else:
            ejected_gas_mass_list += [ejected_gas_mass_till_this_time]

        ejected_gas_Mg_over_Fe_list += [ejected_gas_Mg_over_Fe_till_this_time]
        instant_ejected_gas_Mg_over_Fe_list += [ejected_gas_Mg_over_Fe_at_this_time]

        if M_tot_of_this_time > 0:
            ejected_metal_mass_list += [ejected_metal_mass_till_this_time/M_tot_of_this_time]
        else:
            ejected_metal_mass_list += [0]

        if expansion_factor_instantaneous_list == []:
            expansion_factor_instantaneous_list += [1]
        else:
            expansion_factor_instantaneous_list += [expansion_factor_instantaneous * expansion_factor_instantaneous_list[-1]]

        if expansion_factor_adiabat_list == []:
            expansion_factor_adiabat_list += [1]
        else:
            expansion_factor_adiabat_list += [expansion_factor_adiabat * expansion_factor_adiabat_list[-1]]

        if expansion_factor_list == []:
            expansion_factor_list += [1]
        else:
            expansion_factor_list += [expansion_factor * expansion_factor_list[-1]]

        if stellar_mass_at_this_time == 0:
            stellar_mass_list += [0.1]
        else:
            stellar_mass_list += [stellar_mass_at_this_time]

        SNIa_energy_release_list += [SNIa_energy_release]
        SNIa_number_list += [SNIa_number]

        if len(SNIa_number_per_century) == 0:
            SNIa_number_per_century += [SNIa_number_list[0]]
        else:
            SNIa_number_per_century += [(SNIa_number_list[-1] - SNIa_number_list[-2])/(time_axis[time_step]-time_axis[time_step-1])*100]

        SNII_energy_release_list += [SNII_energy_release]
        SNII_number_list += [SNII_number]

        if len(SNII_number_per_century) == 0:
            SNII_number_per_century += [SNII_number_list[0]]
        else:
            SNII_number_per_century += [(SNII_number_list[-1] - SNII_number_list[-2])/(time_axis[time_step]-time_axis[time_step-1])*100]

        if total_energy_release == 0:
            total_energy_release_list += [0]
        else:
            total_energy_release_list += [total_energy_release]#[math.log((total_energy_release), 10)]

        # if binding_energy == 0:
        #     binding_energy_list += [0]
        # else:
        #     binding_energy_list += [binding_energy]#[math.log((binding_energy), 10)]

        if total_gas_kinetic_energy_list == 0:
            total_gas_kinetic_energy_list += [0]
        else:
            total_gas_kinetic_energy_list += [total_gas_kinetic_energy]#[math.log((total_gas_kinetic_energy), 10)]

        SN_number_per_century += [SNIa_number_per_century[-1] + SNII_number_per_century[-1]]

        # go to next time step
        if time_step/50 > gc_collect_check:
            gc_collect_check += 1
            print("gc_collect:", gc_collect_check)
            gc.collect()
        (time_step)=(time_step + 1)

    ###################
    ### output data ###
    ###################

    # Remnant_Star_ratio = [0]*len(stellar_mass_list)
    # for i in range(len(remnant_mass_list)):
    #     Remnant_Star_ratio[i] = remnant_mass_list[i]/stellar_mass_list[i]
    # import csv
    # with open('GalEvo_time.txt', 'w') as f:
    #     writer = csv.writer(f, delimiter=' ')
    #     f.write("# galaxy_evol.py output file.\n# time\n")
    #     writer.writerows(
    #         zip(time_axis))
    # with open('GalEvo_ratio.txt', 'w') as f:
    #     writer = csv.writer(f, delimiter=' ')
    #     f.write("# galaxy_evol.py output file.\n# Remnant_Star_ratio\n")
    #     writer.writerows(
    #         zip(Remnant_Star_ratio))


    ###################
    ### output plot ###
    ###################

    text_output(imf, SFE, round(math.log(max(SFH_input), 10), 1), SFEN, original_gas_mass, Z_0, Z_solar)

    print(" - Run time: %s -" % round((time.time() - start_time), 2))

    # if output plot applies
    plot_output(plot_show, plot_save, imf, igimf)

    ###################
    ###     end     ###
    ###################
    return

# def function_update_element_gas_infall():
#     return


# # calculate the diet_Salpeter_number_to_mass_ratio:
# Bell & de Jong (2001). Salpeter IMF x = 1.35 with a flat x = 0 slope below 0.35
def function_xi_diet_Salpeter_IMF(mass):
    # integrate this function's output xi result in the number of stars in mass limits.
    xi = diet_Salpeter_IMF.custom_imf(mass, 0)
    return xi

def function_mass_diet_Salpeter_IMF(mass):
    # integrate this function's output m result in the total stellar mass for stars in mass limits.
    m = mass * diet_Salpeter_IMF.custom_imf(mass, 0)
    return m

integrate_all_for_function_mass_SNIa = quad(function_mass_diet_Salpeter_IMF, 0.1, 100, limit=40)[0]
integrate_28_for_function_number_SNIa = quad(function_xi_diet_Salpeter_IMF, 1.5, 8, limit=40)[0]
diet_Salpeter_number_to_mass_ratio = integrate_all_for_function_mass_SNIa / integrate_28_for_function_number_SNIa

def function_number_SNIa(last_delay_time, this_delay_time, stellar_number_in_SNIa_boundary, S_F_R_of_this_epoch):
    # It is commonly assumed that the maximum stellar mass able to produce a degenerate C–O white dwarf is 8 M⊙,
    # The minimum possible binary mass is assumed to be 3 M⊙ in order to ensure that the
    # smallest possible white dwarf can accrete enough mass from the secondary star to reach the Chandrasekhar mass.
    # see Greggio, L., & Renzini, A. 1983, A & A, 118, 217
    # Thus we should normalize the DTD according to the number (but currently, mass) of stars between 1.5 and 8 solar mass
    # normalized with a SNIa assuming fixed diet-Salpeter IMF (Bell et al. 149:289–312, 2003)
    # See Dan Maoz and Filippo Mannucci 2012 review
    global diet_Salpeter_number_to_mass_ratio
    SNIa_normalization_parameter = funtion_SNIa_DTD_normalization_parameter(S_F_R_of_this_epoch)
    # integrate SNIa number from last_delay_time to this_delay_time
    diet_Salpeter_SNIa_number_per_solar_mass =quad(function_SNIa_DTD, last_delay_time, this_delay_time, limit=40)[0]
    # calculate actual SNIa event number
    SNIa_number = stellar_number_in_SNIa_boundary *  SNIa_normalization_parameter * diet_Salpeter_SNIa_number_per_solar_mass * diet_Salpeter_number_to_mass_ratio
    # if this_delay_time == 1 * 10 ** 9:
    #     print("stellar_number_in_SNIa_boundary ===", stellar_number_in_SNIa_boundary)
    #     print("diet_Salpeter_SNIa_number_per_solar_mass", diet_Salpeter_SNIa_number_per_solar_mass)
    #     print("SNIa_number ===", SNIa_number)

    # SNIa_number = S_F_R_of_this_epoch * 10**7 * 0.0023 * diet_Salpeter_SNIa_number_per_solar_mass / quad(function_SNIa_DTD, 0, 10**10, limit=40)[0]
    return SNIa_number

def funtion_SNIa_DTD_normalization_parameter(SFR):
    SFR = 0.0001  # uncommen this line to enable the renormalizaion of SNIa number as a function of SFR
    x = 2
    SFRmatter = SFR + x
    logSFR = math.log(SFRmatter, 10)
    gamma = 2
    SNIa_normalization_parameter = 1 * (logSFR**gamma + 4)
    # This is a toy model relation. Roughly keep the SNIa-rate/stellar-mass-formed unchanged for different SFRs (IMFs).
    # When SFR is high, top-heavy IMF reduce the number possible SNIa progenitor within mass range 1.5-8 solar mass.
    # the dense stellar region pump the SNIa rate as the stars have a larger chance to meet with another star.
    # The lower SFR do not further reduce the SNIa number as the low SFR happens after the star burst epoch
    # thus the newly formed star can still meet with stars formed at ealier epochs regredless of its current SFR.
    return SNIa_normalization_parameter

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
# fig = plt.figure(200, figsize=(4, 3.5))
# fig.add_subplot(1, 1, 1)
# plt.plot(xxx, yyy)
# plt.xlabel(r'log$_{10}$(SFR [M$_\odot$/yr])')
# plt.ylabel(r'SNIa Number renormalization')
# plt.legend()
# plt.tight_layout()
# plt.show()




def function_SNIa_DTD(delay_time):
    # The delay time distribution (DTD) in the unit of per year per total stellar mass [solar]
    # DTD for SNIa is adopted from Maoz & Mannucci 2012, 29, 447–465, their equation 13
    # with a consistent assumed IMF – the Bell et al. 2003 diet-Salpeter IMF
    if delay_time < 4 * 10 **7: # [yr] #  2.3 * 10 ** 7 for a burst of star formation from Greggio 1983
        number = 0
    else:
        number = 10 ** (-4) * delay_time**(-1)  # Normalized such the DTD integral over 10Gyr for IGIMF is N_SN/M_sun = 2.23* 10^-3 (M_sun^-1)
        # number = 2 * 10 ** (-1) * delay_time**(-1.3)  #
        # number = 9.3200246746 * 10 ** (-4) * delay_time**(-1)  # Normalized such the DTD integral over 10Gyr for a diet-Salpeter IMF is N_SN/M_sun = 3* 10^-3 (M_sun^-1).
        # This value changes with igimf where top-heavy and bottom-heavy IGIMF will have lower number of SNIa
        # as the number of stars within the mass range 3.0001 to 8 solar mass is smaller.
        # The observational uncertainty being +-25%. See Maoz & Mannucci 2012 their Table 1
    return number

def function_read_lifetime(str_evo_table, Z_select_in_table):
    file_lifetime = open(
        'yield_tables/rearranged/setllar_lifetime_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table),
        'r')
    data = file_lifetime.readlines()
    metallicity = data[1]
    mass_1 = data[3]
    lifetime_ = data[5]
    file_lifetime.close()
    mass = [float(x) for x in mass_1.split()]
    lifetime_table = [float(x) for x in lifetime_.split()]
    return (mass_1, mass, lifetime_table)

def function_read_Mfinal(str_evo_table, Z_select_in_table):
    file_final_mass = open(
        "yield_tables/rearranged/setllar_final_mass_from_portinari98/portinari98_Z={}.txt".format(Z_select_in_table),
        'r')
    data = file_final_mass.readlines()
    metallicity2 = data[1]
    mass_2 = data[3]
    Mfinal_ = data[5]
    file_final_mass.close()
    Mfinal_table = [float(x) for x in Mfinal_.split()]
    return (mass_2, Mfinal_table)

def lindexsplit(List,*lindex):
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
    lastindexval = index[(len(index)-1)]
    finalcounttrigger = (totalitems-(lastindexval+1))
    for item in List:
        itemcounter += 1
        indexofitem = itemcounter - 1
        nextbreakindex = index[breakcounter]
        #Less than the last cut
        if breakcounter <= numberofbreaks:
            if indexofitem < nextbreakindex:
                templist1.append(item)
            elif breakcounter < (numberofbreaks - 1):
                templist1.append(item)
                templist2.append(templist1)
                templist1 = []
                breakcounter +=1
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

def function_read_Mmetal(str_evo_table, Z_select_in_table_2, Z_select_in_table_3):
    global mm, zz
    if str_evo_table == "portinari98":
        file_Metal_eject = open(
            'yield_tables/rearranged/setllar_Metal_eject_mass_from_{}/{}_Z={}.txt'.format(str_evo_table, str_evo_table, Z_select_in_table_2),
            'r')
        data = file_Metal_eject.readlines()
        metallicity = data[1]
        mass_2 = data[3]
        Metal_eject_ = data[5]
        file_Metal_eject.close()
        mass = [float(x) for x in mass_2.split()]
        Metal_eject_table = [float(x) for x in Metal_eject_.split()]
    elif str_evo_table == "WW95":
        file_Metal_eject = open(
            'yield_tables/rearranged/setllar_Metal_eject_mass_from_{}/{}_Z={}.txt'.format(str_evo_table, str_evo_table, Z_select_in_table_2),
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
            'yield_tables/rearranged/setllar_Metal_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
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
    return (mass_2, mass, Metal_eject_table)

def function_read_M_element(element, str_evo_table, Z_select_in_table_2, Z_select_in_table_3):
    if str_evo_table == "portinari98":
        if element == "H":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_H_eject_mass_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "He":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_He_eject_mass_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "C":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_C_eject_mass_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "N":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_N_eject_mass_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "O":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_O_eject_mass_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Mg":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Mg_eject_mass_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Ne":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Ne_eject_mass_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Si":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Si_eject_mass_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "S":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_S_eject_mass_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Ca":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Ca_eject_mass_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Fe":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Fe_eject_mass_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        else:
            file_M_eject = 0
            print("Error: element parameter for function_read_M_element do not exsit.")
        data = file_M_eject.readlines()
        M_eject_ = data[5]
        file_M_eject.close()
        M_eject_table = [float(x) for x in M_eject_.split()]
    elif str_evo_table == "WW95":
        if element == "H":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_H_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "He":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_He_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "C":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_C_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "N":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_N_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "O":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_O_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Mg":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Mg_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Ne":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Ne_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Si":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Si_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "S":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_S_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Ca":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Ca_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
                'r')
        elif element == "Fe":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Fe_eject_mass_from_WW95/WW95_Z={}.txt'.format(Z_select_in_table_2),
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
                'yield_tables/rearranged/setllar_H_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
                'r')
        elif element == "He":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_He_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
                'r')
        elif element == "C":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_C_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
                'r')
        elif element == "N":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_N_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
                'r')
        elif element == "O":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_O_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
                'r')
        elif element == "Mg":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Mg_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
                'r')
        elif element == "Ne":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Ne_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
                'r')
        elif element == "Si":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Si_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
                'r')
        elif element == "S":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_S_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
                'r')
        elif element == "Ca":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Ca_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
                'r')
        elif element == "Fe":
            file_M_eject = open(
                'yield_tables/rearranged/setllar_Fe_eject_mass_from_marigo01/marigo01_Z={}.txt'.format(Z_select_in_table_3),
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

def get_BH_mass(mass_boundary, mass_grid_table_number, Mtarget_table_number, mass_calibration_factor, steller_mass_upper_bound):
    if mass_boundary < steller_mass_upper_bound:
        BH_mass = function_get_target_mass_in_range(max(mass_boundary, 40), steller_mass_upper_bound, mass_grid_table_number,
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

def function_get_target_mass_in_range(lower_mass_limit, upper_mass_limit, mass_grid_table_number, Mtarget_table_number, mass_calibration_factor):
    integrate_in_range = quad(integrator_for_function_get_target_mass_in_range, lower_mass_limit, upper_mass_limit,
                              (mass_grid_table_number, Mtarget_table_number), limit=40)[0]####################################
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
    if Mtarget_table_number ==1:
        Mtarget_table = Mfinal_table
    if Mtarget_table_number ==2:
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
    target_mass = round((Mtarget_table[x] + (Mtarget_table[x+1]-Mtarget_table[x]) * (initial_mass-mass_grid_table_n[x]) /
                           (mass_grid_table_n[x+1]-mass_grid_table_n[x])), 5)
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
    #     formed_mass = SFH * unit_SFR * 10 ** 7
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
    # a = quad(IMF.imf, 0.08, steller_mass_upper_bound, limit=40)[0]
    # b = quad(IMF.imf, mass_boundary_observe[0], mass_boundary_observe[1], limit=40)[0]
    # for mass_in_range in M_in_range_list:
    #     est_mass = mass_in_range * a / b
    #     if est_mass == 0:
    #         M_tot_est_list += [1]
    #     else:
    #         M_tot_est_list += [est_mass]


def function_get_igimf_for_this_epoch(SFR_input, Z_over_X, this_time, this_epoch, check_igimf):
    # this function calculate igimf, write them in directory Generated_IGIMFs, and import the file
    # with igimf = function_get_igimf_for_every_epoch(SFH_input, unit_SFR, Z, Z_solar),
    # the igimf can be called by: igimf.custom_imf(stellar_mass, this_time).
    function_generate_igimf_file(SFR=SFR_input, Z_over_X=Z_over_X, printout=None, sf_epoch=this_epoch, check=check_igimf)
    if SFR_input == 0:
        igimf_file_name = "igimf_SFR_Zero"
    else:
        igimf_file_name = "igimf_SFR_{}_Fe_over_H_{}".format(round(math.log(SFR_input, 10)*100000), round(Z_over_X * 100000))
    igimf = __import__(igimf_file_name)
    # if shows ModuleNotFoundError:
    # No module named 'igimf_SFR_..._Fe_over_H_...',
    # then try clear all (except for the first and last) lines in the file Generated_IGIMGs/all_igimf_list.txt.
    # This will force the program to generate new IGIMF functions for future use,
    # instead of looking for the IGIMF in the old generated ones.
    return igimf

def function_generate_igimf_file(SFR=None, Z_over_X=None, printout=False, sf_epoch=0, check=False):
    # python3 code, last update Wed 10 July 2017


    # --------------------------------------------------------------------------------------------------------------------------------
    # import modules and libraries
    # --------------------------------------------------------------------------------------------------------------------------------

    import galIMF  # galIMF containing IGIMF function and OSGIMF function and additional computational modules

    import matplotlib.pyplot as plt  # matplotlib for plotting
    import numpy as np
    from scipy.integrate import simps  # numpy and scipi for array operations
    import math
    import time
    import sys


    # --------------------------------------------------------------------------------------------------------------------------------
    # check if the required IGIMF has already been generated
    # --------------------------------------------------------------------------------------------------------------------------------


    check_file = open('Generated_IGIMFs/all_igimf_list.txt', 'r')
    igimf_list = check_file.readlines()
    check_file.close()

    exist = 0

    if check == True:
        i = 0
        while i < len(igimf_list):
            data = [float(x) for x in igimf_list[i].split()]
            if SFR == data[0] and Z_over_X == data[1]:
                exist = 1
                break
            (i) = (i + 1)

    if exist == 0:

        print("Generating new IGIMF...")

        # --------------------------------------------------------------------------------------------------------------------------------
        # add new headline into the list file -- all_igimf_list.txt:
        # --------------------------------------------------------------------------------------------------------------------------------

        check_file = open('Generated_IGIMFs/all_igimf_list.txt', 'r')
        igimf_list = check_file.read()
        check_file.close()

        check_file = open('Generated_IGIMFs/all_igimf_list.txt', 'w')
        new_headline = igimf_list + '{} {}\n'.format(SFR, Z_over_X)
        check_file.write(new_headline)
        check_file.close()

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
            bindw = galIMF.resolution_histogram_relative = 10 ** (max((0 - math.log(SFR, 10)), 0) ** (0.2) - 1.9)
        # will change the resolution of histogram for optimall sampling automatically addjusted with SFR value.
        alpha3_model = 2#1  # IMF high-mass-end power-index model, see Function_alpha_3_change in file 'galIMF.py'
        alpha_2 = 2.3  # IMF middle-mass power-index
        alpha_1 = 1.3  # IMF low-mass-end power-index
        alpha2_model = 1#0  # see file 'galIMF.py'
        alpha1_model = 1#0  # see file 'galIMF.py'
        beta_model = 1
        if Z_over_X is None:
            Z_over_X = float(input("\nPlease input the metallicity, [Z/X] = log(M_{metal}/M_{H})-log(M_{metal,sun}/M_{H,sun})"
                                    "\n\n[Z/X] = ..."))
        # ----------------------------------------------------------------

        # Parameters below are internal parameters of the theory.
        # Read Yan et al. 2017 carefully before change them!

        delta_t = 10.  # star formation epoch [Myr]
        I_ecl = 1.  # normalization factor in the Optimal Sampling condition equation
        M_ecl_U = 10**9 # 10**(0.75 * math.log(SFR, 10) + 4.8269) # Recchi 2009
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

        galIMF.function_galimf(
            "I",  # IorS ### "I" for IGIMF; "OS" for OSGIMF
            SFR,  # Star Formation Rate [solar mass / yr]
            alpha3_model,  # IMF high-mass-end power-index model, see file 'alpha3.py'
            delta_t,  # star formation epoch [Myr]
            Z_over_X,
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

        masses = np.array(galIMF.List_M_str_for_xi_str)
        igimf = np.array(galIMF.List_xi)


        #######################################################
        # generated igimf is normalized by default to a total mass formed in 10 Myr given the SFR,
        # i.e., total stellar mass.
        # to change the normalization, follow the below commented part:
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
                file = open('Generated_IGIMFs/igimf_SFR_{}_Fe_over_H_{}.py'.format(round(math.log(SFR, 10)*100000),
                                                                                   round(Z_over_X * 100000)), 'w')
                file.write("# File to define a custom IMF\n"
                           "# The return value represents the chosen IMF value for the input mass\n\n\n")
                file.write("def custom_imf(mass, time):  # there is no time dependence for IGIMF\n")
                file.write("    if mass < 0.08:\n")
                file.write("        return 0\n")
                file.write("    elif mass < %s:\n" % masses[1])
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

def function_element_abundunce(element_1_name, element_2_name, metal_1_mass, metal_2_mass):
    # this function calculate the atom number ratio compare to solar value [metal/H]
    if metal_2_mass == 0:
        if metal_1_mass == 0:
            metal_1_over_2 = -4
        elif metal_1_mass > 0:
            metal_1_over_2 = None
        elif metal_1_mass < 0:
            print("Warning: current {} mass < 0.".format(element_1_name))
            metal_1_over_2 = -4
    elif metal_2_mass < 0:
        print("Warning: current {} mass < 0.".format(element_2_name))
        if metal_1_mass == 0:
            metal_1_over_2 = -4
        elif metal_1_mass > 0:
            metal_1_over_2 = -4
        elif metal_1_mass < 0:
            print("Warning: current {} mass < 0.".format(element_1_name))
            metal_1_over_2 = None
    else:
        if metal_1_mass == 0:
            metal_1_over_2 = None
        elif metal_1_mass < 0:
            print("Warning: current {} mass < 0.".format(element_1_name))
            metal_1_over_2 = None
        else:
            solar_metal_1_logarithmic_abundances = solar_element_abundances.function_solar_element_abundances(element_1_name)
            solar_metal_2_logarithmic_abundances = solar_element_abundances.function_solar_element_abundances(element_2_name)

            metal_1_element_weight = element_weight_table.function_element_weight(element_1_name)
            metal_2_element_weight = element_weight_table.function_element_weight(element_2_name)

            metal_1_over_2 = math.log(metal_1_mass / metal_2_mass / metal_1_element_weight * metal_2_element_weight, 10) \
                         - (solar_metal_1_logarithmic_abundances - solar_metal_2_logarithmic_abundances)
    return metal_1_over_2

def function_get_avaliable_Z(str_evo_table):
    # extract avalible metallicity in the given grid table
    # stellar life-time table and metal production tables have different avalible metal grid.
    import os

    # list 1
    file_names_setllar_lifetime_from_str_evo_table = os.listdir('yield_tables/rearranged/setllar_lifetime_from_portinari98')
    Z_list = []
    for name in file_names_setllar_lifetime_from_str_evo_table:
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
        Z_list += [float(Z)]
    sorted_Z_list = sorted(Z_list)
    # list 2
    file_names_setllar_lifetime_from_str_evo_table = os.listdir('yield_tables/rearranged/setllar_Metal_eject_mass_from_{}'.format(str_evo_table))
    Z_list_2 = []
    for name in file_names_setllar_lifetime_from_str_evo_table:
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
        Z_list_2 += [float(Z)]
    sorted_Z_list_2 = sorted(Z_list_2)
    if str_evo_table != "portinari98":
        # list 3
        file_names_setllar_lifetime_from_str_evo_table = os.listdir(
            'yield_tables/rearranged/setllar_Metal_eject_mass_from_marigo01')
        Z_list_3 = []
        for name in file_names_setllar_lifetime_from_str_evo_table:
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
            Z_list_3 += [float(Z)]
        sorted_Z_list_3 = sorted(Z_list_3)
    else:
        sorted_Z_list_3 = []
    return (sorted_Z_list, sorted_Z_list_2, sorted_Z_list_3)

def function_select_metal(Z, Z_list):
    if Z <= Z_list[0]:
        Z_select_in_table = Z_list[0]
        return Z_select_in_table
    elif Z >= Z_list[-1]:
        Z_select_in_table = Z_list[-1]
        return Z_select_in_table
    else:
        i = 1
        while i < len(Z_list):
            if Z < Z_list[i]:
                if Z <= (Z_list[i] + Z_list[i - 1]) / 2:
                    Z_select_in_table = Z_list[i - 1]
                    return Z_select_in_table
                else:
                    Z_select_in_table = Z_list[i]
                    return Z_select_in_table
            (i)=(i+1)


def fucntion_mass_boundary(time, mass_grid_for_lifetime, lifetime):
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


def function_xi_Kroupa_IMF(mass):
    # integrate this function's output xi result in the number of stars in mass limits.
    xi = Kroupa_IMF.custom_imf(mass, 0)
    return xi

def function_mass_Kroupa_IMF(mass):
    # integrate this function's output m result in the total stellar mass for stars in mass limits.
    m = mass * Kroupa_IMF.custom_imf(mass, 0)
    return m

def text_output(imf, SFE, SFR, SFEN, original_gas_mass, Z_0, Z_solar):
    global time_axis
    # print("time:", time_axis)

    global all_sf_imf
    number_of_sf_epoch = len(all_sf_imf)

    # data = exec(open("simulation_results_from_galaxy_evol/imf:{}-SFE:{}-log_SFR:{}-SFEN:{}.txt".format(IMF, SFE[0], SFR[0], SFEN[0])).read())
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
    # integrate_Kroupa_stellar_mass_range_1 = quad(function_mass_Kroupa_IMF, mass_range_1[0], mass_range_1[1], limit=50)[
    #     0]
    # integrate_Kroupa_stellar_mass_range_2 = quad(function_mass_Kroupa_IMF, mass_range_2[0], mass_range_2[1], limit=50)[
    #     0]
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
    print("Mass of all alive stars at final time: 10^", stellar_mass)
    downsizing_relation__star_formation_duration = round(10**(2.38-0.24*stellar_mass), 4)  # Recchi 2009
    # print("star formation duration (downsizing relation):", downsizing_relation__star_formation_duration, "Gyr")

    stellar_and_remnant_mass = round(math.log(stellar_mass_list[-1] + remnant_mass_list[-1], 10), 1)
    print("Mass of stars and remnants at final time: 10^", stellar_and_remnant_mass)

    total_mas_in_box = original_gas_mass

    # # Dabringhausen 2008 eq.4
    # Dabringhausen_2008_a = 2.95
    # Dabringhausen_2008_b = 0.596
    # expansion_factor = 5 ################ the expansion_factor should be a funtion of galaxy mass and rise with the mass
    # log_binding_energy = round(
    #     math.log(4.3 * 6 / 5, 10) + 40 + (2 - Dabringhausen_2008_b) * math.log(total_mas_in_box, 10) - math.log(
    #         Dabringhausen_2008_a, 10) + 6 * Dabringhausen_2008_b + math.log(expansion_factor, 10), 1)
    # # print("the gravitational binding energy: 10^", log_binding_energy, "erg")

    global Fe_over_H_list, stellar_Fe_over_H_list
    # print("Gas [Fe/H]:", round(Fe_over_H_list[-1], 3))
    # print("Stellar [Fe/H]:", round(stellar_Fe_over_H_list[-1], 3))

    global Mg_over_Fe_list, stellar_Mg_over_Fe_list
    # print("Gas [Mg/Fe]:", round(Mg_over_Fe_list[-1], 3))
    # print("Stellar [Mg/Fe]:", round(stellar_Mg_over_Fe_list[-1], 3))

    global O_over_Fe_list, stellar_O_over_Fe_list
    # print("Gas [O/Fe]:", round(O_over_Fe_list[-1], 3))
    # print("Stellar [O/Fe]:", round(stellar_O_over_Fe_list[-1], 3))

    global Mg_over_H_list, stellar_Mg_over_H_list
    # print("Gas [Mg/H]:", round(Mg_over_H_list[-1], 3))
    # print("Stellar [Mg/H]:", round(stellar_Mg_over_H_list[-1], 3))

    global O_over_H_list, stellar_O_over_H_list
    # print("Gas [O/H]:", round(O_over_H_list[-1], 3))
    # print("Stellar [O/H]:", round(stellar_O_over_H_list[-1], 3))

    global gas_Z_over_X_list, stellar_Z_over_X_list, stellar_Z_over_H_list
    # print("Gas metallicity:", round(gas_Z_over_X_list[-1], 3))
    # print("Stellar metallicity:", round(stellar_Z_over_X_list[-1], 3))
    # print("Stellar [Z/H]:", round(stellar_Z_over_H_list[-1], 3))

    log_Z_0 = round(math.log(Z_0/Z_solar, 10), 2)
    file = open('simulation_results_from_galaxy_evol/imf:{}-SFE:{}-log_SFR:{}-SFEN:{}-Z_0:{}.txt'.format(imf, SFE, SFR, SFEN, log_Z_0), 'w')

    print("simulation results saved in the file:\n"
          "simulation_results_from_galaxy_evol/imf:{}-SFE:{}-log_SFR:{}-SFEN:{}-Z_0:{}.txt".format(imf, SFE, SFR, SFEN, log_Z_0))

    file.write("# Number of star formation event epoch (10^7 yr):\n")
    file.write("%s\n" % number_of_sf_epoch)

    file.write("# Modeled star formation duration (Gyr):\n")
    file.write("{}\n".format(number_of_sf_epoch/100))

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

    file.write("# Time step list:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % time_axis[i])
        (i) = (i+1)
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

    file.write("# Stellar [Fe/H]:\n")
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

    file.write("# Stellar [Mg/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Mg_over_Fe_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [O/Fe]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % O_over_Fe_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar [O/Fe]:\n")
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

    file.write("# Stellar [Mg/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Mg_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Gas [O/H]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % O_over_H_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Stellar [O/H]:\n")
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

    file.write("# Stellar metallicity, [Z/X]:\n")
    i = 0
    while i < length_of_time_axis:
        file.write("%s " % stellar_Z_over_X_list[i])
        (i) = (i + 1)
    file.write("\n")

    file.write("# Total number of SNIa (log_10):\n")
    file.write("%s\n" % round(math.log(SNIa_energy_release_list[-1], 10), 1))

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

    file.close()

    return

def plot_output(plot_show, plot_save, imf, igimf):
    # plot SFH
    global all_sfr
    SFR_list = []
    age_list = []
    age_list.append(0)
    SFR_list.append(-5)
    for i in range(len(all_sfr)):
        age_list.append(i / 100)
        SFR_list.append(math.log(all_sfr[i][0], 10))
        age_list.append((i+1) / 100)
        SFR_list.append(math.log(all_sfr[i][0], 10))
    age_list.append((i+1) / 100)
    SFR_list.append(-5)
    age_list.append(5)
    SFR_list.append(-5)

    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(0, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(age_list, SFR_list)
        plt.xlabel('Time [Gyr]')
        plt.ylabel(r'log$_{10}({\rm SFR} [{\rm M}_\odot$/yr])')
        # plt.xlim(6.4, 1.01 * log_time_axis[-1])
        # plt.ylim(-5, 1)
        plt.tight_layout()
        # plt.legend()
        if plot_save is True:
            plt.savefig('galaxy_evolution_fig_SFH_{}.pdf'.format(imf), dpi=250)

    length_of_SFH_list = len(SFR_list)
    file = open('simulation_results_from_galaxy_evol/plots/SFH.txt', 'w')
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
            while j < i+1:
                xi_each_time__ += all_sf_imf[j][0].custom_imf(mass, 0)
                (j) = (j+1)
            if xi_each_time__ == 0:
                xi_each_time_log[i] += [-10]
                xi_each_time[i] += [0]
            else:
                xi_each_time_log[i] += [math.log(xi_each_time__, 10)]
                xi_each_time[i] += [xi_each_time__]
            (mass) = (mass * 0.99)
        (i)=(i+1)

    j = 0
    xi_1_last_time = 0
    while j < number_of_sf_epoch:
        xi_1_last_time += all_sf_imf[j][0].custom_imf(1, 0)
        (j) = (j + 1)
    normal = xi_1_last_time / Kroupa_IMF.custom_imf(1, 0)

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
        xi_Kroupa__ = Kroupa_IMF.custom_imf(mass, 0) * normal
        if xi_Kroupa__ == 0:
            xi_Kroupa += [-10]
        else:
            xi_Kroupa += [math.log(xi_Kroupa__, 10)]
        xi_Kroupa_not_log += [xi_Kroupa__]
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
        j = mass_list_length-1
        while j > 0:
            file.write("{} {}\n".format(mass_list[j], xi_each_time[i][j]))
            (j) = (j-1)
        file.close()
        (i) = (i + 1)

    i = 0
    while i < number_of_sf_epoch:
        time = round(all_sf_imf[i][2] / 10 ** 6)
        length_of_xi = len(mass_list)
        file = open('simulation_results_from_galaxy_evol/plots/imf_at_time_{}_Myr.txt'.format(time), 'w')
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
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(1, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(mass_list, xi_Kroupa, linestyle='dashed', color='r', label='Kroupa IMF')
        # i = 0
        # while i < number_of_sf_epoch:
        #     time = round(all_sf_imf[i][2] / 10**6)
        #     plt.plot(mass_list, xi_each_time_log[i], label='IMF at {} Myr'.format(time))
        #     (i) = (i + 1)
        plt.plot(mass_list, xi_observe, label='IGIMF')
        plt.xlabel(r'log$_{10}(M_\star [M_\odot$])')
        plt.ylabel(r'log$_{10}(\xi_\star)$')
        plt.ylim(0, 11)
        plt.legend(prop={'size': 6})
        plt.tight_layout()
        #
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(2, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(mass_list, xi_Kroupa, linestyle='dashed', color='r', label='Kroupa IMF')
        i = 0
        while i < number_of_sf_epoch:
            plt.plot(mass_list, xi_each_epoch[i], label='SF epoch {}'.format(i))
            (i) = (i + 1)
        plt.plot(mass_list, xi_observe, label='final observed IMF')
        plt.xlabel(r'log$_{10}(M_\star)$ [M$_{\odot}$]')
        plt.ylabel(r'$\log_(\xi_\star)$')
        plt.legend(prop={'size': 6})
        plt.tight_layout()

    #
    global time_axis
    length_of_time_axis = len(time_axis)
    log_time_axis = []
    for i in range(length_of_time_axis):
        if time_axis[i] != 0:
            log_time_axis += [math.log((time_axis[i]), 10)]
        else:
            log_time_axis += [0]

    # global mm, zz
    # fig = plt.figure(0, figsize=(4, 3.5))
    # plt.plot(mm, zz)

    global Fe_over_H_list, stellar_Fe_over_H_list, stellar_Fe_over_H_list_luminosity_weighted

    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(3, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(log_time_axis, Fe_over_H_list, label='gas')
        plt.plot(log_time_axis, stellar_Fe_over_H_list, label='stellar MW')
        plt.plot(log_time_axis, stellar_Fe_over_H_list_luminosity_weighted, label='stellar LW')
        plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
        plt.xlabel(r'log$_{10}$(age) [yr]')
        plt.ylabel('[Fe/H]')
        # if imf == 'igimf':
        #     plt.title('IGIMF')
        # elif imf == 'Kroupa':
        #     plt.title('Kroupa IMF')
        # plt.legend(scatterpoints=1, numpoints=1, loc=0, prop={'size': 7.5}, ncol=2)
        # plt.xlim(6.4, 1.01 * log_time_axis[-1])
        # plt.ylim(-5, 1)
        plt.tight_layout()
        plt.legend()
        if plot_save is True:
            plt.savefig('galaxy_evolution_fig_FeH_{}.pdf'.format(imf), dpi=250)

    file = open('simulation_results_from_galaxy_evol/plots/Fe_over_H_time.txt', 'w')
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
            file.write("{} ".format(-4))
        else:
            file.write("{} ".format(stellar_Fe_over_H_list[i]))
        (i) = (i + 1)
    file.write("\n# stellar_Fe_over_H_list_luminosity_weighted\n")
    i = 0
    while i < length_of_time_axis:
        if stellar_Fe_over_H_list_luminosity_weighted[i] is None:
            file.write("{} ".format(-4))
        else:
            file.write("{} ".format(stellar_Fe_over_H_list_luminosity_weighted[i]))
        (i) = (i + 1)
    file.write("\n")
    file.close()

    file = open('simulation_results_from_galaxy_evol/plots/Fe_over_H_mass.txt', 'w')
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
    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(4, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(log_time_axis, O_over_H_list, label='gas')
        plt.plot(log_time_axis, stellar_O_over_H_list, label='stellar MW')
        plt.plot(log_time_axis, stellar_O_over_H_list_luminosity_weighted, label='stellar LW')
        plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
        plt.xlabel(r'log$_{10}$(age) [yr]')
        plt.ylabel('[O/H]')
        # plt.xlim(6.4, 1.01 * log_time_axis[-1])
        # plt.ylim(-5, 1)
        plt.tight_layout()
        plt.legend()
        if plot_save is True:
            plt.savefig('galaxy_evolution_fig_OH_{}.pdf'.format(imf), dpi=250)
    #
    global Mg_over_H_list, stellar_Mg_over_H_list, stellar_Mg_over_H_list_luminosity_weighted
    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(5, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(log_time_axis, Mg_over_H_list, label='gas')
        # print(stellar_Mg_over_H_list)
        plt.plot(log_time_axis, stellar_Mg_over_H_list, label='stellar MW')
        plt.plot(log_time_axis, stellar_Mg_over_H_list_luminosity_weighted, label='stellar LW')
        plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
        plt.xlabel(r'log$_{10}$(age) [yr]')
        plt.ylabel('[Mg/H]')
        # plt.xlim(6.4, 1.01 * log_time_axis[-1])
        # plt.ylim(-5, 1)
        plt.tight_layout()
        plt.legend()
        if plot_save is True:
            plt.savefig('galaxy_evolution_fig_MgH_{}.pdf'.format(imf), dpi=250)
    #
    global stellar_Y_list, Y_list, Y_list_luminosity_weighted
    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(55, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(log_time_axis, Y_list, label='gas')
        plt.plot(log_time_axis, stellar_Y_list, label='stellar MW')
        plt.plot(log_time_axis, stellar_Y_list_luminosity_weighted, label='stellar LW')
        plt.plot([log_time_axis[0], log_time_axis[-1]], [0.2485, 0.2485], color='red', ls='dashed', label='solar')
        plt.xlabel(r'log$_{10}$(age) [yr]')
        plt.ylabel('Y')
        plt.tight_layout()
        plt.legend()
        if plot_save is True:
            plt.savefig('galaxy_evolution_fig_Y_{}.pdf'.format(imf), dpi=250)

    file = open('simulation_results_from_galaxy_evol/plots/Y_time.txt', 'w')
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
    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(7, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(log_time_axis, Mg_over_Fe_list, label='gas')
        plt.plot(log_time_axis, stellar_Mg_over_Fe_list, label='stellar MW')
        plt.plot(log_time_axis, stellar_Mg_over_Fe_list_luminosity_weighted, label='stellar LW')
        plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
        plt.xlabel(r'log$_{10}$(age) [yr]')
        plt.ylabel('[Mg/Fe]')
        # plt.xlim(6.4, 1.01 * log_time_axis[-1])
        # plt.ylim(-1, 3.5)
        plt.tight_layout()
        plt.legend()
        if plot_save is True:
            plt.savefig('galaxy_evolution_fig_MgFe_{}.pdf'.format(imf), dpi=250)

    file = open('simulation_results_from_galaxy_evol/plots/Mg_over_Fe_time.txt', 'w')
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
            file.write("-4 ")
        else:
            file.write("{} ".format(stellar_Mg_over_Fe_list[i]))
        (i) = (i + 1)
    file.write("\n# stellar_Mg_over_Fe_list_luminosity_weighted\n")
    i = 0
    while i < length_of_time_axis:
        if stellar_Mg_over_Fe_list_luminosity_weighted[i] is None:
            file.write("-4 ")
        else:
            file.write("{} ".format(stellar_Mg_over_Fe_list_luminosity_weighted[i]))
        (i) = (i + 1)
    file.write("\n")
    file.close()

    file = open('simulation_results_from_galaxy_evol/plots/Mg_over_Fe_mass.txt', 'w')
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
    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(8, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(log_time_axis, O_over_Fe_list, label='gas')
        plt.plot(log_time_axis, stellar_O_over_Fe_list, label='stellar MW')
        plt.plot(log_time_axis, stellar_O_over_Fe_list_luminosity_weighted, label='stellar LW')
        plt.plot([log_time_axis[0], log_time_axis[-1]], [0, 0], color='red', ls='dashed', label='solar')
        plt.xlabel(r'log$_{10}$(age) [yr]')
        plt.ylabel('[O/Fe]')
        # plt.xlim(6.4, 1.01 * log_time_axis[-1])
        # plt.ylim(-1, 3.5)
        plt.tight_layout()
        plt.legend()
        if plot_save is True:
            plt.savefig('galaxy_evolution_fig_MgFe_{}.pdf'.format(imf), dpi=250)
    #
    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(9, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(Fe_over_H_list, Mg_over_Fe_list, label='gas')
        plt.plot(stellar_Fe_over_H_list, stellar_Mg_over_Fe_list, label='stellar MW')
        plt.plot(stellar_Fe_over_H_list, stellar_Mg_over_Fe_list_luminosity_weighted, label='stellar LW')
        plt.plot([-5, 1], [0, 0], color='red', ls='dashed', label='solar')
        plt.plot([0, 0], [-1, 3.5], color='red', ls='dashed')
        plt.xlabel('[Fe/H]')
        plt.ylabel('[Mg/Fe]')
        # plt.xlim(-5, 1)
        # plt.ylim(-1, 3.5)
        plt.legend()
        plt.tight_layout()
        if plot_save is True:
            plt.savefig('galaxy_evolution_fig_MgFe-FeH_{}.pdf'.format(imf), dpi=250)
    #
    global SNIa_number_per_century, SNII_number_per_century
    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(10, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.loglog(time_axis, SNIa_number_per_century, label='SNIa')  # Number per century
        plt.loglog(time_axis, SNII_number_per_century, label='SNII')  # Number per century
        plt.loglog(time_axis, SN_number_per_century, ls="dotted", label='total')
        plt.xlabel(r'Time [yr]')
        plt.ylabel(r'Number of SN per century')
        # plt.xlim(6, 1.01 * log_time_axis[-1])
        # plt.ylim(8.5, 11.6)
        plt.legend()
        plt.tight_layout()
        if plot_save is True:
            plt.savefig('galaxy_evolution_SN_number.pdf'.format(imf), dpi=250)

    file = open('simulation_results_from_galaxy_evol/plots/SN_number_evolution.txt', 'w')
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
    final_crossing_time = 1 / (gravitational_constant * 10 ** (-11) * finial_galaxy_inner_mass * 2 * 10 ** 30 / (final_galaxy_radii * 3 * 10 ** 16) ** 3) ** (0.5) / (60 * 60 * 24 * 365 * 10 ** 6)  # in Myr

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
        crossing_time = 1/(gravitational_constant*10**(-11)*current_galaxy_inner_mass*2*10**30/(current_galaxy_radii*3*10**16)**3)**(0.5)/(60*60*24*365*10**6) # in Myr
        SN_energy_per_current_crossing_time = SN_number_per_century[i] * crossing_time * 10**4 * 10**51
        SN_energy_per_final_crossing_time = SN_number_per_century[i] * final_crossing_time * 10**4 * 10**51
        SN_energy_per_current_crossing_time_list.append(SN_energy_per_current_crossing_time)
        SN_energy_per_final_crossing_time_list.append(SN_energy_per_final_crossing_time)
        (i) = (i + 1)

    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(11, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.loglog(time_axis, SN_energy_per_current_crossing_time_list, label='SN/CC')
        plt.loglog(time_axis, SN_energy_per_final_crossing_time_list, label='SN/FC')
        plt.loglog(time_axis, SNIa_energy_release_list, label='tot SNIa')
        plt.loglog(time_axis, SNII_energy_release_list, label='tot SNII')
        plt.loglog(time_axis, total_energy_release_list, ls="dotted", label='tot SNIa+II')
        plt.loglog(time_axis, binding_energy_list, label='binding')
        plt.loglog(time_axis, total_gas_kinetic_energy_list, label='gas kinetic')
        plt.xlabel(r'Time [yr]')
        plt.ylabel(r'Energy')
        # plt.xlim(6, 1.01 * log_time_axis[-1])
        # plt.ylim(8.5, 11.6)
        plt.legend()
        plt.tight_layout()

    file = open('simulation_results_from_galaxy_evol/plots/energy_evolution.txt', 'w')
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

    file = open('simulation_results_from_galaxy_evol/plots/energy_mass.txt', 'w')
    file.write("# final SN_energy_release\n")
    file.write("{}\n".format(total_energy_release_list[-1]))
    file.write("# final binding_energy\n")
    file.write("{}\n".format(binding_energy_list[-1]))
    file.write("# final gas_kinetic_energy\n")
    file.write("{}\n".format(total_gas_kinetic_energy_list[-1]))
    file.close()

    #
    global gas_Z_over_X_list, stellar_Z_over_X_list, stellar_Z_over_X_list_luminosity_weighted
    if plot_show is True or plot_save is True:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(12, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        # time_axis[0] = 1
        # time_axis_G = [0]*length_of_time_axis
        # for i in range(length_of_time_axis):
        #     time_axis_G[i] = time_axis[i]/10**9
        # gas_Z_over_X_list[i]=math.log(gas_Z_over_X_list[i], 10)
        plt.plot(log_time_axis, gas_Z_over_X_list, label='gas')
        plt.plot(log_time_axis, stellar_Z_over_X_list, label='stellar MW')
        plt.plot(log_time_axis, stellar_Z_over_X_list_luminosity_weighted, label='stellar LW')
        plt.plot([6, 11], [0, 0], color='red', ls='dashed', label='solar')
        # The [Z/X]s where the applied portinari98 stellar yield table will be changed for Z=0.0127, 0.008, 0.004, 0.0004.
        plt.plot([6, 11], [-1.173, -1.173], color='red', ls='dotted', lw=0.5)
        plt.plot([6, 11], [-0.523, -0.523], color='red', ls='dotted', lw=0.5)
        plt.plot([6, 11], [-0.272, -0.272], color='red', ls='dotted', lw=0.5)
        plt.xlabel(r'log(age [Gyr])')
        plt.ylabel('[Z/X]')
        # plt.ylim(-2, 1)
        # if imf == 'igimf':
        #     plt.title('IGIMF')
        # elif imf == 'Kroupa':
        #     plt.title('Kroupa IMF')
        # plt.legend(scatterpoints=1, numpoints=1, loc=0, prop={'size': 7.5}, ncol=2)
        # plt.xlim(6.4, 1.01*time_axis[-1])
        # plt.ylim(-0.4, 0.2)
        plt.legend()
        plt.tight_layout()
        if plot_save is True:
            plt.savefig('galaxy_evolution_fig_Z_{}.pdf'.format(imf), dpi=250)

    for i in range(length_of_time_axis):
        time_axis[i] = math.log(time_axis[i], 10)

    file = open('simulation_results_from_galaxy_evol/plots/Z_over_X_time.txt', 'w')
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

    file = open('simulation_results_from_galaxy_evol/plots/Z_over_X_mass.txt', 'w')
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
        remnant_mass_list[i] = math.log(remnant_mass_list[i], 10)
        total_gas_mass_list[i] = math.log(total_gas_mass_list[i], 10)
        stellar_mass_list[i] = math.log(stellar_mass_list[i], 10)
        ejected_gas_mass_list[i] = math.log(ejected_gas_mass_list[i], 10)
        if WD_mass_list[i]<0:
            print("Warrning: WD_mass_list[{}] < 0".format(i))
            WD_mass_list[i]=1
        else:
            WD_mass_list[i] = math.log(WD_mass_list[i], 10)
        NS_mass_list[i] = math.log(NS_mass_list[i], 10)
        BH_mass_list[i] = math.log(BH_mass_list[i], 10)
    # time_axis[0] = time_axis[1]
    if plot_show is True or plot_save is True:
        fig = plt.figure(13, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(time_axis, total_gas_mass_list, lw=2, label='all gas')
        plt.plot(time_axis, ejected_gas_mass_list, lw=2, label='ejected gas')
        plt.plot(time_axis, stellar_mass_list, lw=2, label='alive stars')
        plt.plot(time_axis, remnant_mass_list, lw=2, label='all remnants')
        plt.plot(time_axis, BH_mass_list, lw=2, label='black holes')
        plt.plot(time_axis, NS_mass_list, lw=2, label='neutron stars')
        plt.plot(time_axis, WD_mass_list, lw=2, label='white dwarfs')

    file = open('simulation_results_from_galaxy_evol/plots/mass_evolution.txt', 'w')
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
    final_alive_and_remnant_stellar_mass = math.log((10**final_alive_stellar_mass+10**final_remnant_stellar_mass), 10)
    file = open('simulation_results_from_galaxy_evol/plots/mass_ratio.txt', 'w')
    file.write("# final alive stellar mass in log\n")
    file.write("{}\n".format(final_alive_stellar_mass))
    file.write("# final alive + remnant mass in log\n")
    file.write("{}\n".format(final_alive_and_remnant_stellar_mass))
    file.write("# (alive + remnant) / alive IN LOG\n")
    file.write("{}\n".format(final_alive_and_remnant_stellar_mass-final_alive_stellar_mass))
    file.close()

    global total_star_formed
    file = open('simulation_results_from_galaxy_evol/plots/SN_number_mass.txt', 'w')
    file.write("# final SNIa_number per stellar mass formed\n")
    file.write("{}\n".format(SNIa_number_list[-1]/total_star_formed))
    print(SNIa_number_list[-1]/total_star_formed)
    file.write("# final SNII_number per stellar mass formed\n")
    file.write("{}\n".format(SNII_number_list[-1]/total_star_formed))
    file.close()

    if plot_show is True or plot_save is True:
        plt.xlabel(r'log$_{10}$(age [yr])')
        plt.ylabel(r'log$_{10}$(Mass [M$_\odot$])')
        # if imf == 'igimf':
        #     plt.title('IGIMF')
        # elif imf == 'Kroupa':
        #     plt.title('Kroupa IMF')
        plt.legend(prop={'size': 7.5}, loc='lower right')
        # plt.xlim(6.4, 1.01 * time_axis[-1])
        plt.xlim(6.4, 10.1)
        # plt.ylim(6.4, 10.6)
        plt.tight_layout()

    global expansion_factor_instantaneous_list, expansion_factor_adiabat_list
    if plot_show is True or plot_save is True:
        fig = plt.figure(14, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(time_axis, expansion_factor_instantaneous_list, label='instantaneous')
        plt.plot(time_axis, expansion_factor_adiabat_list, label='slow')
        plt.plot(time_axis, expansion_factor_list, label='average')
        plt.xlabel(r'log$_{10}$(age) [yr]')
        plt.ylabel(r'expansion_factor')
        plt.legend(prop={'size': 7.5}, loc='lower right')
        # plt.xlim(6.4, 1.01 * time_axis[-1])
        # plt.ylim(7.3, 12.2)
        # plt.ylim(6, 12)
        plt.tight_layout()

    file = open('simulation_results_from_galaxy_evol/plots/expansion_factor.txt', 'w')
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

    global ejected_gas_Mg_over_Fe_list, instant_ejected_gas_Mg_over_Fe_list
    if plot_show is True or plot_save is True:
        fig = plt.figure(15, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(time_axis, ejected_gas_Mg_over_Fe_list, label='total')
        plt.plot(time_axis, instant_ejected_gas_Mg_over_Fe_list, label='instant')
        plt.xlabel(r'log$_{10}$(age) [yr]')
        plt.ylabel(r'[Mg/Fe]')
        # plt.xlim(6.4, 1.01 * time_axis[-1])
        # plt.ylim(7.3, 12.2)
        # plt.ylim(6, 12)
        plt.legend(prop={'size': 7.5}, loc='best')
        plt.tight_layout()

    global ejected_metal_mass_list
    if plot_show is True or plot_save is True:
        fig = plt.figure(16, figsize=(4, 3.5))
        fig.add_subplot(1, 1, 1)
        plt.plot(time_axis, ejected_metal_mass_list, label='total')
        plt.xlabel(r'log$_{10}$(age) [yr]')
        plt.ylabel(r'ejected metal mass')
        # plt.xlim(6.4, 1.01 * time_axis[-1])
        # plt.ylim(7.3, 12.2)
        # plt.ylim(6, 12)
        plt.legend(prop={'size': 7.5}, loc='best')
        plt.tight_layout()

    if plot_show is True:
        plt.show()
    return


def generate_SFH(distribution, Log_SFR, SFEN):
    if distribution == "skewnorm":
        generate_sfh_skewnorm(Log_SFR, SFEN)
    elif distribution == "flat":
        generate_sfh_flat(Log_SFR, SFEN)
    return


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
    while j < 1001 - SFEN:
        file.write("0\n")
        (j) = (j + 1)
    file.write("# The value in each line stand for the SFR [solar mass / yr]\n")
    file.write("# SFR = the value * unit_SFR (given in file galaxy_evol.py)\n")
    file.write("# in a star formation epoch (10 Myr)\n")
    file.write("# start from time 0 for the first line.\n")
    file.write("# Warning! Effective line number must be larger than 1.\n")
    file.write("# Add a '0' in the next line if there is only one line.\n")

    file.close()
    return


def generate_sfh_skewnorm(Log_SFR, SFEN):
    global sfr_tail
    tot_sf_set = 10 ** Log_SFR * SFEN
    tot_sf = 0
    while tot_sf < tot_sf_set:
        SFEN += 1
        result_cal_tot_sf = cal_tot_sf(Log_SFR, SFEN)
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
        if sfr_for_this_epoch > 10 ** Log_SFR/2:
            result_starburst_sf += sfr_for_this_epoch
        else:
            result_tail_sf += sfr_for_this_epoch
        (j) = (j + 1)
    sfr_for_the_tail_epoch = sfr_for_this_epoch / 2
    if sfr_tail == 0:
        j = 0
        while j < 1001 - SFEN:
            file.write("0\n")
            (j) = (j + 1)
    elif sfr_tail == 1:
        j = 0
        while j < 101:
            file.write("{}\n".format(sfr_for_the_tail_epoch))
            result_tail_sf += sfr_for_the_tail_epoch
            (j) = (j + 1)
        while j < 1001 - SFEN:
            file.write("0\n")
            (j) = (j + 1)
    file.write("# The value in each line stand for the SFR [solar mass / yr]\n")
    file.write("# SFR = the value * unit_SFR (given in file galaxy_evol.py)\n")
    file.write("# in a star formation epoch (10 Myr)\n")
    file.write("# start from time 0 for the first line.\n")
    file.write("# Warning! Effective line number must be larger than 1.\n")
    file.write("# Add a '0' in the next line if there is only one line.\n")

    if sfr_tail == 1:
        print("star formation tail (after the SFR is lower than half of the maximum value) contributes {}% "
              "of the total star formation.".format(round(result_tail_sf/(result_starburst_sf+result_tail_sf)*100, 2)))

    file.close()
    return


def cal_tot_sf(SFR, SFEN):
    # Skew normal distribution for star formation history
    # took input: maximum star formation rate, star formation event number
    # from scipy.stats import f
    global skewness, location
    from scipy.stats import skewnorm
    x = np.linspace(skewnorm.ppf(0.01, skewness, location, 1), skewnorm.ppf(0.99, skewness, location, 1), SFEN)
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
    Log_SFR = 3
    SFEN = 10
    location = 0
    skewness = 10
    sfr_tail = 0
    # generate_SFH("flat", Log_SFR, SFEN)  # "skewnorm" or "flat"
    # galaxy_evol(unit_SFR=1e5, Z_0=0.012, IMF_name='Salpeter', steller_mass_upper_bound=150, time_resolution_in_Myr=1,
    #                  mass_boundary_observe_low=0.5, mass_boundary_observe_up=8)
    # stellar evolution table being "WW95" or "portinari98"
    # imf='igimf' or 'diet_Salpeter'
    galaxy_evol(imf='igimf', unit_SFR=1, SFE=0.45, SFEN=SFEN, Z_0=0.00000001886, Z_solar=0.01886,
                str_evo_table='portinari98', IMF_name='Kroupa', steller_mass_upper_bound=150,
                time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
                SNIa_ON=True, high_time_resolution=None, plot_show=None, plot_save=None, outflow=None, check_igimf=True)
