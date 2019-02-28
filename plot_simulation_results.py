import matplotlib.pyplot as plt
import math
import numpy as np


def plot_galaxy_evol():
    # Figure 1 SFH

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/SFH.txt', 'r')
    SFH = file.readlines()
    file.close()
    age_list = [float(x) for x in SFH[1].split()]
    SFR_list = [10**float(x) for x in SFH[3].split()]
    age_list = [-1] + age_list
    SFR_list = [0] + SFR_list

    # file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/SFH.txt', 'r')
    # SFH0 = file.readlines()
    # file.close()
    # age_list0 = [float(x) for x in SFH0[1].split()]
    # SFR_list0 = [10**float(x) for x in SFH0[3].split()]
    # age_list0 = [-1] + age_list0
    # SFR_list0 = [0] + SFR_list0

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(1, figsize=(4, 2.5))
    fig.add_subplot(1, 1, 1)
    plt.plot(age_list, SFR_list, c='k')
    # plt.plot(age_list, SFR_list, ls='dashed', label='Flat')
    # plt.plot(age_list0, SFR_list0, label='Skewed-normal')
    plt.xlabel('Time [Gyr]')
    plt.ylabel(r'log$_{10}({\rm SFR} [{\rm M}_\odot$/yr])')
    plt.legend(prop={'size': 7}, loc='best')
    plt.xlim(-0.3, 2.6)
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/SFH.pdf', dpi=250)

    # Figure 6 1000_IMF_evolution

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/imf_at_time_10_Myr.txt', 'r')
    imf_at_time_ = file.readlines()
    file.close()
    mass_list = [math.log(float(x), 10) for x in imf_at_time_[1].split()]
    xi_Kroupa = [math.log(float(x)+1e-10, 10) for x in imf_at_time_[5].split()]
    xi_each_time = []
    time_length = 100
    i = time_length
    while i > 0:
        file = open('simulation_results_from_galaxy_evol/plots_IGIMF/imf_at_time_{}0_Myr.txt'.format(i), 'r')
        imf_at_time_ = file.readlines()
        file.close()
        xi_each_time.append([math.log(float(x)+1e-10, 10) for x in imf_at_time_[3].split()])
        (i) = (i-1)

    colors = plt.cm.rainbow(np.linspace(0, 1, time_length))

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(6, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    i = time_length-1
    while i > -1:
        plt.plot(mass_list, xi_each_time[i], lw=2, color=colors[time_length-1-i])
        (i) = (i - 1)
    plt.plot(mass_list, xi_Kroupa, lw=2, c='k', label='Kroupa IMF', ls='dashed')
    # plt.plot([], [])
    plt.xlabel(r'log$_{10}$(M$_*$ [M$_\odot$])')
    plt.ylabel(r'log$_{10}$($\xi$ [M$_\odot^{-1}$])')
    plt.legend(prop={'size': 7}, loc='best')
    # plt.xlim(6.7, 10.1)
    plt.ylim(4.5, 13.5)
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/IMF_evolution.pdf', dpi=250)


    # Figure 7 mass_evolution

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/mass_evolution.txt', 'r')
    mass_evolution = file.readlines()
    file.close()
    time_axis = [float(x) for x in mass_evolution[1].split()]
    total_gas_mass_list = [float(x) for x in mass_evolution[3].split()]
    ejected_gas_mass_list = [float(x) for x in mass_evolution[5].split()]
    stellar_mass_list = [float(x) for x in mass_evolution[7].split()]
    remnant_mass_list = [float(x) for x in mass_evolution[9].split()]

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/mass_evolution.txt', 'r')
    mass_evolution_fixedIMF = file.readlines()
    file.close()
    time_axis_fixedIMF = [float(x) for x in mass_evolution_fixedIMF[1].split()]
    total_gas_mass_list_fixedIMF = [float(x) for x in mass_evolution_fixedIMF[3].split()]
    ejected_gas_mass_list_fixedIMF = [float(x) for x in mass_evolution_fixedIMF[5].split()]
    stellar_mass_list_fixedIMF = [float(x) for x in mass_evolution_fixedIMF[7].split()]
    remnant_mass_list_fixedIMF = [float(x) for x in mass_evolution_fixedIMF[9].split()]

    # , color="tab:blue"
    # , color="tab:orange"

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(7, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)

    plt.plot([], [], color="k", lw=1, label='alive stars')
    plt.plot([], [], ls='dashed', color="k", lw=2, label='all remnants')
    plt.plot([], [], ls='dotted', color="k", lw=3, label='all gas')
    plt.plot([], [], ls='dotted', color="k", lw=1, label='ejected gas')


    plt.plot(time_axis, total_gas_mass_list, lw=3, ls='dotted', color="tab:orange")
    plt.plot(time_axis, ejected_gas_mass_list, lw=1, ls='dotted', color="tab:orange")
    plt.plot(time_axis, stellar_mass_list, lw=1, color="tab:orange", label='IGIMF')
    plt.plot(time_axis, remnant_mass_list, lw=2, ls='dashed', color="tab:orange")

    plt.plot(time_axis_fixedIMF, total_gas_mass_list_fixedIMF, lw=3, ls='dotted', color="tab:blue")
    plt.plot(time_axis_fixedIMF, ejected_gas_mass_list_fixedIMF, lw=1, ls='dotted', color="tab:blue")
    plt.plot(time_axis_fixedIMF, stellar_mass_list_fixedIMF, lw=1, color="tab:blue", label='Fixed IMF')
    plt.plot(time_axis_fixedIMF, remnant_mass_list_fixedIMF, lw=2, ls='dashed', color="tab:blue")
    plt.xlabel(r'log$_{10}$(age [yr])')
    plt.ylabel(r'log$_{10}$(Mass [M$_\odot$])')
    plt.legend(prop={'size': 7}, loc='best')
    plt.xlim(7, 10.1)
    plt.ylim(8, 12.3)
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/mass_evolution.pdf', dpi=250)

    # Figure 8 SN_number_evolution

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/SN_number_evolution.txt', 'r')
    SN_number_evolution = file.readlines()
    file.close()
    time_axis = [float(x) for x in SN_number_evolution[1].split()]
    SNIa_number_per_century = [float(x) for x in SN_number_evolution[3].split()]
    SNII_number_per_century = [float(x) for x in SN_number_evolution[5].split()]

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/SN_number_evolution.txt', 'r')
    SN_number_evolution_fixedIMF = file.readlines()
    file.close()
    time_axis_fixedIMF = [float(x) for x in SN_number_evolution_fixedIMF[1].split()]
    SNIa_number_per_century_fixedIMF = [float(x) for x in SN_number_evolution_fixedIMF[3].split()]
    SNII_number_per_century_fixedIMF = [float(x) for x in SN_number_evolution_fixedIMF[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(8, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], label='SNIa', color="k", lw=2.5)
    plt.plot([], [], label='SNII', color="k", lw=0.7)
    plt.loglog(time_axis_fixedIMF, SNIa_number_per_century_fixedIMF, color="tab:blue", lw=2.5)  # Number per century
    plt.loglog(time_axis_fixedIMF, SNII_number_per_century_fixedIMF, color="tab:blue", label='Fixed IMF', lw=1)  # Number per century
    plt.loglog(time_axis, SNIa_number_per_century, color="tab:orange", lw=2.5)  # Number per century
    plt.loglog(time_axis, SNII_number_per_century, color="tab:orange", label='IGIMF', lw=1)  # Number per century
    plt.xlabel(r'Age [yr]')
    plt.ylabel(r'# of SN per century')
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/SN_number_evolution.pdf', dpi=250)

    # Figure 42 energy_evolution

    file = open('simulation_results_from_galaxy_evol/plots/energy_evolution.txt', 'r')
    energy_evolution = file.readlines()
    file.close()

    time_axis = [float(x) for x in energy_evolution[1].split()]
    SN_energy_per_current_crossing_time_list = [float(x) for x in energy_evolution[7].split()]
    SN_energy_per_final_crossing_time_list = [float(x) for x in energy_evolution[9].split()]
    total_energy_release_list = [float(x) for x in energy_evolution[11].split()]
    binding_energy_list = [float(x) for x in energy_evolution[13].split()]
    total_gas_kinetic_energy_list = [float(x) for x in energy_evolution[15].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(42, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.loglog(time_axis, SN_energy_per_current_crossing_time_list, lw=2, label='SN/CC')
    plt.loglog(time_axis, SN_energy_per_final_crossing_time_list, lw=2, label='SN/FC')
    plt.loglog(time_axis, total_energy_release_list, lw=2, label='Total SN')
    plt.loglog(time_axis, binding_energy_list, lw=2, ls="dashed", label='Gravitational binding')
    plt.loglog(time_axis, total_gas_kinetic_energy_list, lw=2, ls="dotted", label='Gas kinetic')
    plt.xlabel(r'Age [yr]')
    plt.ylabel(r'Energy [erg]')
    plt.legend(prop={'size': 6}, loc='best')
    plt.tight_layout()

    # # Figure 14 expansion_factor
    #
    # file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/expansion_factor.txt', 'r')
    # expansion_factor = file.readlines()
    # file.close()
    #
    # time_axis = [float(x) for x in expansion_factor[1].split()]
    # expansion_factor_list = [math.log(float(x), 10) for x in expansion_factor[3].split()]
    # expansion_factor_instantaneous_list = [math.log(float(x), 10) for x in expansion_factor[5].split()]
    # expansion_factor_adiabatic_list = [math.log(float(x), 10) for x in expansion_factor[7].split()]
    #
    # file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/expansion_factor.txt', 'r')
    # expansion_factor0 = file.readlines()
    # file.close()
    #
    # time_axis0 = [float(x) for x in expansion_factor0[1].split()]
    # expansion_factor0_list = [math.log(float(x), 10) for x in expansion_factor0[3].split()]
    # expansion_factor0_instantaneous_list = [math.log(float(x), 10) for x in expansion_factor0[5].split()]
    # expansion_factor0_adiabatic_list = [math.log(float(x), 10) for x in expansion_factor0[7].split()]
    #
    #
    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(14, figsize=(4, 3.5))
    # fig.add_subplot(1, 1, 1)
    # plt.plot(time_axis, expansion_factor_list, color="tab:red")
    # plt.plot(time_axis0, expansion_factor0_list, color="tab:blue")
    # plt.plot([], [], label='Log-average', color="k")
    # plt.plot(time_axis, expansion_factor_instantaneous_list, ls="dashed", color="tab:red")
    # plt.plot(time_axis0, expansion_factor0_instantaneous_list, ls="dashed", color="tab:blue")
    # plt.plot([], [], ls="dashed", label='Instantaneous', color="k")
    # plt.plot(time_axis, expansion_factor_adiabatic_list, ls="dotted", color="tab:red")
    # plt.plot(time_axis0, expansion_factor0_adiabatic_list, ls="dotted", color="tab:blue")
    # plt.plot([], [], ls="dotted", label='Adiabatic', color="k")
    # plt.plot([], [], lw=3, label=r'SFR=10$^3$ M$_\odot$/yr', color="tab:red")
    # plt.plot([], [], lw=3, label=r'SFR=1 M$_\odot$/yr', color="tab:blue")
    # plt.xlabel(r'log$_{10}$(Age [yr])')
    # plt.ylabel(r'log$_{10}$(Expansion factor)')
    # plt.xlim(6, time_axis[-1])
    # plt.legend(prop={'size': 7}, loc='best')
    # plt.tight_layout()

    # Figure 9 Z_over_X_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Z_over_X_time.txt', 'r')
    Z_over_X_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Z_over_X_time[1].split()]
    gas_Z_over_X_list = [float(x) for x in Z_over_X_time[3].split()]
    stellar_Z_over_X_list = [float(x) for x in Z_over_X_time[5].split()]

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Z_over_X_time.txt', 'r')
    Z_over_X_time_fixedIMF = file.readlines()
    file.close()
    log_time_axis_fixedIMF = [float(x) for x in Z_over_X_time_fixedIMF[1].split()]
    gas_Z_over_X_list_fixedIMF = [float(x) for x in Z_over_X_time_fixedIMF[3].split()]
    stellar_Z_over_X_list_fixedIMF = [float(x) for x in Z_over_X_time_fixedIMF[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(9, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], ls='dashed', label='Gas', color="k")
    plt.plot([], [], label='Stellar', color="k")
    plt.plot(log_time_axis, gas_Z_over_X_list, color="tab:orange", ls='dashed', lw=2)
    plt.plot(log_time_axis, stellar_Z_over_X_list, color="tab:orange", lw=2, label='IGIMF')
    plt.plot(log_time_axis_fixedIMF, gas_Z_over_X_list_fixedIMF, color="tab:blue", ls='dashed', lw=2)
    plt.plot(log_time_axis_fixedIMF, stellar_Z_over_X_list_fixedIMF, color="tab:blue", lw=2, label='Fixed IMF')
    plt.plot([7, 10.1], [0, 0], color='red', ls='dotted', label='Solar')
    # The [Z/X]s where the applied portinari98 stellar yield table will be changed for Z=0.0127, 0.008, 0.004, 0.0004.
    plt.plot([7, 10.1], [-1.173, -1.173], color='k', lw=0.4, label='Yield Switch')
    plt.plot([7, 10.1], [-0.523, -0.523], color='k', lw=0.4)
    plt.plot([7, 10.1], [-0.272, -0.272], color='k', lw=0.4)
    plt.xlabel(r'log$_{10}$(age [yr])')
    plt.ylabel('[Z/X]')
    plt.xlim(7, 10.1)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Z_over_X_time.pdf', dpi=250)

    # Figure 10 Y_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Y_time.txt', 'r')
    Y_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Y_time[1].split()]
    gas_Y_list = [float(x) for x in Y_time[3].split()]
    stellar_Y_list = [float(x) for x in Y_time[5].split()]

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Y_time.txt', 'r')
    Y_time_fixedIMF = file.readlines()
    file.close()
    log_time_axis_fixedIMF = [float(x) for x in Y_time_fixedIMF[1].split()]
    gas_Y_list_fixedIMF = [float(x) for x in Y_time_fixedIMF[3].split()]
    stellar_Y_list_fixedIMF = [float(x) for x in Y_time_fixedIMF[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(10, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], ls='dashed', label='Gas', color="k")
    plt.plot([], [], label='Stellar', color="k")
    plt.plot(log_time_axis, gas_Y_list, color='tab:orange', ls='dashed', lw=2)
    plt.plot(log_time_axis, stellar_Y_list, color='tab:orange', lw=2, label='IGIMF')
    plt.plot(log_time_axis_fixedIMF, gas_Y_list_fixedIMF, color="tab:blue", ls='dashed', lw=2)
    plt.plot(log_time_axis_fixedIMF, stellar_Y_list_fixedIMF, color="tab:blue", lw=2, label='Fixed IMF')
    plt.plot([6.5, 10.2], [0.2485, 0.2485], color='red', ls='dotted', label='solar')
    plt.xlabel(r'log$_{10}$(age [yr])')
    plt.ylabel('Y')
    plt.xlim(7.5, 10.1)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Y_time.pdf', dpi=250)

    # Figure 11 Fe_over_H_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Fe_over_H_time.txt', 'r')
    Fe_over_H_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Fe_over_H_time[1].split()]
    gas_Fe_over_H_list = [float(x) for x in Fe_over_H_time[3].split()]
    stellar_Fe_over_H_list = [float(x) for x in Fe_over_H_time[5].split()]

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Fe_over_H_time.txt', 'r')
    Fe_over_H_time_fixedIMF = file.readlines()
    file.close()
    log_time_axis_fixedIMF = [float(x) for x in Fe_over_H_time_fixedIMF[1].split()]
    gas_Fe_over_H_list_fixedIMF = [float(x) for x in Fe_over_H_time_fixedIMF[3].split()]
    stellar_Fe_over_H_list_fixedIMF = [float(x) for x in Fe_over_H_time_fixedIMF[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(11, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], ls='dashed', label='Gas', color="k")
    plt.plot([], [], label='Stellar', color="k")
    plt.plot(log_time_axis, gas_Fe_over_H_list, color='tab:orange', ls='dashed', lw=2)
    plt.plot(log_time_axis, stellar_Fe_over_H_list, color='tab:orange', lw=2, label='IGIMF')
    plt.plot(log_time_axis_fixedIMF, gas_Fe_over_H_list_fixedIMF, color='tab:blue', ls='dashed', lw=2)
    plt.plot(log_time_axis_fixedIMF, stellar_Fe_over_H_list_fixedIMF, color='tab:blue', lw=2, label='Fixed IMF')
    plt.plot([7, 10.1], [0, 0], color='red', ls='dotted', label='solar')
    plt.xlabel(r'log$_{10}$(age [yr])')
    plt.ylabel('[Fe/H]')
    plt.xlim(7, 10.1)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Fe_over_H_time.pdf', dpi=250)

    # Figure 12 Mg_over_Fe_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Mg_over_Fe_time.txt', 'r')
    Mg_over_Fe_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Mg_over_Fe_time[1].split()]
    gas_Mg_over_Fe_list = [float(x) for x in Mg_over_Fe_time[3].split()]
    stellar_Mg_over_Fe_list = [float(x) for x in Mg_over_Fe_time[5].split()]

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Mg_over_Fe_time.txt', 'r')
    Mg_over_Fe_time_fixedIMF = file.readlines()
    file.close()
    log_time_axis_fixedIMF = [float(x) for x in Mg_over_Fe_time_fixedIMF[1].split()]
    gas_Mg_over_Fe_list_fixedIMF = [float(x) for x in Mg_over_Fe_time_fixedIMF[3].split()]
    stellar_Mg_over_Fe_list_fixedIMF = [float(x) for x in Mg_over_Fe_time_fixedIMF[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(12, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot(log_time_axis, gas_Mg_over_Fe_list, color='tab:orange', ls='dashed', lw=2)
    plt.plot(log_time_axis, stellar_Mg_over_Fe_list, color='tab:orange', lw=2, label='IGIMF')
    plt.plot(log_time_axis_fixedIMF, gas_Mg_over_Fe_list_fixedIMF, color='tab:blue', ls='dashed', lw=2)
    plt.plot(log_time_axis_fixedIMF, stellar_Mg_over_Fe_list_fixedIMF, color='tab:blue', lw=2, label='Fixed IMF')
    plt.plot([7, 10.1], [0, 0], color='red', ls='dotted', label='solar')
    plt.xlabel(r'log$_{10}$(age [yr])')
    plt.ylabel('[Mg/Fe]')
    plt.xlim(7, 10.1)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Mg_over_Fe_time.pdf', dpi=250)

    plt.show()

    return


if __name__ == '__main__':
    plot_galaxy_evol()
