import matplotlib.pyplot as plt
import math


def plot_galaxy_evol():
    # Figure 1 SFH

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/SFH.txt', 'r')
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
    fig = plt.figure(1, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot(age_list, SFR_list)
    # plt.plot(age_list, SFR_list, ls='dashed', label='Flat')
    # plt.plot(age_list0, SFR_list0, label='Skewed-normal')
    plt.xlabel('Time [Gyr]')
    plt.ylabel(r'log$_{10}({\rm SFR} [{\rm M}_\odot$/yr])')
    plt.legend(prop={'size': 7}, loc='best')
    plt.xlim(-0.3, 2.6)
    plt.tight_layout()
    plt.savefig('SFH.pdf', dpi=250)

    # Figure 6 1000_IMF_evolution

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/imf_at_time_10_Myr.txt', 'r')
    imf_at_time_ = file.readlines()
    file.close()
    mass_list = [math.log(float(x), 10) for x in imf_at_time_[1].split()]
    xi_Kroupa = [math.log(float(x)+1e-10, 10) for x in imf_at_time_[5].split()]
    xi_each_time = []
    i = 100
    while i > 0:
        file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/imf_at_time_{}0_Myr.txt'.format(i), 'r')
        imf_at_time_ = file.readlines()
        file.close()
        xi_each_time.append([math.log(float(x)+1e-10, 10) for x in imf_at_time_[3].split()])
        (i) = (i-1)
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(6, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot(mass_list, xi_Kroupa, lw=3, label='Canonical IMF', ls='dashed')
    i = 99
    while i > -1:
        plt.plot(mass_list, xi_each_time[i], lw=1)
        (i) = (i - 1)
    # plt.plot([], [])
    plt.xlabel(r'log$_{10}$(M$_*$ [M$_\odot$])')
    plt.ylabel(r'log$_{10}$($\xi$ [M$_\odot^{-1}$])')
    plt.legend(prop={'size': 7}, loc='best')
    # plt.xlim(6.7, 10.1)
    # plt.ylim(8.8, 13.3)
    plt.tight_layout()
    plt.savefig('IMF_evolution.pdf', dpi=250)
        

    # Figure 7 mass_evolution

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/mass_evolution.txt', 'r')
    mass_evolution = file.readlines()
    file.close()
    time_axis = [float(x) for x in mass_evolution[1].split()]
    total_gas_mass_list = [float(x) for x in mass_evolution[3].split()]
    ejected_gas_mass_list = [float(x) for x in mass_evolution[5].split()]
    stellar_mass_list = [float(x) for x in mass_evolution[7].split()]
    remnant_mass_list = [float(x) for x in mass_evolution[9].split()]
    BH_mass_list = [float(x) for x in mass_evolution[11].split()]
    NS_mass_list = [float(x) for x in mass_evolution[13].split()]
    WD_mass_list = [float(x) for x in mass_evolution[15].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(7, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot(time_axis, total_gas_mass_list, lw=3, label='all gas', ls='dashed')
    plt.plot(time_axis, ejected_gas_mass_list, lw=1, label='ejected gas', ls='dashed')
    # plt.plot([], [])
    plt.plot(time_axis, stellar_mass_list, lw=3, label='alive stars')
    plt.plot(time_axis, remnant_mass_list, lw=3, label='all remnants')
    plt.plot(time_axis, BH_mass_list, lw=1, label='black holes')
    plt.plot(time_axis, NS_mass_list, lw=1, label='neutron stars')
    plt.plot(time_axis, WD_mass_list, lw=1, label='white dwarfs')
    plt.xlabel(r'log$_{10}$(age [yr])')
    plt.ylabel(r'log$_{10}$(Mass [M$_\odot$])')
    plt.legend(prop={'size': 7}, loc='best')
    plt.xlim(6.7, 10.1)
    # plt.ylim(8.8, 13.3)
    plt.tight_layout()
    plt.savefig('mass_evolution.pdf', dpi=250)

    # Figure 8 SN_number_evolution

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/SN_number_evolution.txt', 'r')
    SN_number_evolution = file.readlines()
    file.close()
    time_axis = [float(x) for x in SN_number_evolution[1].split()]
    SNIa_number_per_century = [float(x) for x in SN_number_evolution[3].split()]
    SNII_number_per_century = [float(x) for x in SN_number_evolution[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(8, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.loglog(time_axis, SNIa_number_per_century, label='SNIa')  # Number per century
    plt.loglog(time_axis, SNII_number_per_century, label='SNII')  # Number per century
    plt.xlabel(r'Age [yr]')
    plt.ylabel(r'# of SN per century')
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()

    # # Figure 12 energy_evolution
    #
    # file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/energy_evolution.txt', 'r')
    # energy_evolution = file.readlines()
    # file.close()
    #
    # time_axis = [float(x) for x in energy_evolution[1].split()]
    # total_energy_release_list = [float(x) for x in energy_evolution[7].split()]
    # binding_energy_list = [float(x) for x in energy_evolution[9].split()]
    # total_gas_kinetic_energy_list = [float(x) for x in energy_evolution[11].split()]
    #
    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(12, figsize=(4, 3.5))
    # fig.add_subplot(1, 1, 1)
    # plt.loglog(time_axis, total_energy_release_list, lw=2, label='Total SN')
    # plt.loglog(time_axis, binding_energy_list, lw=2, ls="dashed", label='Gravitational binding')
    # plt.loglog(time_axis, total_gas_kinetic_energy_list, lw=2, ls="dotted", label='Gas kinetic')
    # plt.xlabel(r'Age [yr]')
    # plt.ylabel(r'Energy [erg]')
    # plt.legend(prop={'size': 7}, loc='best')
    # plt.tight_layout()

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

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Z_over_X_time.txt', 'r')
    Z_over_X_time = file.readlines()
    file.close()

    log_time_axis = [float(x) for x in Z_over_X_time[1].split()]
    gas_Z_over_X_list = [float(x) for x in Z_over_X_time[3].split()]
    stellar_Z_over_X_list = [float(x) for x in Z_over_X_time[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(9, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot(log_time_axis, gas_Z_over_X_list, label='gas', color='k', ls='dashed', lw=2)
    plt.plot(log_time_axis, stellar_Z_over_X_list, label='stellar', color='k', lw=2)
    plt.plot([6.5, 10], [0, 0], color='red', ls='dashed', label='solar')
    # The [Z/X]s where the applied portinari98 stellar yield table will be changed for Z=0.0127, 0.008, 0.004, 0.0004.
    plt.plot([6.5, 10], [-1.173, -1.173], color='red', lw=0.5)
    plt.plot([6.5, 10], [-0.523, -0.523], color='red', lw=0.5)
    plt.plot([6.5, 10], [-0.272, -0.272], color='red', lw=0.5)
    plt.xlabel(r'log$_{10}$(age [yr])')
    plt.ylabel('[Z/X]')
    plt.xlim(6.5, 10)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()

    # Figure 10 Y_time

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Y_time.txt', 'r')
    Y_time = file.readlines()
    file.close()

    log_time_axis = [float(x) for x in Y_time[1].split()]
    gas_Y_list = [float(x) for x in Y_time[3].split()]
    stellar_Y_list = [float(x) for x in Y_time[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(10, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot(log_time_axis, gas_Y_list, label='gas', color='k', ls='dashed', lw=2)
    plt.plot(log_time_axis, stellar_Y_list, label='stellar', color='k', lw=2)
    plt.plot([6.5, 10], [0.2485, 0.2485], color='red', ls='dashed', label='solar')
    plt.xlabel(r'log$_{10}$(age [yr])')
    plt.ylabel('Y')
    plt.xlim(6.5, 10)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()

    # Figure 11 Fe_over_H_time

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Fe_over_H_time.txt', 'r')
    Fe_over_H_time = file.readlines()
    file.close()

    log_time_axis = [float(x) for x in Fe_over_H_time[1].split()]
    gas_Fe_over_H_list = [float(x) for x in Fe_over_H_time[3].split()]
    stellar_Fe_over_H_list = [float(x) for x in Fe_over_H_time[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(11, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot(log_time_axis, gas_Fe_over_H_list, label='gas', color='k', ls='dashed', lw=2)
    plt.plot(log_time_axis, stellar_Fe_over_H_list, label='stellar', color='k', lw=2)
    plt.plot([6.5, 10], [0, 0], color='red', ls='dashed', label='solar')
    plt.xlabel(r'log$_{10}$(age [yr])')
    plt.ylabel('[Fe/H]')
    plt.xlim(6.5, 10)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()

    # Figure 12 Mg_over_Fe_time

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Mg_over_Fe_time.txt', 'r')
    Mg_over_Fe_time = file.readlines()
    file.close()

    log_time_axis = [float(x) for x in Mg_over_Fe_time[1].split()]
    gas_Mg_over_Fe_list = [float(x) for x in Mg_over_Fe_time[3].split()]
    stellar_Mg_over_Fe_list = [float(x) for x in Mg_over_Fe_time[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(12, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot(log_time_axis, gas_Mg_over_Fe_list, label='gas', color='k', ls='dashed', lw=2)
    plt.plot(log_time_axis, stellar_Mg_over_Fe_list, label='stellar', color='k', lw=2)
    plt.plot([6.5, 10], [0, 0], color='red', ls='dashed', label='solar')
    plt.xlabel(r'log$_{10}$(age [yr])')
    plt.ylabel('[Mg/Fe]')
    plt.xlim(6.5, 10)
    plt.ylim(-0.1, 1.7)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()

    plt.show()

    return


if __name__ == '__main__':
    plot_galaxy_evol()
