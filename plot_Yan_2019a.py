import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.stats import lognorm
from matplotlib.collections import LineCollection

def multiline(xs, ys, c, ax=None, **kwargs):
    """Plot lines with different colorings
    Code copy from https://stackoverflow.com/a/50029441

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    c : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """
    ax = plt.gca() if ax is None else ax
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = LineCollection(segments, **kwargs)
    lc.set_array(np.asarray(c))
    ax.add_collection(lc)
    ax.autoscale()
    return lc

def plot_galaxy_evol():

    # Figure 1 SFH

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/SFH.txt', 'r')
    SFH = file.readlines()
    file.close()
    age_list = [float(x) for x in SFH[1].split()]
    SFR_list = [float(x) for x in SFH[3].split()]
    # age_list = [-1] + age_list
    # SFR_list = [-10] + SFR_list

    file = open('simulation_results_from_galaxy_evol/plots_extended/SFH.txt', 'r')
    SFH = file.readlines()
    file.close()
    age_list_extended = [float(x) for x in SFH[1].split()]
    SFR_list_extended = [float(x) for x in SFH[3].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(1, figsize=(4, 2.5))
    fig.add_subplot(1, 1, 1)
    plt.plot(age_list, SFR_list, label='boxy', ls='dashed', color="tab:orange")
    plt.plot(age_list_extended, SFR_list_extended, color="tab:red", label='log-normal')
    plt.xlabel('time [Gyr]')
    plt.ylabel(r'log$_{10}({\rm SFR} [{\rm M}_\odot$/yr])')
    plt.legend(prop={'size': 7}, loc='best')
    plt.xlim(-0.5, 13)
    plt.ylim(-1, 3.5)
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/SFH.pdf', dpi=250)

    # Figure 4 1000_IMF_evolution

    # the metallicity (range from -6 to 1 dex) is used for color code
    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Z_over_X_time.txt', 'r')
    Z_over_X_time_ = file.readlines()
    file.close()
    stellar_Z_over_X_color_list_ = [float(x) for x in Z_over_X_time_[3].split()]

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/imf_at_time_10_Myr.txt', 'r')
    imf_at_time_ = file.readlines()
    file.close()
    mass_list = [math.log(float(x), 10) for x in imf_at_time_[1].split()]
    xi_Kroupa = [math.log(float(x)+1e-10, 10) for x in imf_at_time_[5].split()]
    xi_each_time = []
    time_length = 100

    i = 1
    while i < time_length + 1:
        file = open('simulation_results_from_galaxy_evol/plots_IGIMF/imf_at_time_{}0_Myr.txt'.format(i), 'r')
        imf_at_time_ = file.readlines()
        file.close()
        xi_each_time.append([math.log(float(x)+1e-10, 10) for x in imf_at_time_[3].split()])
        (i) = (i + 1)

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(4, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    i = time_length - 1

    xs = []
    ys = []
    color = []
    while i > -1:
        xs.append(mass_list)
        ys.append(xi_each_time[i])
        color.append(stellar_Z_over_X_color_list_[i])
        (i) = (i - 1)
    lc = multiline(xs, ys, color, cmap='rainbow', lw=1)
    axcb = fig.colorbar(lc)
    axcb.set_label('[Z]')
    plt.plot(mass_list, xi_Kroupa, c='k', label='Canonical IMF', ls='dashed')
    # plt.plot([], [])
    plt.xlabel(r'log$_{10}$(M$_*$ [M$_\odot$])')
    plt.ylabel(r'log$_{10}$($\xi$ [M$_\odot^{-1}$])')
    plt.legend(prop={'size': 7}, loc='best')
    plt.xlim(-1.2, 2.26)
    plt.ylim(3.5, 13)
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/IMF_evolution.pdf', dpi=250)


    # Figure 5 mass_evolution

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/mass_evolution.txt', 'r')
    mass_evolution = file.readlines()
    file.close()
    time_axis = [float(x) for x in mass_evolution[1].split()]
    # print(time_axis)
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

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(5, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)

    plt.plot([], [], color="k", label='alive stars')
    plt.plot([], [], ls='dashed', color="k", label='all remnants')
    plt.plot([], [], ls='dotted', color="k", lw=3, label='all gas')
    plt.plot([], [], ls='dotted', color="k", label='ejected gas')
    plt.plot([], [], label='Canonical IMF', color="tab:blue", lw=3)
    plt.plot([], [], label='IGIMF', color="tab:orange", lw=3)


    plt.plot(time_axis, total_gas_mass_list, lw=3, ls='dotted', color="tab:orange")
    plt.plot(time_axis, ejected_gas_mass_list, ls='dotted', color="tab:orange")
    plt.plot(time_axis, stellar_mass_list, color="tab:orange")
    plt.plot(time_axis, remnant_mass_list, ls='dashed', color="tab:orange")

    plt.plot(time_axis_fixedIMF, total_gas_mass_list_fixedIMF, lw=3, ls='dotted', color="tab:blue")
    plt.plot(time_axis_fixedIMF, ejected_gas_mass_list_fixedIMF, ls='dotted', color="tab:blue")
    plt.plot(time_axis_fixedIMF, stellar_mass_list_fixedIMF, color="tab:blue")
    plt.plot(time_axis_fixedIMF, remnant_mass_list_fixedIMF, ls='dashed', color="tab:blue")

    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel(r'log$_{10}$(Mass [M$_\odot$])')
    plt.legend(prop={'size': 7}, loc='best')
    plt.xlim(7, 10.114)
    plt.ylim(7.8, 12.3)
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/mass_evolution.pdf', dpi=250)

    # Figure 6 SN_number_evolution

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/SN_number_evolution.txt', 'r')
    SN_number_evolution = file.readlines()
    file.close()
    time_axis = [float(x) for x in SN_number_evolution[1].split()]
    SNIa_number_per_century = [float(x)+10**-10 for x in SN_number_evolution[3].split()]
    SNII_number_per_century = [float(x)+10**-10 for x in SN_number_evolution[5].split()]

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/SN_number_evolution.txt', 'r')
    SN_number_evolution_fixedIMF = file.readlines()
    file.close()
    time_axis_fixedIMF = [float(x) for x in SN_number_evolution_fixedIMF[1].split()]
    SNIa_number_per_century_fixedIMF = [float(x)+10**-10 for x in SN_number_evolution_fixedIMF[3].split()]
    SNII_number_per_century_fixedIMF = [float(x)+10**-10 for x in SN_number_evolution_fixedIMF[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(6, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], label='SNIa', color="k", ls='dotted')
    plt.plot([], [], label='SNII', color="k", lw=0.5)
    plt.plot([], [], label='Canonical IMF', color="tab:blue", lw=3)
    plt.plot([], [], label='IGIMF', color="tab:orange", lw=3)
    plt.loglog(time_axis_fixedIMF, SNIa_number_per_century_fixedIMF, color="tab:blue", ls='dotted')  # Number per century
    plt.loglog(time_axis_fixedIMF, SNII_number_per_century_fixedIMF, color="tab:blue", lw=0.5)  # Number per century
    plt.loglog(time_axis, SNIa_number_per_century, color="tab:orange", ls='dotted')  # Number per century
    plt.loglog(time_axis, SNII_number_per_century, color="tab:orange", lw=0.5)  # Number per century
    plt.xlim(10**7, 14*10**9)
    plt.ylim(0.5, 10**5.5)
    plt.xlabel(r'time [yr]')
    plt.ylabel(r'# of SN per century')
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/SN_number_evolution.pdf', dpi=250)

    # Figure 7 Z_over_X_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Z_over_X_time.txt', 'r')
    Z_over_X_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Z_over_X_time[1].split()]
    gas_Z_over_X_list = [float(x) for x in Z_over_X_time[3].split()]
    stellar_Z_over_X_list = [float(x) for x in Z_over_X_time[5].split()]
    stellar_Z_over_X_list_luminosity_weighted = [float(x) for x in Z_over_X_time[7].split()]

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Z_over_X_time.txt', 'r')
    Z_over_X_time_fixedIMF = file.readlines()
    file.close()
    log_time_axis_fixedIMF = [float(x) for x in Z_over_X_time_fixedIMF[1].split()]
    gas_Z_over_X_list_fixedIMF = [float(x) for x in Z_over_X_time_fixedIMF[3].split()]
    stellar_Z_over_X_list_fixedIMF = [float(x) for x in Z_over_X_time_fixedIMF[5].split()]
    stellar_Z_over_X_list_luminosity_weighted_fixedIMF = [float(x) for x in Z_over_X_time_fixedIMF[7].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(7, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], label='Gas', lw=0.7, color="k")
    plt.plot([], [], ls='dashed', lw=0.7, label='Stellar mass-weighted', color="k")
    plt.plot([], [], ls='dotted', lw=0.7, label='Stellar luminosity-weighted', color="k")
    plt.plot([], [], label='Canonical IMF', color="tab:blue", lw=3)
    plt.plot([], [], label='IGIMF', color="tab:orange", lw=3)
    plt.plot(log_time_axis, gas_Z_over_X_list, color="tab:orange", lw=0.7)
    plt.plot(log_time_axis, stellar_Z_over_X_list, color="tab:orange", lw=0.7, ls='dashed')
    plt.plot(log_time_axis, stellar_Z_over_X_list_luminosity_weighted, ls='dotted', lw=0.7, color="tab:orange")
    plt.plot(log_time_axis_fixedIMF, gas_Z_over_X_list_fixedIMF, color="tab:blue", lw=0.7)
    plt.plot(log_time_axis_fixedIMF, stellar_Z_over_X_list_fixedIMF, color="tab:blue", lw=0.7, ls='dashed')
    plt.plot(log_time_axis_fixedIMF, stellar_Z_over_X_list_luminosity_weighted_fixedIMF, ls='dotted', lw=0.7, color="tab:blue")
    plt.plot([7, 10.114], [0, 0], color='red', ls='dotted', label='Solar')
    plt.plot([7, 10.114], [-0.933, -0.933], color='k', lw=0.5, ls='-.', label='Yield Switch')
    plt.plot([7, 10.114], [-0.497, -0.497], color='k', lw=0.5, ls='-.')
    plt.plot([7, 10.114], [-0.26, -0.26], color='k', lw=0.5, ls='-.')
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('[Z/X]')
    plt.xlim(7, 10.114)
    plt.ylim(-3.1, 1.1)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Z_over_X_time.pdf', dpi=250)

    # Figure 8 Y_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Y_time.txt', 'r')
    Y_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Y_time[1].split()]
    gas_Y_list = [float(x) for x in Y_time[3].split()]
    stellar_Y_list = [float(x) for x in Y_time[5].split()]
    stellar_Y_list_luminosity_weighted = [float(x) for x in Y_time[7].split()]

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Y_time.txt', 'r')
    Y_time_fixedIMF = file.readlines()
    file.close()
    log_time_axis_fixedIMF = [float(x) for x in Y_time_fixedIMF[1].split()]
    gas_Y_list_fixedIMF = [float(x) for x in Y_time_fixedIMF[3].split()]
    stellar_Y_list_fixedIMF = [float(x) for x in Y_time_fixedIMF[5].split()]
    stellar_Y_list_luminosity_weighted_fixedIMF = [float(x) for x in Y_time_fixedIMF[7].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(8, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], label='Gas', lw=0.7, color="k")
    plt.plot([], [], ls='dashed', lw=0.7, label='Stellar mass-weighted', color="k")
    plt.plot([], [], ls='dotted', lw=0.7, label='Stellar luminosity-weighted', color="k")
    plt.plot([], [], label='Canonical IMF', color="tab:blue", lw=3)
    plt.plot([], [], label='IGIMF', color="tab:orange", lw=3)
    plt.plot(log_time_axis, gas_Y_list, color='tab:orange', lw=0.7)
    plt.plot(log_time_axis, stellar_Y_list, color='tab:orange', lw=0.7, ls='dashed')
    plt.plot(log_time_axis, stellar_Y_list_luminosity_weighted, ls='dotted', lw=0.7, color="tab:orange")
    plt.plot(log_time_axis_fixedIMF, gas_Y_list_fixedIMF, color="tab:blue", lw=0.7)
    plt.plot(log_time_axis_fixedIMF, stellar_Y_list_fixedIMF, color="tab:blue", lw=0.7, ls='dashed')
    plt.plot(log_time_axis_fixedIMF, stellar_Y_list_luminosity_weighted_fixedIMF, ls='dotted', lw=0.7, color="tab:blue")
    plt.plot([6.5, 10.2], [0.273, 0.273], color='red', ls='dotted', label='solar')
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('Y')
    plt.xlim(7, 10.114)
    plt.ylim(0.24, 0.38)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Y_time.pdf', dpi=250)

    # Figure 9 Fe_over_H_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Fe_over_H_time.txt', 'r')
    Fe_over_H_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Fe_over_H_time[1].split()]
    gas_Fe_over_H_list = [float(x) for x in Fe_over_H_time[3].split()]
    stellar_Fe_over_H_list = [float(x) for x in Fe_over_H_time[5].split()]
    stellar_Fe_over_H_list_luminosity_weighted = [float(x) for x in Fe_over_H_time[7].split()]

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Fe_over_H_time.txt', 'r')
    Fe_over_H_time_fixedIMF = file.readlines()
    file.close()
    log_time_axis_fixedIMF = [float(x) for x in Fe_over_H_time_fixedIMF[1].split()]
    gas_Fe_over_H_list_fixedIMF = [float(x) for x in Fe_over_H_time_fixedIMF[3].split()]
    stellar_Fe_over_H_list_fixedIMF = [float(x) for x in Fe_over_H_time_fixedIMF[5].split()]
    stellar_Fe_over_H_list_luminosity_weighted_fixedIMF = [float(x) for x in Fe_over_H_time_fixedIMF[7].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(9, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], label='Gas', lw=0.7, color="k")
    plt.plot([], [], ls='dashed', lw=0.7, label='Stellar mass-weighted', color="k")
    plt.plot([], [], ls='dotted', lw=0.7, label='Stellar luminosity-weighted', color="k")
    plt.plot([], [], label='Canonical IMF', color="tab:blue", lw=3)
    plt.plot([], [], label='IGIMF', color="tab:orange", lw=3)
    plt.plot(log_time_axis, gas_Fe_over_H_list, color='tab:orange', lw=0.7)
    plt.plot(log_time_axis, stellar_Fe_over_H_list, color='tab:orange', lw=0.7, ls='dashed')
    plt.plot(log_time_axis, stellar_Fe_over_H_list_luminosity_weighted, ls='dotted', lw=0.7, color="tab:orange")
    plt.plot(log_time_axis_fixedIMF, gas_Fe_over_H_list_fixedIMF, color='tab:blue', lw=0.7)
    plt.plot(log_time_axis_fixedIMF, stellar_Fe_over_H_list_fixedIMF, color='tab:blue', lw=0.7, ls='dashed')
    plt.plot(log_time_axis_fixedIMF, stellar_Fe_over_H_list_luminosity_weighted_fixedIMF, ls='dotted', lw=0.7, color="tab:blue")
    plt.plot([7, 10.114], [0, 0], color='red', ls='dotted', label='solar')
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('[Fe/H]')
    plt.xlim(7, 10.114)
    plt.ylim(-4.5, 0.5)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Fe_over_H_time.pdf', dpi=250)

    # Figure 10 Mg_over_Fe_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Mg_over_Fe_time.txt', 'r')
    Mg_over_Fe_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Mg_over_Fe_time[1].split()]
    gas_Mg_over_Fe_list = [float(x) for x in Mg_over_Fe_time[3].split()]
    stellar_Mg_over_Fe_list = [float(x) for x in Mg_over_Fe_time[5].split()]
    stellar_Mg_over_Fe_list_luminosity_weighted = [float(x) for x in Mg_over_Fe_time[7].split()]

    file = open('simulation_results_from_galaxy_evol/plots_fixedIMF/Mg_over_Fe_time.txt', 'r')
    Mg_over_Fe_time_fixedIMF = file.readlines()
    file.close()
    log_time_axis_fixedIMF = [float(x) for x in Mg_over_Fe_time_fixedIMF[1].split()]
    gas_Mg_over_Fe_list_fixedIMF = [float(x) for x in Mg_over_Fe_time_fixedIMF[3].split()]
    stellar_Mg_over_Fe_list_fixedIMF = [float(x) for x in Mg_over_Fe_time_fixedIMF[5].split()]
    stellar_Mg_over_Fe_list_luminosity_weighted_fixedIMF = [float(x) for x in Mg_over_Fe_time_fixedIMF[7].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(10, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], label='Gas', lw=0.7, color="k")
    plt.plot([], [], ls='dashed', lw=0.7, label='Stellar mass-weighted', color="k")
    plt.plot([], [], ls='dotted', lw=0.7, label='Stellar luminosity-weighted', color="k")
    plt.plot([], [], label='Canonical IMF', color="tab:blue", lw=3)
    plt.plot([], [], label='IGIMF', color="tab:orange", lw=3)
    plt.plot(log_time_axis, gas_Mg_over_Fe_list, color='tab:orange', lw=0.7)
    plt.plot(log_time_axis, stellar_Mg_over_Fe_list, color='tab:orange', ls='dashed', lw=0.7)
    plt.plot(log_time_axis, stellar_Mg_over_Fe_list_luminosity_weighted, ls='dotted', lw=0.7, color="tab:orange")
    plt.plot(log_time_axis_fixedIMF, gas_Mg_over_Fe_list_fixedIMF, color='tab:blue', lw=0.7)
    plt.plot(log_time_axis_fixedIMF, stellar_Mg_over_Fe_list_fixedIMF, color='tab:blue', lw=0.7, ls='dashed')
    plt.plot(log_time_axis_fixedIMF, stellar_Mg_over_Fe_list_luminosity_weighted_fixedIMF, ls='dotted', lw=0.7, color="tab:blue")
    plt.plot([7, 10.114], [0, 0], color='red', ls='dotted', label='solar')
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('[Mg/Fe]')
    plt.xlim(7, 10.114)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Mg_over_Fe_time.pdf', dpi=250)

    return


def plot_galaxy_evol_extended():

    # Figure A1 1000_IMF_evolution

    # the metallicity (range from -6 to 1 dex) is used for color code
    file = open('simulation_results_from_galaxy_evol/plots_extended/Z_over_X_time.txt', 'r')
    Z_over_X_time_extended = file.readlines()
    file.close()
    stellar_Z_over_X_color_list_extended = [float(x) for x in Z_over_X_time_extended[3].split()]

    file = open('simulation_results_from_galaxy_evol/plots_extended/imf_at_time_10_Myr.txt', 'r')
    imf_at_time_ = file.readlines()
    file.close()
    mass_list = [math.log(float(x), 10) for x in imf_at_time_[1].split()]
    xi_Kroupa = [math.log(float(x) + 1e-10, 10) for x in imf_at_time_[5].split()]
    xi_each_time = []
    time_length = 1299

    i = 1
    while i < time_length+1:
        file = open('simulation_results_from_galaxy_evol/plots_extended/imf_at_time_{}0_Myr.txt'.format(i), 'r')
        imf_at_time_ = file.readlines()
        file.close()
        xi_each_time.append([math.log(float(x)+1e-10, 10) for x in imf_at_time_[3].split()])
        (i) = (i+1)


    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(101, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    i = time_length-1
    xs = []
    ys = []
    color = []
    while i > -1:
        xs.append(mass_list)
        ys.append(xi_each_time[i])
        color.append(stellar_Z_over_X_color_list_extended[i])
        (i) = (i - 1)
    lc = multiline(xs, ys, color, cmap='rainbow', lw=1)
    axcb = fig.colorbar(lc)
    axcb.set_label('[Z]')
    plt.plot(mass_list, xi_Kroupa, c='k', label='Canonical IMF', ls='dashed')
    # plt.plot([], [])
    plt.xlabel(r'log$_{10}$(M$_*$ [M$_\odot$])')
    plt.ylabel(r'log$_{10}$($\xi$ [M$_\odot^{-1}$])')
    plt.legend(prop={'size': 7}, loc='best')
    plt.xlim(-1.2, 2.26)
    plt.ylim(3.5, 13)
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/IMF_evolution_extended.pdf', dpi=250)


    # Figure A2 mass_evolution

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/mass_evolution.txt', 'r')
    mass_evolution = file.readlines()
    file.close()
    time_axis = [float(x) for x in mass_evolution[1].split()]
    total_gas_mass_list = [float(x) for x in mass_evolution[3].split()]
    ejected_gas_mass_list = [float(x) for x in mass_evolution[5].split()]
    stellar_mass_list = [float(x) for x in mass_evolution[7].split()]
    remnant_mass_list = [float(x) for x in mass_evolution[9].split()]

    file = open('simulation_results_from_galaxy_evol/plots_extended/mass_evolution.txt', 'r')
    mass_evolution_extended = file.readlines()
    file.close()
    time_axis_extended = [float(x) for x in mass_evolution_extended[1].split()]
    total_gas_mass_list_extended = [float(x) for x in mass_evolution_extended[3].split()]
    ejected_gas_mass_list_extended = [float(x) for x in mass_evolution_extended[5].split()]
    stellar_mass_list_extended = [float(x) for x in mass_evolution_extended[7].split()]
    remnant_mass_list_extended = [float(x) for x in mass_evolution_extended[9].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(102, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)

    plt.plot([], [], color="k", label='alive stars')
    plt.plot([], [], ls='dashed', color="k", label='all remnants')
    plt.plot([], [], ls='dotted', color="k", lw=3, label='all gas')
    plt.plot([], [], ls='dotted', color="k", label='ejected gas')
    plt.plot([], [], label='boxy SFH', lw=3, color="tab:orange")
    plt.plot([], [], label='log-norm SFH', lw=3, color="tab:red")

    plt.plot(time_axis, total_gas_mass_list, lw=3, ls='dotted', color="tab:orange")
    plt.plot(time_axis, ejected_gas_mass_list, ls='dotted', color="tab:orange")
    plt.plot(time_axis, stellar_mass_list, color="tab:orange")
    plt.plot(time_axis, remnant_mass_list, ls='dashed', color="tab:orange")

    plt.plot(time_axis_extended, total_gas_mass_list_extended, lw=3, ls='dotted', color="tab:red")
    plt.plot(time_axis_extended, ejected_gas_mass_list_extended, ls='dotted', color="tab:red")
    plt.plot(time_axis_extended, stellar_mass_list_extended, color="tab:red")
    plt.plot(time_axis_extended, remnant_mass_list_extended, ls='dashed', color="tab:red")
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel(r'log$_{10}$(Mass [M$_\odot$])')
    plt.legend(prop={'size': 7}, loc='best')
    plt.xlim(7, 10.114)
    plt.ylim(6.5, 12.3)
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/mass_evolution_extended.pdf', dpi=250)

    # Figure A3 SN_number_evolution

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/SN_number_evolution.txt', 'r')
    SN_number_evolution = file.readlines()
    file.close()
    time_axis = [float(x) for x in SN_number_evolution[1].split()]
    SNIa_number_per_century = [float(x)+10**-10 for x in SN_number_evolution[3].split()]
    SNII_number_per_century = [float(x)+10**-10 for x in SN_number_evolution[5].split()]

    file = open('simulation_results_from_galaxy_evol/plots_extended/SN_number_evolution.txt', 'r')
    SN_number_evolution_extended = file.readlines()
    file.close()
    time_axis_extended = [float(x) for x in SN_number_evolution_extended[1].split()]
    SNIa_number_per_century_extended = [float(x)+10**-10 for x in SN_number_evolution_extended[3].split()]
    SNII_number_per_century_extended = [float(x)+10**-10 for x in SN_number_evolution_extended[5].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(103, figsize=(4, 3.5))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], label='SNIa', color="k", ls='dotted')
    plt.plot([], [], label='SNII', color="k", lw=0.5)
    plt.plot([], [], label='boxy SFH', lw=3, color="tab:orange")
    plt.plot([], [], label='log-norm SFH', lw=3, color="tab:red")
    plt.loglog(time_axis_extended, SNIa_number_per_century_extended, color="tab:red", ls='dotted')  # Number per century
    plt.loglog(time_axis_extended, SNII_number_per_century_extended, color="tab:red", lw=0.5)  # Number per century
    plt.loglog(time_axis, SNIa_number_per_century, color="tab:orange", ls='dotted')  # Number per century
    plt.loglog(time_axis, SNII_number_per_century, color="tab:orange", lw=0.5)  # Number per century
    plt.xlim(10**7, 14*10**9)
    plt.ylim(0.5, 10**5.5)
    plt.xlabel(r'time [yr]')
    plt.ylabel(r'# of SN per century')
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/SN_number_evolution_extended.pdf', dpi=250)

    # Figure A4 Z_over_X_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Z_over_X_time.txt', 'r')
    Z_over_X_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Z_over_X_time[1].split()]
    gas_Z_over_X_list = [float(x) for x in Z_over_X_time[3].split()]
    stellar_Z_over_X_list = [float(x) for x in Z_over_X_time[5].split()]
    stellar_Z_over_X_list_luminosity_weighted = [float(x) for x in Z_over_X_time[7].split()]

    file = open('simulation_results_from_galaxy_evol/plots_extended/Z_over_X_time.txt', 'r')
    Z_over_X_time_extended = file.readlines()
    file.close()
    log_time_axis_extended = [float(x) for x in Z_over_X_time_extended[1].split()]
    gas_Z_over_X_list_extended = [float(x) for x in Z_over_X_time_extended[3].split()]
    stellar_Z_over_X_list_extended = [float(x) for x in Z_over_X_time_extended[5].split()]
    stellar_Z_over_X_list_luminosity_weighted_extended = [float(x) for x in Z_over_X_time_extended[7].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(104, figsize=(4, 4.8))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], label='Gas', lw=0.7, color="k")
    plt.plot([], [], ls='dashed', lw=0.7, label='Stellar mass-weighted', color="k")
    plt.plot([], [], ls='dotted', lw=0.7, label='Stellar luminosity-weighted', color="k")
    plt.plot([], [], label='boxy SFH', lw=3, color="tab:orange")
    plt.plot([], [], label='log-norm SFH', lw=3, color="tab:red")
    plt.plot(log_time_axis, gas_Z_over_X_list, color="tab:orange", lw=0.7)
    plt.plot(log_time_axis, stellar_Z_over_X_list, color="tab:orange", ls='dashed', lw=0.7)
    plt.plot(log_time_axis, stellar_Z_over_X_list_luminosity_weighted, ls='dotted', lw=0.7, color="tab:orange")
    plt.plot(log_time_axis_extended, gas_Z_over_X_list_extended, color="tab:red", lw=0.7)
    plt.plot(log_time_axis_extended, stellar_Z_over_X_list_extended, color="tab:red", ls='dashed', lw=0.7)
    plt.plot(log_time_axis_extended, stellar_Z_over_X_list_luminosity_weighted_extended, ls='dotted', lw=0.7, color="tab:red")
    plt.plot([7, 10.114], [0, 0], color='red', ls='dotted', label='Solar')
    # The [Z/X]s where the applied portinari98 stellar yield table will be changed for Z=0.0127, 0.008, 0.004, 0.0004.
    plt.plot([7, 10.114], [-0.933, -0.933], color='k', lw=0.5, ls='-.', label='Yield Switch')
    plt.plot([7, 10.114], [-0.497, -0.497], color='k', lw=0.5, ls='-.')
    plt.plot([7, 10.114], [-0.26, -0.26], color='k', lw=0.5, ls='-.')
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('[Z/X]')
    plt.xlim(7, 10.114)
    plt.ylim(-3.1, 1.1)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Z_over_X_time_extended.pdf', dpi=250)

    # Figure A5 Y_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Y_time.txt', 'r')
    Y_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Y_time[1].split()]
    gas_Y_list = [float(x) for x in Y_time[3].split()]
    stellar_Y_list = [float(x) for x in Y_time[5].split()]
    stellar_Y_list_luminosity_weighted = [float(x) for x in Y_time[7].split()]

    file = open('simulation_results_from_galaxy_evol/plots_extended/Y_time.txt', 'r')
    Y_time_extended = file.readlines()
    file.close()
    log_time_axis_extended = [float(x) for x in Y_time_extended[1].split()]
    gas_Y_list_extended = [float(x) for x in Y_time_extended[3].split()]
    stellar_Y_list_extended = [float(x) for x in Y_time_extended[5].split()]
    stellar_Y_list_luminosity_weighted_extended = [float(x) for x in Y_time_extended[7].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(105, figsize=(4, 4.8))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], label='Gas', lw=0.7, color="k")
    plt.plot([], [], ls='dashed', lw=0.7, label='Stellar mass-weighted', color="k")
    plt.plot([], [], ls='dotted', lw=0.7, label='Stellar luminosity-weighted', color="k")
    plt.plot([], [], label='boxy SFH', lw=3, color="tab:orange")
    plt.plot([], [], label='log-norm SFH', lw=3, color="tab:red")
    plt.plot(log_time_axis, gas_Y_list, color='tab:orange', lw=0.7)
    plt.plot(log_time_axis, stellar_Y_list, color='tab:orange', lw=0.7, ls='dashed')
    plt.plot(log_time_axis, stellar_Y_list_luminosity_weighted, ls='dotted', lw=0.7, color="tab:orange")
    plt.plot(log_time_axis_extended, gas_Y_list_extended, color="tab:red", lw=0.7)
    plt.plot(log_time_axis_extended, stellar_Y_list_extended, color="tab:red", lw=0.7, ls='dashed')
    plt.plot(log_time_axis_extended, stellar_Y_list_luminosity_weighted_extended, ls='dotted', lw=0.7, color="tab:red")
    plt.plot([6.5, 10.2], [0.273, 0.273], color='red', ls='dotted', label='solar')
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('Y')
    plt.xlim(7, 10.114)
    plt.ylim(0.24, 0.38)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Y_time_extended.pdf', dpi=250)

    # Figure A6 Fe_over_H_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Fe_over_H_time.txt', 'r')
    Fe_over_H_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Fe_over_H_time[1].split()]
    gas_Fe_over_H_list = [float(x) for x in Fe_over_H_time[3].split()]
    stellar_Fe_over_H_list = [float(x) for x in Fe_over_H_time[5].split()]
    stellar_Fe_over_H_list_luminosity_weighted = [float(x) for x in Fe_over_H_time[7].split()]

    file = open('simulation_results_from_galaxy_evol/plots_extended/Fe_over_H_time.txt', 'r')
    Fe_over_H_time_extended = file.readlines()
    file.close()
    log_time_axis_extended = [float(x) for x in Fe_over_H_time_extended[1].split()]
    gas_Fe_over_H_list_extended = [float(x) for x in Fe_over_H_time_extended[3].split()]
    stellar_Fe_over_H_list_extended = [float(x) for x in Fe_over_H_time_extended[5].split()]
    stellar_Fe_over_H_list_luminosity_weighted_extended = [float(x) for x in Fe_over_H_time_extended[7].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(106, figsize=(4, 4.8))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], label='Gas', lw=0.7, color="k")
    plt.plot([], [], ls='dashed', lw=0.7, label='Stellar mass-weighted', color="k")
    plt.plot([], [], ls='dotted', lw=0.7, label='Stellar luminosity-weighted', color="k")
    plt.plot([], [], label='boxy SFH', lw=3, color="tab:orange")
    plt.plot([], [], label='log-norm SFH', lw=3, color="tab:red")
    plt.plot(log_time_axis, gas_Fe_over_H_list, color='tab:orange', lw=0.7)
    plt.plot(log_time_axis, stellar_Fe_over_H_list, color='tab:orange', lw=0.7, ls='dashed')
    plt.plot(log_time_axis, stellar_Fe_over_H_list_luminosity_weighted, ls='dotted', lw=0.7, color="tab:orange")
    plt.plot(log_time_axis_extended, gas_Fe_over_H_list_extended, color='tab:red', lw=0.7)
    plt.plot(log_time_axis_extended, stellar_Fe_over_H_list_extended, color='tab:red', lw=0.7, ls='dashed')
    plt.plot(log_time_axis_extended, stellar_Fe_over_H_list_luminosity_weighted_extended, ls='dotted', lw=0.7, color="tab:red")
    plt.plot([7, 10.114], [0, 0], color='red', ls='dotted', label='solar')
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('[Fe/H]')
    plt.xlim(7, 10.114)
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Fe_over_H_time_extended.pdf', dpi=250)

    # Figure A7 Mg_over_Fe_time

    file = open('simulation_results_from_galaxy_evol/plots_IGIMF/Mg_over_Fe_time.txt', 'r')
    Mg_over_Fe_time = file.readlines()
    file.close()
    log_time_axis = [float(x) for x in Mg_over_Fe_time[1].split()]
    gas_Mg_over_Fe_list = [float(x) for x in Mg_over_Fe_time[3].split()]
    stellar_Mg_over_Fe_list = [float(x) for x in Mg_over_Fe_time[5].split()]
    stellar_Mg_over_Fe_list_luminosity_weighted = [float(x) for x in Mg_over_Fe_time[7].split()]

    file = open('simulation_results_from_galaxy_evol/plots_extended/Mg_over_Fe_time.txt', 'r')
    Mg_over_Fe_time_extended = file.readlines()
    file.close()
    log_time_axis_extended = [float(x) for x in Mg_over_Fe_time_extended[1].split()]
    gas_Mg_over_Fe_list_extended = [float(x) for x in Mg_over_Fe_time_extended[3].split()]
    stellar_Mg_over_Fe_list_extended = [float(x) for x in Mg_over_Fe_time_extended[5].split()]
    stellar_Mg_over_Fe_list_luminosity_weighted_extended = [float(x) for x in Mg_over_Fe_time_extended[7].split()]

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(107, figsize=(4, 4.8))
    fig.add_subplot(1, 1, 1)
    plt.plot([], [], label='Gas', lw=0.7, color="k")
    plt.plot([], [], ls='dashed', lw=0.7, label='Stellar mass-weighted', color="k")
    plt.plot([], [], ls='dotted', lw=0.7, label='Stellar luminosity-weighted', color="k")
    plt.plot([], [], label='boxy SFH', lw=3, color="tab:orange")
    plt.plot([], [], label='log-norm SFH', lw=3, color="tab:red")
    plt.plot(log_time_axis, gas_Mg_over_Fe_list, color='tab:orange', lw=0.7)
    plt.plot(log_time_axis, stellar_Mg_over_Fe_list, color='tab:orange', lw=0.7, ls='dashed')
    plt.plot(log_time_axis, stellar_Mg_over_Fe_list_luminosity_weighted, ls='dotted', lw=0.7, color="tab:orange")
    plt.plot(log_time_axis_extended, gas_Mg_over_Fe_list_extended, color='tab:red', lw=0.7)
    plt.plot(log_time_axis_extended, stellar_Mg_over_Fe_list_extended, color='tab:red', lw=0.7, ls='dashed')
    plt.plot(log_time_axis_extended, stellar_Mg_over_Fe_list_luminosity_weighted_extended, ls='dotted', lw=0.7, color="tab:red")
    plt.plot([7, 10.114], [0, 0], color='red', ls='dotted', label='solar')
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('[Mg/Fe]')
    plt.xlim(7, 10.114)
    plt.yticks(np.arange(0, 2.21, step=0.2))
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    plt.savefig('simulation_results_from_galaxy_evol/Mg_over_Fe_time_extended.pdf', dpi=250)

    return

if __name__ == '__main__':
    # plots for two different set of simulations, i.e., boxy SFH and lognormal SFH:
    plot_galaxy_evol()
    # plot_galaxy_evol_extended()
    plt.show()
