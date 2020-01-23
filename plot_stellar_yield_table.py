import time
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import element_abundances_solar

reference_name = 'Anders1989'
H_abundances_solar = element_abundances_solar.function_solar_element_abundances(reference_name, 'H')
# He_abundances_solar = element_abundances_solar.function_solar_element_abundances(reference_name, 'He')
C_abundances_solar = element_abundances_solar.function_solar_element_abundances(reference_name, 'C')
# N_abundances_solar = element_abundances_solar.function_solar_element_abundances(reference_name, 'N')
O_abundances_solar = element_abundances_solar.function_solar_element_abundances(reference_name, 'O')
Mg_abundances_solar = element_abundances_solar.function_solar_element_abundances(reference_name, 'Mg')
Fe_abundances_solar = element_abundances_solar.function_solar_element_abundances(reference_name, 'Fe')
Si_abundances_solar = element_abundances_solar.function_solar_element_abundances(reference_name, 'Si')
Ca_abundances_solar = element_abundances_solar.function_solar_element_abundances(reference_name, 'Ca')


def plot_lifetime_and_finalmass():
    Z2_list = [0.0004, 0.004, 0.008, 0.012]

    file = open('yield_tables/rearranged/setllar_final_mass_from_portinari98/portinari98_Z=0.004.txt', 'r')
    data = file.readlines()
    file.close()
    list2 = str.split(data[3])
    list_ini_mass = []
    for j in list2:
        list_ini_mass.append(math.log(float(j), 10))

    list_fin_mass = []
    i = len(Z2_list) - 1
    while i > -1:
        file = open('yield_tables/rearranged/setllar_final_mass_from_portinari98/portinari98_Z={}.txt'.format(Z2_list[i]), 'r')
        data = file.readlines()
        file.close()
        list2 = str.split(data[5])
        list3 = []
        for j in list2:
            list3.append(math.log(float(j), 10))
        list = [float(data[1]), list3]
        list_fin_mass.append(list)
        (i) = (i - 1)

    color_list_ = []
    for i in range(len(list_fin_mass)):
        ZZZ = list_fin_mass[i][0]
        Z_box = math.log(ZZZ, 10) - math.log(0.01886, 10)
        color_list_.append(round(((Z_box+7)**4.001 - (-6.001 + 7) ** 4.001) / ((1 + 7) ** 4.001 - (-6.001 + 7) ** 4.001) * 1000))

    colors = plt.cm.hsv_r(np.linspace(0, 1, 1000))

    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(21, figsize=(4, 3.5))
    # plt.xlim(-1.5, 2.5)
    # plt.ylim(-1.5, 1.5)
    # i = len(Z2_list) - 1
    # while i > -1:
    #     plt.plot(list_ini_mass, list_fin_mass[i][1], label='Z={}'.format(list_fin_mass[i][0]))
    #     (i) = (i - 1)
    # plt.plot([-2, 3], [-2, 3], ls='dashed', c='k', lw=0.7)
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$ [$M_\odot$])')
    # plt.ylabel(r'log$_{10}$($M_{\rm *, final}$ [$M_\odot$])')
    # plt.tight_layout()
    # plt.savefig('Interpolated_stellar_final_mass.pdf', dpi=250)

    list_lifetime = []
    i = len(Z2_list) - 1
    while i > -1:
        file = open(
            'yield_tables/rearranged/setllar_lifetime_from_portinari98/portinari98_Z={}.txt'.format(Z2_list[i]), 'r')
        data = file.readlines()
        file.close()
        list2 = str.split(data[5])
        list3 = []
        for j in list2:
            list3.append(math.log(float(j), 10))
        list = [float(data[1]), list3]
        list_lifetime.append(list)
        (i) = (i - 1)

    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(22, figsize=(4, 3.5))
    # plt.xlim(-1.5, 2.5)
    # plt.ylim(6, 15)
    # i = len(Z2_list) - 1
    # while i > -1:
    #     plt.plot(list_ini_mass, list_lifetime[i][1], label='Z={}'.format(list_fin_mass[i][0]))
    #     (i) = (i - 1)
    # # plt.plot([-2, 3], [-2, 3], ls='dashed', c='k', lw=0.7)
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$ [$M_\odot$])')
    # plt.ylabel(r'log$_{10}$(life time [yr])')
    # plt.tight_layout()
    # plt.savefig('Interpolated_stellar_lifetime.pdf', dpi=250)

    ##########
    Metallicity_origen = [0.008, 0.02]
    Age_origen = [
        [6.47E+10, 3.54E+10, 2.09E+10, 1.30E+10, 8.46E+09, 5.72E+09, 4.12E+09, 2.92E+09, 2.36E+09, 2.18E+09, 1.82E+09,
         1.58E+09, 1.41E+09, 1.25E+09, 1.23E+09, 6.86E+08, 4.12E+08, 1.93E+08, 1.15E+08, 7.71E+07, 5.59E+07, 3.44E+07,
         2.10E+07, 1.49E+07, 1.01E+07, 6.65E+06, 5.30E+06, 4.15E+06, 3.44E+06, 3.32E+06],
        [7.92E+10, 4.45E+10, 2.61E+10, 1.59E+10, 1.03E+10, 6.89E+09, 4.73E+09, 3.59E+09, 2.87E+09, 2.64E+09, 2.18E+09,
         1.84E+09, 1.59E+09, 1.38E+09, 1.21E+09, 7.64E+08, 4.56E+08, 2.03E+08, 1.15E+08, 7.45E+07, 5.31E+07, 3.17E+07,
         1.89E+07, 1.33E+07, 9.15E+06, 6.13E+06, 5.12E+06, 4.12E+06, 3.39E+06, 3.23E+06]]
    Age_012 = []
    for i in range(len(Age_origen[0])):
        Age_012.append((Age_origen[0][i]*2+Age_origen[1][i])/3)

    Remnant_mass_origen = [
        [1.35, 1.48, 1.84, 2.04, 6.9, 12.5, 5.69, 9.89],
        [1.31, 1.44, 1.87, 2.11, 7.18, 2.06, 2.09, 2.11]
    ]
    Remnant_mass_012 = []
    for i in range(len(Remnant_mass_origen[0])):
        Remnant_mass_012.append((Remnant_mass_origen[0][i]*2+Remnant_mass_origen[1][i])/3)


    Mass = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5,
            1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0,
            7.0, 9.0, 12., 15., 20., 30., 40., 60., 100, 120]
    Metallicity = [0.0004, 0.004, 0.008, 0.12]
    Age = [
        [4.28E+10, 2.37E+10, 1.41E+10, 8.97E+09, 6.03E+09, 4.23E+09, 3.08E+09, 2.34E+09, 1.92E+09, 1.66E+09, 1.39E+09,
         1.18E+09, 1.11E+09, 9.66E+08, 8.33E+08, 4.64E+08, 3.03E+08, 1.61E+08, 1.01E+08, 7.15E+07, 5.33E+07, 3.42E+07,
         2.13E+07, 1.54E+07, 1.06E+07, 6.90E+06, 5.45E+06, 4.20E+06, 3.32E+06, 3.11E+06],
        [5.35E+10, 2.95E+10, 1.73E+10, 1.09E+10, 7.13E+09, 4.93E+09, 3.52E+09, 2.64E+09, 2.39E+09, 1.95E+09, 1.63E+09,
         1.28E+09, 1.25E+09, 1.23E+09, 1.08E+09, 5.98E+08, 3.67E+08, 1.82E+08, 1.11E+08, 7.62E+07, 5.61E+07, 3.51E+07,
         2.14E+07, 1.52E+07, 1.05E+07, 6.85E+06, 5.44E+06, 4.19E+06, 3.38E+06, 3.23E+06],
        [6.47E+10, 3.54E+10, 2.09E+10, 1.30E+10, 8.46E+09, 5.72E+09, 4.12E+09, 2.92E+09, 2.36E+09, 2.18E+09, 1.82E+09,
         1.58E+09, 1.41E+09, 1.25E+09, 1.23E+09, 6.86E+08, 4.12E+08, 1.93E+08, 1.15E+08, 7.71E+07, 5.59E+07, 3.44E+07,
         2.10E+07, 1.49E+07, 1.01E+07, 6.65E+06, 5.30E+06, 4.15E+06, 3.44E+06, 3.32E+06],
        Age_012]

    len_mass = len(Mass)
    log_Mass = []
    for i in range(len_mass):
        log_Mass.append(math.log(Mass[i], 10))

    len_metal = len(Metallicity)
    log_Metallicity = []
    for i in range(len_metal):
        log_Metallicity.append(math.log(Metallicity[i], 10))

    log_Age = []
    for i in range(len_metal):
        log_Age.append([])
        for j in range(len_mass):
            log_Age[i].append(math.log(Age[i][j], 10))


    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4, 4))

    i = 0
    while i < len(Z2_list):
        ZZZ = list_fin_mass[i][0]
        Z_box = round(math.log(ZZZ, 10)-math.log(0.01886, 10), 2)
        axs[0].plot(list_ini_mass, list_lifetime[i][1], lw=(6-i)/2, label='Z={}, [Z]={}'.format(ZZZ, Z_box), color=colors[color_list_[i]])
        (i) = (i + 1)
    i = len_metal-1
    # while i > -1:
    #     axs[0].scatter(log_Mass, log_Age[i], s=3, marker='*', edgecolors='w', linewidth='0.1', zorder=10)
    #     (i) = (i - 1)
    axs[0].plot([-1, 2], [7, 7])
    axs[0].plot([math.log(17, 10), math.log(17, 10)], [6, 15])
    # axs[0].set_yticks(np.arange(6, 16, 2))
    axs[0].set_ylim(6, 15)
    axs[0].set_ylabel(r'log$_{10}$(life time [yr])')
    axs[0].legend(prop={'size': 6}, loc='best')

    Mass = [
        [9, 12, 15, 20, 30, 40, 60, 100, 120],
        [9, 12, 15, 20, 30, 40, 100, 120],
        [9, 12, 15, 20, 30, 40, 60, 120],
        [9, 12, 15, 20, 30, 40, 60, 120]
    ]
    Metallicity = [0.0004, 0.004, 0.008, 0.12]
    Remnant_mass = [
        [1.35, 1.5, 1.8, 2.07, 6.98, 14.91, 24.58, 32.06, 30.6],
        [1.35, 1.5, 1.82, 2.04, 6.98, 12.6, 36.7, 35.2],
        [1.35, 1.48, 1.84, 2.04, 6.9, 12.5, 5.69, 9.89],
        Remnant_mass_012
    ]

    #################################################################
    # WW95_solar = 0.01886
    # Metallicity_WW95 = [0, WW95_solar*10**-4, WW95_solar*0.01, WW95_solar*0.1, WW95_solar]
    # Mass_WW95 = [12, 13, 15, 18, 20, 22, 25, 30, 35, 40]
    # Remnant_mass_WW95_B = [
    #     [1.32, 1.46, 1.43, 1.76, 2.06, 2.02, 2.07, 1.94, 3.86, 5.45],
    #     [1.38, 1.31, 1.49, 1.69, 1.97, 2.12, 1.99, 2.01, 3.39, 4.45],
    #     [1.40, 1.44, 1.56, 1.58, 1.98, 2.04, 1.87, 2.21, 2.42, 4.42],
    #     [1.28, 1.44, 1.63, 1.61, 1.97, 2.01, 1.87, 2.08, 3.03, 4.09],
    #     [1.35, 1.28, 1.53, 3.40, 4.12, 1.49, 1.90, 1.54, 7.62, 12.2]
    # ]
    # Interpolation_remnant_mass_WW95_B = interpolate.interp2d(Mass_WW95, Metallicity_WW95, Remnant_mass_WW95_B)
    # Remnant_mass_WW95_B_new = []
    # for i in range(len(Metallicity)):
    #     Remnant_mass_WW95_B_new.append([])
    #     for j in range(len(Mass_WW95)):
    #         Remnant_mass_WW95_B_new[i].append(Interpolation_remnant_mass_WW95_B(Mass_WW95[j], Metallicity[i]))
    #
    # log_Remnant_mass_WW95_B = []
    # for i in range(len_metal):
    #     log_Remnant_mass_WW95_B.append([])
    #     for j in range(len(Remnant_mass_WW95_B[i])):
    #         log_Remnant_mass_WW95_B[i].append(math.log(Remnant_mass_WW95_B[i][j], 10))
    #
    # log_mass_WW95 = []
    # for i in range(len(Mass_WW95)):
    #     log_mass_WW95.append(math.log(Mass_WW95[i], 10))
    #################################################################

    len_metal = len(Metallicity)
    log_Metallicity = []
    for i in range(len_metal):
        log_Metallicity.append(math.log(Metallicity[i], 10))

    log_Remnant_mass = []
    for i in range(len_metal):
        log_Remnant_mass.append([])
        for j in range(len(Remnant_mass[i])):
            log_Remnant_mass[i].append(math.log(Remnant_mass[i][j], 10))

    log_mass = []
    for i in range(len_metal):
        log_mass.append([])
        for j in range(len(Mass[i])):
            log_mass[i].append(math.log(Mass[i][j], 10))

    # print(log_mass)
    # print(len(log_mass[0]))
    # print(len(log_mass))
    # print(len(log_Remnant_mass[0]))

    i = 0
    while i < len(Z2_list):
        axs[1].plot(list_ini_mass, list_fin_mass[i][1], lw=(6-i)/2, label='Z={}'.format(list_fin_mass[i][0]), color=colors[color_list_[i]])
        (i) = (i + 1)
    i = len_metal-1
    # while i > -1:
    #     axs[1].scatter(log_mass[i], log_Remnant_mass[i], s=10, marker='*', edgecolors='w', linewidth='0.1', zorder=10)
    #     (i) = (i - 1)
    # i = len_metal-1
    # # while i > -1:
    # #     axs[1].scatter(log_mass_WW95, log_Remnant_mass_WW95_B[i], s=10, marker='^', edgecolors='w', linewidth='0.1', zorder=10)
    # #     (i) = (i - 1)
    axs[1].set_yticks(np.arange(-2, 2, 1))
    axs[1].set_ylim(-1.5, 1.5)
    axs[1].set_ylabel(r'log$_{10}(M_{\rm *, final}$ [$M_\odot$])')
    axs[1].set_xlabel(r'log$_{10}(M_{\rm *, initial}$ [$M_\odot$])')

    plt.tight_layout()
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)
    plt.savefig('Interpolated_stellar_lifetime_final_mass.pdf', dpi=250)

    plt.show()

    return

def function_read_file(yield_table_name):

    ####################
    ### read in file ###
    ####################
    if yield_table_name == "portinari98":
        file_yield = open(
            'yield_tables/agb_and_massive_stars_portinari98_marigo01_gce_totalyields.txt', 'r')
            # 'yield_tables/agb_and_massive_stars_portinari98_marigo01.txt', 'r')
        # Use net yields of Portinari and Marigo
        # Net yields with masses up to 7Msun are from Marigo, above those of Portinari are taken.
        # Only isotopes are selected which are available in both yield sets and go up to Fe.
        # Initial masses go from the lowest mass available up to 100Msun.
        # Yield set ID M01P98 in Ritter et al. 2017.
        # References: Marigo et al. 2001, http://ukads.nottingham.ac.uk/abs/2001A%26A...370..194M
        #		      Portinari et al. 1998, http://ukads.nottingham.ac.uk/abs/1998A%26A...334..505P
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Kobayashi06":
        file_yield = open(
            'yield_tables/agb_and_massive_stars_Kobayashi06_marigo01_gce_totalyields.txt', 'r')
        # Use net yields of Woosley S. E., Weaver T. A., 1995, ApJS, 101, 181 (WW95)
        # Use WW95 model B which has the highest [Mg/Fe].
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "WW95":
        file_yield = open(
            'yield_tables/massive_stars_WW95_totalyields.txt', 'r')
        # Use net yields of Woosley S. E., Weaver T. A., 1995, ApJS, 101, 181 (WW95)
        # Use WW95 model B which has the highest [Mg/Fe].
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "marigo01":
        file_yield = open(
            'yield_tables/agb_marigo01_totalyields.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()


    ###########################
    ### extract information ###
    ###########################
    #
    H_relative_line_number = function_get_element_line_number(data, 'H-1')
    He_relative_line_number = function_get_element_line_number(data, 'He-4')
    C_relative_line_number = function_get_element_line_number(data, 'C-12')
    N_relative_line_number = function_get_element_line_number(data, 'N-14')
    O_relative_line_number = function_get_element_line_number(data, 'O-16')
    Ne_relative_line_number = function_get_element_line_number(data, 'Ne-20')
    Mg_relative_line_number = function_get_element_line_number(data, 'Mg-24')
    Si_relative_line_number = function_get_element_line_number(data, 'Si-28')
    S_relative_line_number = function_get_element_line_number(data, 'S-32')
    Ca_relative_line_number = function_get_element_line_number(data, 'Ca-40')
    Fe_relative_line_number = function_get_element_line_number(data, 'Fe-56')
    #
    global M_list, Z_list, eject_mass_list, H_eject_mass_list, He_eject_mass_list, C_eject_mass_list, \
        N_eject_mass_list, O_eject_mass_list, Ne_eject_mass_list, Mg_eject_mass_list, Si_eject_mass_list, \
        S_eject_mass_list, Ca_eject_mass_list,  Fe_eject_mass_list, Metal_eject_mass_list
    global O_over_Mg_list, Mg_over_Fe_list, Ca_over_Fe_list, Si_over_Fe_list, C_over_H_list, Mg_over_H_list, \
        Si_over_H_list, Fe_over_H_list, O_over_H_list, Z_over_H_list, \
        Z_over_X_list, Z_over_Z0_list, XXX_list, YYY_list, ZZZ_list, O_over_Fe_list
    #
    i = len(data)-1
    while i > -1:
        line_i = str.split(data[i])
        if line_i[1] == 'Table:':
            line_H = str.split(data[i + H_relative_line_number])
            line_He = str.split(data[i + He_relative_line_number])
            line_C = str.split(data[i + C_relative_line_number])
            line_N = str.split(data[i + N_relative_line_number])
            line_O = str.split(data[i + O_relative_line_number])
            line_Ne = str.split(data[i + Ne_relative_line_number])
            line_Mg = str.split(data[i + Mg_relative_line_number])
            line_Si = str.split(data[i + Si_relative_line_number])
            line_S = str.split(data[i + S_relative_line_number])
            line_Ca = str.split(data[i + Ca_relative_line_number])
            line_Fe = str.split(data[i + Fe_relative_line_number])
            line_Mfinal = str.split(data[i + 2])
            (Z, M) = function_get_Z_M(line_i[2])  # metallicity and mass of the star
            ejecta_mass = round((M - function_get_Mfinal(line_Mfinal[2])), 5)  ####################
            H_mass = function_get_element_mass(line_H[1])
            He_mass = function_get_element_mass(line_He[1])
            C_mass = function_get_element_mass(line_C[1])
            N_mass = function_get_element_mass(line_N[1])
            O_mass = function_get_element_mass(line_O[1])
            Ne_mass = function_get_element_mass(line_Ne[1])
            Mg_mass = function_get_element_mass(line_Mg[1])
            Si_mass = function_get_element_mass(line_Si[1])
            S_mass = function_get_element_mass(line_S[1])
            Ca_mass = function_get_element_mass(line_Ca[1])
            Fe_mass = function_get_element_mass(line_Fe[1])
            H_num = H_mass/1.0079
            C_num = C_mass/12.011
            N_num = N_mass/14.007
            O_num = O_mass/15.9994
            Ne_num = Ne_mass/20.18
            Mg_num = Mg_mass/24.305
            Si_num = Si_mass/28.085
            S_num = S_mass/32.06
            Ca_num = Ca_mass/40.078
            Fe_num = Fe_mass/55.845
            Metal_num = C_num+N_num+O_num+Ne_num+Mg_num+Si_num+S_num+Ca_num+Fe_num
            O_over_Mg = math.log(O_num/Mg_num, 10) - O_abundances_solar + Mg_abundances_solar
            Mg_over_H = math.log(Mg_num/H_num, 10) - Mg_abundances_solar + H_abundances_solar
            Si_over_H = math.log(Si_num/H_num, 10) - Si_abundances_solar + H_abundances_solar
            C_over_H = math.log(C_num/H_num, 10) - C_abundances_solar + H_abundances_solar
            Fe_over_H = math.log(Fe_num/H_num, 10) - Fe_abundances_solar + H_abundances_solar
            O_over_H = math.log(O_num/H_num, 10) - O_abundances_solar + H_abundances_solar
            Mg_over_Fe = math.log(Mg_num/Fe_num, 10) - Mg_abundances_solar + Fe_abundances_solar
            Ca_over_Fe = math.log(Ca_num/Fe_num, 10) - Ca_abundances_solar + Fe_abundances_solar
            Si_over_Fe = math.log(Si_num/Fe_num, 10) - Si_abundances_solar + Fe_abundances_solar
            O_over_Fe = math.log(O_num/Fe_num, 10) - O_abundances_solar + Fe_abundances_solar
            Metal_mass = round((ejecta_mass - H_mass - He_mass), 5)  ####################
            # Metal_mass = round((C_mass+N_mass+O_mass+Ne_mass+Mg_mass+Si_mass+S_mass+Ca_mass+Fe_mass), 5)  ###### the same ######
            if Metal_mass<0:
                print("Warning: Metal_mass=", Metal_mass, "<0")
                print("check stellar yield table with metallicity and mass being:", Z, "&", M)
                Metal_mass = 0
            Z_over_X = math.log(Metal_mass / H_mass, 10) - math.log(0.01886 / 0.7381, 10)
            Z_over_Z0 = math.log(Metal_mass / ejecta_mass, 10) - math.log(0.01886, 10)
            Z_over_H = math.log(Metal_num / H_num, 10) - math.log(0.01886 / 18 / 0.7381, 10)  # where 18 is the estimated average atomic weight over the weight of hydrogen.
            XXX = H_mass / ejecta_mass
            YYY = He_mass / ejecta_mass
            ZZZ = Metal_mass / ejecta_mass
            if len(Z_list) == 0:
                Z_list.append(Z)
                Z_n = 0
                M_list.append([])
                eject_mass_list.append([])
                H_eject_mass_list.append([])
                He_eject_mass_list.append([])
                C_eject_mass_list.append([])
                N_eject_mass_list.append([])
                O_eject_mass_list.append([])
                Ne_eject_mass_list.append([])
                Mg_eject_mass_list.append([])
                Si_eject_mass_list.append([])
                S_eject_mass_list.append([])
                Ca_eject_mass_list.append([])
                Fe_eject_mass_list.append([])
                Metal_eject_mass_list.append([])
                Z_over_H_list.append([])
                Z_over_X_list.append([])
                Z_over_Z0_list.append([])
                XXX_list.append([])
                YYY_list.append([])
                ZZZ_list.append([])
                O_over_Mg_list.append([])
                Mg_over_Fe_list.append([])
                Si_over_Fe_list.append([])
                Ca_over_Fe_list.append([])
                Mg_over_H_list.append([])
                Si_over_H_list.append([])
                C_over_H_list.append([])
                Fe_over_H_list.append([])
                O_over_H_list.append([])
                O_over_Fe_list.append([])
            if Z != Z_list[-1]:
                Z_list.append(Z)
                Z_n += 1
                M_list.append([])
                eject_mass_list.append([])
                H_eject_mass_list.append([])
                He_eject_mass_list.append([])
                C_eject_mass_list.append([])
                N_eject_mass_list.append([])
                O_eject_mass_list.append([])
                Ne_eject_mass_list.append([])
                Mg_eject_mass_list.append([])
                Si_eject_mass_list.append([])
                S_eject_mass_list.append([])
                Ca_eject_mass_list.append([])
                Fe_eject_mass_list.append([])
                Metal_eject_mass_list.append([])
                O_over_Mg_list.append([])
                Mg_over_Fe_list.append([])
                Ca_over_Fe_list.append([])
                Si_over_Fe_list.append([])
                Mg_over_H_list.append([])
                Si_over_H_list.append([])
                C_over_H_list.append([])
                Fe_over_H_list.append([])
                O_over_H_list.append([])
                Z_over_H_list.append([])
                Z_over_X_list.append([])
                Z_over_Z0_list.append([])
                XXX_list.append([])
                YYY_list.append([])
                ZZZ_list.append([])
                O_over_Fe_list.append([])
            M_list[Z_n].append(M)
            eject_mass_list[Z_n].append(ejecta_mass)
            H_eject_mass_list[Z_n].append(H_mass)
            He_eject_mass_list[Z_n].append(He_mass)
            C_eject_mass_list[Z_n].append(C_mass)
            N_eject_mass_list[Z_n].append(N_mass)
            O_eject_mass_list[Z_n].append(O_mass)
            Ne_eject_mass_list[Z_n].append(Ne_mass)
            Mg_eject_mass_list[Z_n].append(Mg_mass)
            Si_eject_mass_list[Z_n].append(Si_mass)
            S_eject_mass_list[Z_n].append(S_mass)
            Ca_eject_mass_list[Z_n].append(Ca_mass)
            Fe_eject_mass_list[Z_n].append(Fe_mass)
            Metal_eject_mass_list[Z_n].append(Metal_mass)
            O_over_Mg_list[Z_n].append(O_over_Mg)
            Mg_over_Fe_list[Z_n].append(Mg_over_Fe)
            Ca_over_Fe_list[Z_n].append(Ca_over_Fe)
            Si_over_Fe_list[Z_n].append(Si_over_Fe)
            Mg_over_H_list[Z_n].append(Mg_over_H)
            Si_over_H_list[Z_n].append(Si_over_H)
            C_over_H_list[Z_n].append(C_over_H)
            O_over_H_list[Z_n].append(O_over_H)
            Z_over_H_list[Z_n].append(Z_over_H)
            Z_over_X_list[Z_n].append(Z_over_X)
            Z_over_Z0_list[Z_n].append(Z_over_Z0)
            XXX_list[Z_n].append(XXX)
            YYY_list[Z_n].append(YYY)
            ZZZ_list[Z_n].append(ZZZ)
            Fe_over_H_list[Z_n].append(Fe_over_H)
            O_over_Fe_list[Z_n].append(O_over_Fe)
        (i) = (i - 1)

    return

def function_get_Mfinal(Mfinal_string):
    i_end = len(Mfinal_string)
    i = 0
    mass_str = ''
    while i < i_end:
        mass_str += Mfinal_string[i]
        (i) = (i + 1)
    mass = float(mass_str)
    return mass

def function_get_element_mass(element_mass_string):
    i_end = len(element_mass_string)
    i = 1
    mass_str = ''
    while i < i_end:
        mass_str += element_mass_string[i]
        (i) = (i + 1)
    mass = float(mass_str)
    return mass

def function_get_element_line_number(data, element):
    i = 0
    while i < len(data):
        line_i = str.split(data[i])
        if line_i[1] == 'Table:':
            start = i
            j = 0
            while j < 100:
                line_j = str.split(data[j])
                if line_j[0] == '&'+element:
                    end = j
                    element_relative_line_number = j - i
                    break
                (j) = (j+1)
            break
        (i) = (i + 1)
    return element_relative_line_number

def function_get_Z_M(M_Z_string):
    i = 0
    i_M_start = 0
    i_M_end = 0
    i_Z_start = 0
    i_Z_end = 0
    while i < len(M_Z_string):
        if M_Z_string[i] == 'M':
            i_M_start = i+2
        if M_Z_string[i] == ',':
            i_M_end = i
            i_Z_start = i+3
        if M_Z_string[i] == ')':
            i_Z_end = i
        (i) = (i+1)
    i = i_Z_start
    Z_str = ''
    while i < i_Z_end:
        Z_str += M_Z_string[i]
        (i) = (i + 1)
    Z = float(Z_str)
    i = i_M_start
    M_str = ''
    while i < i_M_end:
        M_str += M_Z_string[i]
        (i) = (i + 1)
    M = float(M_str)
    return (Z, M)

def funtion_plot_yields():
    global O_over_Mg_list, Mg_over_Fe_list, C_over_H_list, Mg_over_H_list, Si_over_H_list, Fe_over_H_list, O_over_H_list, Z_over_X_list, Z_over_Z0_list, \
        Z_over_H_list, O_over_Fe_list, M_list, Z_list, XXX_list, YYY_list, ZZZ_list
    color_list_ = []
    for i in range(len(Z_list)):
        ZZZ = Z_list[i]
        if ZZZ > 0:
            Z_box = math.log(ZZZ, 10) - math.log(0.01886, 10)
        else:
            Z_box = -6
        color_list_.append(round(((Z_box+7)**4.001 - (-6.001 + 7) ** 4.001) / ((1 + 7) ** 4.001 - (-6.001 + 7) ** 4.001) * 1000))

    colors = plt.cm.hsv_r(np.linspace(0, 1, 1000))


    j = 0
    while j < len(M_list):
        i = 0
        while i < len(M_list[j]):
            M_list[j][i] = math.log(M_list[j][i], 10)
            (i) = (i+1)
        (j) = (j+1)

    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(1, figsize=(4, 3.5))
    # plt.xlim(-0.5, 2.2)
    # # plt.ylim(0, 2)
    # i = len(M_list) - 1
    # while i > -1:
    #     plt.plot(M_list[i], O_over_Mg_list[i], label='Z={}'.format(Z_list[i]))
    #     (i) = (i - 1)
    # O_mass_eject_SNIa = 0.148  # TNH93 0.148 i99CDD1 0.09, i99CDD2 0.06, i99W7 0.14, ivo12/13 0.09-0.1, t03 0.14, t86 0.13
    # Mg_mass_eject_SNIa = 0.009  # TNH93 0.009 i99CDD1 0.0077, i99CDD2 0.0042, i99W7 0.0085, ivo12/13 0.015-0.029, t03 0.013, t86 0.016
    # O_num = O_mass_eject_SNIa / 15.9994
    # Mg_num = Mg_mass_eject_SNIa / 24.305
    # O_over_Mg_SNIa = math.log(O_num / Mg_num, 10) - O_abundances_solar + Mg_abundances_solar
    # plt.plot([-0.3, 0.9], [O_over_Mg_SNIa, O_over_Mg_SNIa], ls="--", lw=2, label="SNIa")
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$/[$M_\odot$])')
    # plt.ylabel(r'[O/Mg]')
    # plt.tight_layout()


    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(2, figsize=(4, 3.5))
    # i = len(M_list) - 1
    # while i > -1:
    #     plt.plot(M_list[i], Mg_over_Fe_list[i], label='Z={}'.format(Z_list[i]))
    #     (i) = (i - 1)
    Mg_mass_eject_SNIa = 0.0158  # TNH93 0.148 i99CDD1 0.09, i99CDD2 0.06, i99W7 0.14, ivo12/13 0.09-0.1, t03 0.14, t86 0.13
    Fe_mass_eject_SNIa = 0.68  #0.63 # Recchi2009 halfed to 0.372  # TNH93 0.744 i99CDD1 0.56, i99CDD2 0.76, i99W7 0.63, ivo12/13 0.62-0.67, t03 0.74, t86 0.63
    Ca_mass_eject_SNIa = 0.0181
    Si_mass_eject_SNIa = 0.142
    Ca_num = Ca_mass_eject_SNIa / 40.078
    Si_num = Si_mass_eject_SNIa / 28.085
    Mg_num = Mg_mass_eject_SNIa / 24.305
    Fe_num = Fe_mass_eject_SNIa / 55.845
    Mg_over_Fe_SNIa = math.log(Mg_num / Fe_num, 10) - Mg_abundances_solar + Fe_abundances_solar
    Si_over_Fe_SNIa = math.log(Si_num / Fe_num, 10) - Si_abundances_solar + Fe_abundances_solar
    Ca_over_Fe_SNIa = math.log(Ca_num / Fe_num, 10) - Ca_abundances_solar + Fe_abundances_solar
    # plt.plot([-0.3, 0.9], [Mg_over_Fe_SNIa, Mg_over_Fe_SNIa], ls="--", lw=2, label="SNIa")
    # plt.plot([-2, 3], [0, 0], lw=0.5, ls='dotted')
    # plt.xlim(-0.5, 2.2)
    # plt.ylim(-2, 3.5)
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$ [$M_\odot$])')
    # plt.ylabel(r'[Mg/Fe]')
    # plt.tight_layout()
    # plt.savefig('steller_yield_Mg_over_Fe.pdf', dpi=250)
    #
    #
    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(3, figsize=(4, 3.5))
    # plt.xlim(-0.5, 2.2)
    # plt.ylim(-2, 7)
    # i = len(M_list) - 1
    # while i > -1:
    #     plt.plot(M_list[i], O_over_Fe_list[i], label='Z={}'.format(Z_list[i]))
    #     (i) = (i - 1)
    # O_over_Fe_SNIa = math.log(O_num / Fe_num, 10) - O_abundances_solar + Fe_abundances_solar
    # plt.plot([-0.3, 0.9], [O_over_Fe_SNIa, O_over_Fe_SNIa], ls="--", lw=2, label="SNIa")
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$/[$M_\odot$])')
    # plt.ylabel(r'[O/Fe]')
    # plt.tight_layout()
    # #
    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(4, figsize=(4, 3.5))
    # plt.xlim(-0.5, 2.2)
    # plt.ylim(-2, 2)
    # i = len(M_list) - 1
    # while i > -1:
    #     plt.plot(M_list[i], Mg_over_H_list[i], label='Z={}'.format(Z_list[i]))
    #     (i) = (i - 1)
    # plt.plot([-2, 3], [0, 0], lw=0.1)
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$/[$M_\odot$])')
    # plt.ylabel(r'[Mg/H]')
    # plt.tight_layout()
    # #
    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(42, figsize=(4, 3.5))
    # plt.xlim(-0.5, 2.2)
    # plt.ylim(-2, 2)
    # i = len(M_list) - 1
    # while i > -1:
    #     plt.plot(M_list[i], Si_over_H_list[i], label='Z={}'.format(Z_list[i]))
    #     (i) = (i - 1)
    # plt.plot([-2, 3], [0, 0], lw=0.1)
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$/[$M_\odot$])')
    # plt.ylabel(r'[Si/H]')
    # plt.tight_layout()
    # #
    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(41, figsize=(4, 3.5))
    # plt.xlim(-0.5, 2.2)
    # # plt.ylim(-2, 2)
    # i = len(M_list) - 1
    # while i > -1:
    #     plt.plot(M_list[i], C_over_H_list[i], label='Z={}'.format(Z_list[i]))
    #     (i) = (i - 1)
    # plt.plot([-2, 3], [0, 0], lw=0.1)
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$/[$M_\odot$])')
    # plt.ylabel(r'[C/H]')
    # plt.tight_layout()
    # plt.savefig('steller_yield_Mg.pdf', dpi=250)
    #
    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(5, figsize=(4, 3.5))
    # plt.xlim(-0.5, 2.2)
    # plt.ylim(-2, 2)
    # i = len(M_list) - 1
    # while i > -1:
    #     plt.plot(M_list[i], O_over_H_list[i], label='Z={}'.format(Z_list[i]))
    #     (i) = (i - 1)
    # plt.plot([-2, 3], [0, 0], lw=0.1)
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$/[$M_\odot$])')
    # plt.ylabel(r'[O/H]')
    # plt.tight_layout()
    # # plt.savefig('steller_yield_O.pdf', dpi=250)
    #
    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(6, figsize=(4, 3.5))
    # plt.xlim(-0.5, 2.2)
    # plt.ylim(-2, 2)
    # i = len(M_list)-1
    # while i > -1:
    #     plt.plot(M_list[i], Fe_over_H_list[i], label='Z={}'.format(Z_list[i]))
    #     (i) = (i - 1)
    # plt.plot([-2, 3], [0, 0], lw=0.1)
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$/[$M_\odot$])')
    # plt.ylabel(r'[Fe/H]')
    # plt.tight_layout()
    # # plt.savefig('steller_yield_Fe.pdf', dpi=250)
    #
    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(7, figsize=(4, 3.5))
    # plt.xlim(-0.5, 2.2)
    # plt.ylim(-2, 2)
    # i = len(M_list) - 1
    # while i > -1:
    #     plt.plot(M_list[i], Z_over_H_list[i], label='Z={}'.format(Z_list[i]))
    #     (i) = (i - 1)
    # plt.plot([-2, 3], [0, 0], lw=0.1)
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$/[$M_\odot$])')
    # plt.ylabel(r'[Z/H]')
    # plt.title("Number ratio")
    # plt.tight_layout()
    #
    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(8, figsize=(4, 3.5))
    # plt.xlim(-0.5, 2.2)
    # plt.ylim(-2, 2)
    # i = len(M_list) - 1
    # while i > -1:
    #     plt.plot(M_list[i], Z_over_X_list[i], label='Z={}'.format(Z_list[i]))
    #     (i) = (i - 1)
    # plt.plot([-2, 3], [0, 0], lw=0.1)
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$/[$M_\odot$])')
    # plt.ylabel(r'[Z/X]')
    # plt.title("Mass ratio")
    # plt.tight_layout()
    #
    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(11, figsize=(4, 3.5))
    # plt.xlim(-0.5, 2.2)
    # plt.ylim(0.23, 0.6)
    # i = len(M_list) - 1
    # while i > -1:
    #     plt.plot(M_list[i], YYY_list[i], label='Z={}'.format(Z_list[i]))
    #     (i) = (i - 1)
    # # plt.plot([-2, 3], [0.25, 0.25], lw=0.5)
    # plt.legend(prop={'size': 6}, loc='best')
    # plt.xlabel(r'log$_{10}$($M_{\rm *, initial}$/[$M_\odot$])')
    # plt.ylabel('Y')
    # plt.tight_layout()
    # # plt.savefig('steller_yield_Y.pdf', dpi=250)


    ##########
    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(3, 4))

    # i = len(M_list) - 1
    # while i > -1:
    #     axs[0].plot(M_list[i], Z_over_Z0_list[i], lw=(i+2)/2, color=colors[color_list_[i]])
    #     (i) = (i - 1)
    # axs[0].plot([-2, 3], [0, 0], lw=0.7, ls='dotted')
    # # axs[0].set_yticks(np.arange(-1, 2.1, 1))
    # axs[0].set_ylim(-2, 1.6)
    # axs[0].set_ylabel(r'[Z]')
    # 
    # i = len(M_list) - 1
    # while i > -1:
    #     # axs[1].plot(M_list[i], XXX_list[i], lw=(i+2)/2, color=colors[color_list_[i]])
    #     axs[1].plot(M_list[i], YYY_list[i], lw=(i+2)/2, color=colors[color_list_[i]])
    #     # axs[1].plot(M_list[i], ZZZ_list[i], lw=(i+2)/2, color=colors[color_list_[i]])
    #     (i) = (i - 1)
    # axs[1].plot([-2, 3], [0.273, 0.273], lw=0.7, ls='dotted')
    # # axs[1].set_yticks(np.arange(0.2, 0.61, 0.1))
    # axs[1].set_ylim(0.24, 0.605)
    # axs[1].set_xlim(-0.5, 2.2)
    # axs[1].set_ylabel('Y')

    # axs[0].plot([1.3073, 1.3073], [-0.1, 1.7], lw=0.2)
    axs[0].axvspan(1.3073, 3, alpha=0.2, color='red')
    i = len(M_list) - 1
    while i > -1:
        ZZZ = Z_list[i]
        if ZZZ > 0:
            Z_box = round(math.log(ZZZ, 10) - math.log(0.01886, 10), 2)
        else:
            Z_box = -6
        M_list[i].insert(0, math.log(150, 10))
        Mg_over_Fe_list[i].insert(0, Mg_over_Fe_list[i][0])
        axs[0].plot(M_list[i], Mg_over_Fe_list[i], lw=2**i*0.7, label=r'$Z={}$'.format(ZZZ), color='k', ls=['-', 'dashed', 'dotted'][i])
        (i) = (i - 1)
    # axs[0].plot([-0.3, 0.9], [Mg_over_Fe_SNIa, Mg_over_Fe_SNIa], ls="--", lw=1, label="SNIa", c='k')
    # axs[0].plot([-2, 3], [0, 0], lw=0.7, ls='dotted')
    # axs[0].set_yticks(np.arange(-2, 2.1, 2))
    axs[0].set_xlim(0.7, 1.7)
    axs[0].set_ylim(-0.1, 1.7)
    axs[0].set_ylabel(r'[Mg/Fe]')
    axs[0].set_xlabel(r'log$_{10}(M_{\rm *, initial}$ [$M_\odot$])')
    axs[0].legend(prop={'size': 6}, loc='best')

    axs[1].axvspan(1.3073, 3, alpha=0.2, color='red')
    i = len(M_list) - 1
    while i > -1:
        Si_over_Fe_list[i].insert(0, Si_over_Fe_list[i][0])
        axs[1].plot(M_list[i], Si_over_Fe_list[i], lw=2**i*0.7, label=r'$Z={}$'.format(ZZZ),
                    color='k', ls=['-', 'dashed', 'dotted'][i])
        (i) = (i - 1)
    # axs[1].plot([-0.3, 0.9], [Si_over_Fe_SNIa, Si_over_Fe_SNIa], ls="--", lw=1, label="SNIa", c='k')
    # axs[1].plot([-2, 3], [0, 0], lw=0.7, ls='dotted')
    # axs[1].set_yticks(np.arange(-2, 2.1, 2))
    axs[1].set_ylim(-0.1, 1.7)
    axs[1].set_ylabel(r'[Si/Fe]')
    axs[1].set_xlabel(r'log$_{10}(M_{\rm *, initial}$ [$M_\odot$])')
    # axs[1].legend(prop={'size': 6}, loc='best')

    axs[2].axvspan(1.3073, 3, alpha=0.2, color='red')
    i = len(M_list) - 1
    while i > -1:
        Ca_over_Fe_list[i].insert(0, Ca_over_Fe_list[i][0])
        axs[2].plot(M_list[i], Ca_over_Fe_list[i], lw=2**i*0.7, label=r'$Z={}$'.format(ZZZ),
                    color='k', ls=['-', 'dashed', 'dotted'][i])
        (i) = (i - 1)
    # axs[2].plot([-0.3, 0.9], [Ca_over_Fe_SNIa, Ca_over_Fe_SNIa], ls="--", lw=1, label="SNIa", c='k')
    # axs[2].plot([-2, 3], [0, 0], lw=0.7, ls='dotted')
    # axs[2].set_yticks(np.arange(-2, 2.1, 2))
    axs[2].set_ylim(-0.1, 1.7)
    axs[2].set_ylabel(r'[Ca/Fe]')
    axs[2].set_xlabel(r'log$_{10}(M_{\rm *, initial}$ [$M_\odot$])')
    # axs[2].legend(prop={'size': 6}, loc='best')

    plt.tight_layout()
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)
    plt.savefig('stellar_yields.pdf', dpi=250)



    plt.show()
    return


if __name__ == '__main__':
    start_time = time.time()

    Z_list = []
    M_list = []
    eject_mass_list = []
    H_eject_mass_list = []
    He_eject_mass_list = []
    C_eject_mass_list = []
    N_eject_mass_list = []
    O_eject_mass_list = []
    Ne_eject_mass_list = []
    Mg_eject_mass_list = []
    Si_eject_mass_list = []
    S_eject_mass_list = []
    Ca_eject_mass_list = []
    Fe_eject_mass_list = []
    Metal_eject_mass_list = []
    O_over_Mg_list = []
    Mg_over_H_list = []
    Si_over_H_list = []
    C_over_H_list = []
    Fe_over_H_list = []
    O_over_H_list = []
    Z_over_H_list = []
    Z_over_X_list = []
    Z_over_Z0_list = []
    XXX_list = []
    YYY_list = []
    ZZZ_list = []
    Mg_over_Fe_list = []
    Si_over_Fe_list = []
    Ca_over_Fe_list = []
    O_over_Fe_list = []
    yield_table_name = "Kobayashi06" # being "WW95" or "portinari98" or "marigo01"
    function_read_file(yield_table_name)
    funtion_plot_yields()

    plot_lifetime_and_finalmass()
    print(" - Run time: %s -" % round((time.time() - start_time), 2))