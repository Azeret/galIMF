import re
import time
import numpy as np
import math
import matplotlib.pyplot as plt
import element_abundances_solar as solar
import os
from yield_tables import SNIa_yield
import element_weight_table
from scipy.interpolate import UnivariateSpline

SNIa_yield_table = "Iwamoto1999_WDD3"  # "Iwamoto1999_W70"-0.39 #"Seitenzahl2013"-0.16 "Iwamoto1999"-0.24
C_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'C')
O_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'O')
Mg_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Mg')
Al_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Al')
Si_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Si')
S_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'S')
Ar_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Ar')
Ne_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Ne')
Ca_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Ca')
Ti_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Ti')
Cr_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Cr')
Mn_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Mn')
Fe_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Fe')
Ni_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Ni')

solar_H = solar.function_solar_element_abundances("Asplund2009", "H")
solar_He = solar.function_solar_element_abundances("Asplund2009", "He")
solar_C = solar.function_solar_element_abundances("Asplund2009", "C")
solar_C13 = solar.function_solar_element_abundances("Asplund2009", "C13")
solar_N = solar.function_solar_element_abundances("Asplund2009", "N")
solar_O = solar.function_solar_element_abundances("Asplund2009", "O")
solar_O17 = solar.function_solar_element_abundances("Asplund2009", "O17")
solar_O18 = solar.function_solar_element_abundances("Asplund2009", "O18")
solar_Ne = solar.function_solar_element_abundances("Asplund2009", "Ne")
solar_Mg = solar.function_solar_element_abundances("Asplund2009", "Mg")
solar_Al = solar.function_solar_element_abundances("Asplund2009", "Al")
solar_Si = solar.function_solar_element_abundances("Asplund2009", "Si")
solar_S = solar.function_solar_element_abundances("Asplund2009", "S")
solar_Ar = solar.function_solar_element_abundances("Asplund2009", "Ar")
solar_Ca = solar.function_solar_element_abundances("Asplund2009", "Ca")
solar_Ti = solar.function_solar_element_abundances("Asplund2009", "Ti")
solar_Cr = solar.function_solar_element_abundances("Asplund2009", "Cr")
solar_Mn = solar.function_solar_element_abundances("Asplund2009", "Mn")
solar_Fe = solar.function_solar_element_abundances("Asplund2009", "Fe")
solar_Ni = solar.function_solar_element_abundances("Asplund2009", "Ni")
solar_Na = solar.function_solar_element_abundances("Asplund2009", "Na")
solar_Y = solar.function_solar_element_abundances("Asplund2009", "Y")
solar_Ba = solar.function_solar_element_abundances("Asplund2009", "Ba")
solar_Ce = solar.function_solar_element_abundances("Asplund2009", "Ce")
solar_Eu = solar.function_solar_element_abundances("Asplund2009", "Eu")


def function_read_file(yield_table_name):
    ####################
    ### read in file ###
    ####################
    if yield_table_name == "portinari98":
        file_yield = open('yield_tables/agb_and_massive_stars_portinari98_marigo01_gce_totalyields.txt', 'r')
        # Use net yields of Portinari and Marigo
        # Net yields with masses up to 7Msun are from Marigo, above those of Portinari are taken.
        # Only isotopes are selected which are available in both yield sets and go up to Fe.
        # Initial masses go from the lowest mass available up to 100Msun. (100!!!)
        # Yield set ID M01P98 in Ritter et al. 2017.
        # References: Marigo et al. 2001, http://ukads.nottingham.ac.uk/abs/2001A%26A...370..194M
        #		      Portinari et al. 1998, http://ukads.nottingham.ac.uk/abs/1998A%26A...334..505P
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Kobayashi06":
        file_yield = open('yield_tables/agb_and_massive_stars_Kobayashi06_marigo01_gce_totalyields.txt', 'r')
        # Use net yields of Kobayashi C., Umeda H., Nomoto K., Tominaga N., Ohkubo T., 2006, ApJ, 653, 1145
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "WW95":
        file_yield = open('yield_tables/massive_stars_WW95_totalyields.txt', 'r')
        # Use net yields of Woosley S. E., Weaver T. A., 1995, ApJS, 101, 181 (WW95)
        # Use WW95 model B which has the highest [Mg/Fe].
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "marigo01":
        file_yield = open('yield_tables/agb_marigo01_totalyields.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Karakas10":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_K06_0.0HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_R000":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_LC18_R000.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_R":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_LC18_R_mix.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_M000":
        file_yield = open('yield_tables/agb_and_massive_stars_K10C20_lc18_r0_M.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_M150":
        file_yield = open('yield_tables/agb_and_massive_stars_K10C20_lc18_r150_M.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_M300":
        file_yield = open('yield_tables/agb_and_massive_stars_K10C20_lc18_r300_M.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_R150":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_LC18_R150.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_R300":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_LC18_R300.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_N13_0_0_HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_HNe":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_N13_1_0_HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_HNe_05":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_N13_0_5_HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "nugrid_FRUITY":
        file_yield = open('yield_tables/agb_and_massive_stars_nugrid_FRUITY.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "nugrid_N13":
        file_yield = open('yield_tables/agb_and_massive_stars_nugrid_N13.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "nugrid_K06":
        file_yield = open('yield_tables/agb_and_massive_stars_nugrid_K06.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "nugrid_MESAonly_ye":
        file_yield = open('yield_tables/agb_and_massive_stars_nugrid_MESAonly_ye.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_ZY_hypernova":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_2.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_ZY_CCSN":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_4.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_ZY_CCSN_popIII":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_4_correctPopIII.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_ZY_CCSN_ECap":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_5.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "K10_N13_HegerPopIII":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_HegerPopIII.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_ZY_6":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_6.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "C15_N13_HNe10":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_N13_1_0_HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "K10_K06_HNe10":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_K06_1.0HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "C15_N13_HNe00":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_N13_0_0_HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "K10_K06_HNe00":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_K06_0.0HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "K10_N13_HNe00":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_0.0HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "popIII_heger10":
        file_yield = open('yield_tables/popIII_heger10.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "popIII_N13":
        file_yield = open('yield_tables/popIII_N13.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()

    ###########################
    ### extract information ###
    ###########################
    #
    H_relative_line_number = function_get_element_line_number(data, 'H-1')
    He_relative_line_number = function_get_element_line_number(data, 'He-4')
    C12_relative_line_number = function_get_element_line_number(data, 'C-12')
    C13_relative_line_number = function_get_element_line_number(data, 'C-13')
    N14_relative_line_number = function_get_element_line_number(data, 'N-14')
    N15_relative_line_number = function_get_element_line_number(data, 'N-15')
    O16_relative_line_number = function_get_element_line_number(data, 'O-16')
    O17_relative_line_number = function_get_element_line_number(data, 'O-17')
    O18_relative_line_number = function_get_element_line_number(data, 'O-18')
    Ne20_relative_line_number = function_get_element_line_number(data, 'Ne-20')
    Ne21_relative_line_number = function_get_element_line_number(data, 'Ne-21')
    Na_relative_line_number = function_get_element_line_number(data, 'Na-23')
    Mg24_relative_line_number = function_get_element_line_number(data, 'Mg-24')
    Mg25_relative_line_number = function_get_element_line_number(data, 'Mg-25')
    Mg26_relative_line_number = function_get_element_line_number(data, 'Mg-26')
    Al_relative_line_number = function_get_element_line_number(data, 'Al-27')
    Si28_relative_line_number = function_get_element_line_number(data, 'Si-28')
    Si29_relative_line_number = function_get_element_line_number(data, 'Si-29')
    Si30_relative_line_number = function_get_element_line_number(data, 'Si-30')
    S32_relative_line_number = function_get_element_line_number(data, 'S-32')
    S33_relative_line_number = function_get_element_line_number(data, 'S-33')
    S34_relative_line_number = function_get_element_line_number(data, 'S-34')
    S36_relative_line_number = function_get_element_line_number(data, 'S-36')
    Ar36_relative_line_number = function_get_element_line_number(data, 'Ar-36')
    Ar38_relative_line_number = function_get_element_line_number(data, 'Ar-38')
    Ar40_relative_line_number = function_get_element_line_number(data, 'Ar-40')
    Ca40_relative_line_number = function_get_element_line_number(data, 'Ca-40')
    Ca42_relative_line_number = function_get_element_line_number(data, 'Ca-42')
    Ca43_relative_line_number = function_get_element_line_number(data, 'Ca-43')
    Ca44_relative_line_number = function_get_element_line_number(data, 'Ca-44')
    Ca46_relative_line_number = function_get_element_line_number(data, 'Ca-46')
    Ca48_relative_line_number = function_get_element_line_number(data, 'Ca-48')
    Ti46_relative_line_number = function_get_element_line_number(data, 'Ti-46')
    Ti47_relative_line_number = function_get_element_line_number(data, 'Ti-47')
    Ti48_relative_line_number = function_get_element_line_number(data, 'Ti-48')
    Ti49_relative_line_number = function_get_element_line_number(data, 'Ti-49')
    Ti50_relative_line_number = function_get_element_line_number(data, 'Ti-50')
    Cr50_relative_line_number = function_get_element_line_number(data, 'Cr-50')
    Cr52_relative_line_number = function_get_element_line_number(data, 'Cr-52')
    Cr53_relative_line_number = function_get_element_line_number(data, 'Cr-53')
    Cr54_relative_line_number = function_get_element_line_number(data, 'Cr-54')
    Mn_relative_line_number = function_get_element_line_number(data, 'Mn-55')
    Fe54_relative_line_number = function_get_element_line_number(data, 'Fe-54')
    Fe56_relative_line_number = function_get_element_line_number(data, 'Fe-56')
    Fe57_relative_line_number = function_get_element_line_number(data, 'Fe-57')
    Fe58_relative_line_number = function_get_element_line_number(data, 'Fe-58')
    Ni58_relative_line_number = function_get_element_line_number(data, 'Ni-58')
    Ni60_relative_line_number = function_get_element_line_number(data, 'Ni-60')
    Ni61_relative_line_number = function_get_element_line_number(data, 'Ni-61')
    Ni62_relative_line_number = function_get_element_line_number(data, 'Ni-62')
    Ni64_relative_line_number = function_get_element_line_number(data, 'Ni-64')
    # Ce_relative_line_number = function_get_element_line_number(data, 'Ce-140')
    #
    global M_list, Z_list, eject_mass_list, lifetime_list, Mfinal_list, H_eject_mass_list, He_eject_mass_list, C_eject_mass_list, C13_eject_mass_list, \
        N_eject_mass_list, O_eject_mass_list, O17_eject_mass_list, O18_eject_mass_list, Ne_eject_mass_list, Na_eject_mass_list, Mg_eject_mass_list, Al_eject_mass_list, Si_eject_mass_list, \
        S_eject_mass_list, Ar_eject_mass_list, Ca_eject_mass_list, Ti_eject_mass_list, Cr_eject_mass_list, Mn_eject_mass_list, Fe_eject_mass_list, Ni_eject_mass_list, Metal_eject_mass_list
    global O_over_Mg_list, C13_over_O17_list, C_over_O_list, C13_over_C_list, O17_over_O_list, Al_over_Mg_list, Al_over_O_list, Mg_over_Fe_list, Al_over_Fe_list, Mg_over_H_list, Al_over_H_list, \
        Cr_over_H_list, Mn_over_H_list, Fe_over_H_list, Ni_over_H_list, \
        C_over_H_list, N_over_H_list, O_over_H_list, C13_over_H_list, O17_over_H_list, O18_over_H_list, Si_over_H_list, Ar_over_H_list, Ca_over_H_list, Ti_over_H_list, \
        Z_over_H_list, O_over_Fe_list, Si_over_Fe_list, Ca_over_Fe_list, Mn_over_Fe_list, Cr_over_Fe_list, Ni_over_Fe_list, Ti_over_Fe_list
    #
    i = 0
    while i < len(data):
        line_i = str.split(data[i])
        if line_i[1] == 'Table:':  # Here select the lines being: ['H', 'Table:', '(M=xxx,Z=xxx)']
            line_H = str.split(data[i + H_relative_line_number])
            line_He = str.split(data[i + He_relative_line_number])
            line_C12 = str.split(data[i + C12_relative_line_number])
            line_C13 = str.split(data[i + C13_relative_line_number])
            line_N14 = str.split(data[i + N14_relative_line_number])
            line_N15 = str.split(data[i + N15_relative_line_number])
            line_O16 = str.split(data[i + O16_relative_line_number])
            line_O17 = str.split(data[i + O17_relative_line_number])
            line_O18 = str.split(data[i + O18_relative_line_number])
            line_Ne20 = str.split(data[i + Ne20_relative_line_number])
            line_Ne21 = str.split(data[i + Ne21_relative_line_number])
            line_Na = str.split(data[i + Na_relative_line_number])
            line_Mg24 = str.split(data[i + Mg24_relative_line_number])
            line_Mg25 = str.split(data[i + Mg25_relative_line_number])
            line_Mg26 = str.split(data[i + Mg26_relative_line_number])
            line_Al = str.split(data[i + Al_relative_line_number])
            line_Si28 = str.split(data[i + Si28_relative_line_number])
            line_Si29 = str.split(data[i + Si29_relative_line_number])
            line_Si30 = str.split(data[i + Si30_relative_line_number])
            line_S32 = str.split(data[i + S32_relative_line_number])
            line_S33 = str.split(data[i + S33_relative_line_number])
            line_S34 = str.split(data[i + S34_relative_line_number])
            line_S36 = str.split(data[i + S36_relative_line_number])
            line_Ar36 = str.split(data[i + Ar36_relative_line_number])
            line_Ar38 = str.split(data[i + Ar38_relative_line_number])
            line_Ar40 = str.split(data[i + Ar40_relative_line_number])
            line_Ca40 = str.split(data[i + Ca40_relative_line_number])
            line_Ca42 = str.split(data[i + Ca42_relative_line_number])
            line_Ca43 = str.split(data[i + Ca43_relative_line_number])
            line_Ca44 = str.split(data[i + Ca44_relative_line_number])
            line_Ca46 = str.split(data[i + Ca46_relative_line_number])
            line_Ca48 = str.split(data[i + Ca48_relative_line_number])
            line_Ti46 = str.split(data[i + Ti46_relative_line_number])
            line_Ti47 = str.split(data[i + Ti47_relative_line_number])
            line_Ti48 = str.split(data[i + Ti48_relative_line_number])
            line_Ti49 = str.split(data[i + Ti49_relative_line_number])
            line_Ti50 = str.split(data[i + Ti50_relative_line_number])
            line_Cr50 = str.split(data[i + Cr50_relative_line_number])
            line_Cr52 = str.split(data[i + Cr52_relative_line_number])
            line_Cr53 = str.split(data[i + Cr53_relative_line_number])
            line_Cr54 = str.split(data[i + Cr54_relative_line_number])
            line_Mn = str.split(data[i + Mn_relative_line_number])
            line_Fe54 = str.split(data[i + Fe54_relative_line_number])
            line_Fe56 = str.split(data[i + Fe56_relative_line_number])
            line_Fe57 = str.split(data[i + Fe57_relative_line_number])
            line_Fe58 = str.split(data[i + Fe58_relative_line_number])
            line_Ni58 = str.split(data[i + Ni58_relative_line_number])
            line_Ni60 = str.split(data[i + Ni60_relative_line_number])
            line_Ni61 = str.split(data[i + Ni61_relative_line_number])
            line_Ni62 = str.split(data[i + Ni62_relative_line_number])
            line_Ni64 = str.split(data[i + Ni64_relative_line_number])
            line_Lifetime = str.split(data[i + 1])
            lifetime = function_get_Mfinal_and_Lifetime(line_Lifetime[2])
            (Z_ini, M_ini) = function_get_Z_M(line_i[2])  # get the initial mass and metallicity of the star
            if yield_table_name == "popIII_heger10":
                Mfinal = 0.1 * M_ini
            else:
                line_Mfinal = str.split(data[i + 2])
                Mfinal = function_get_Mfinal_and_Lifetime(line_Mfinal[2])
            ejecta_mass = round((M_ini - Mfinal), 5)  ####################

            if yield_table_name == "Nomoto_ZY_CCSN" and M_ini == 8.0 and Z_ini == 0.05:
                scale_factor = (8 - 1.12) / (6 - 0.929)
            elif yield_table_name == "Nomoto_ZY_CCSN" and M_ini == 8.0 and Z_ini == 0.02:
                scale_factor = (8 - 1.12) / (6 - 0.929)
            elif yield_table_name == "Nomoto_ZY_CCSN" and M_ini == 8.0 and Z_ini == 0.008:
                scale_factor = (8 - 1.12) / (6 - 0.948)
            elif yield_table_name == "Nomoto_ZY_CCSN" and M_ini == 8.0 and Z_ini == 0.004:
                scale_factor = (8 - 1.12) / (6 - 0.977)
            elif yield_table_name == "Nomoto_ZY_CCSN" and M_ini == 8.0 and Z_ini == 0.001:
                scale_factor = (8 - 1.12) / (6 - 0.9879)
            elif yield_table_name == "Nomoto_ZY_CCSN" and M_ini == 8.0 and Z_ini == 0.0001:
                scale_factor = (8 - 1.12) / (6 - 1.006)
            elif yield_table_name == "Nomoto_ZY_CCSN" and M_ini == 100:
                scale_factor = (100 - 10) / (40 - 2.5)
            else:
                scale_factor = 1
            H_mass = function_get_element_mass(line_H[1]) * scale_factor
            H_num = H_mass / element_weight_table.function_element_weight("H")
            He_mass = function_get_element_mass(line_He[1]) * scale_factor
            C12_mass = function_get_element_mass(line_C12[1]) * scale_factor
            C12_num = C12_mass / element_weight_table.function_element_weight("C12")
            C13_mass = function_get_element_mass(line_C13[1]) * scale_factor
            C13_num = C13_mass / element_weight_table.function_element_weight("C13")
            C_mass = C12_mass + C13_mass
            C_num = C12_num + C13_num
            N14_mass = function_get_element_mass(line_N14[1]) * scale_factor
            N15_mass = function_get_element_mass(line_N15[1]) * scale_factor
            N14_num = N14_mass / element_weight_table.function_element_weight("N14")
            N15_num = N15_mass / element_weight_table.function_element_weight("N15")
            N_mass = N14_mass + N15_mass
            N_num = N14_num + N15_num
            O16_mass = function_get_element_mass(line_O16[1]) * scale_factor
            O17_mass = function_get_element_mass(line_O17[1]) * scale_factor
            O18_mass = function_get_element_mass(line_O18[1]) * scale_factor
            O_mass = O16_mass + O17_mass + O18_mass
            O16_num = O16_mass / element_weight_table.function_element_weight("O16")
            O17_num = O17_mass / element_weight_table.function_element_weight("O17")
            O18_num = O18_mass / element_weight_table.function_element_weight("O18")
            O_num = O16_num + O17_num + O18_num
            Ne20_mass = function_get_element_mass(line_Ne20[1]) * scale_factor
            Ne21_mass = function_get_element_mass(line_Ne21[1]) * scale_factor
            Ne_mass = Ne20_mass + Ne21_mass
            Ne20_num = Ne20_mass / element_weight_table.function_element_weight("Ne20")
            Ne21_num = Ne21_mass / element_weight_table.function_element_weight("Ne21")
            Ne_num = Ne20_num + Ne21_num
            Na_mass = function_get_element_mass(line_Na[1]) * scale_factor
            Na_num = Na_mass / element_weight_table.function_element_weight("Na")
            # if 8 < M_ini:
            #     Mg_mass = function_get_element_mass(line_Mg[1]) * scale_factor * 3
            #     # Francois 2004 is based on WW95, for more recent tables, * 7 is applied for all massive stars.
            # # elif 20 < M_ini:
            # #     Mg_mass = function_get_element_mass(line_Mg[1]) * scale_factor / 2  # Francois 2004 for WW95 table
            # else:
            #     Mg_mass = function_get_element_mass(line_Mg[1]) * scale_factor
            Mg24_mass = function_get_element_mass(line_Mg24[1]) * scale_factor
            Mg25_mass = function_get_element_mass(line_Mg25[1]) * scale_factor
            Mg26_mass = function_get_element_mass(line_Mg26[1]) * scale_factor
            Mg_mass = Mg24_mass + Mg25_mass + Mg26_mass
            Mg24_num = Mg24_mass / element_weight_table.function_element_weight("Mg24")
            Mg25_num = Mg25_mass / element_weight_table.function_element_weight("Mg25")
            Mg26_num = Mg26_mass / element_weight_table.function_element_weight("Mg26")
            Mg_num = Mg24_num + Mg25_num + Mg26_num
            Al_mass = function_get_element_mass(line_Al[1]) * scale_factor
            Al_num = Al_mass / element_weight_table.function_element_weight("Al")
            Si28_mass = function_get_element_mass(line_Si28[1]) * scale_factor
            Si29_mass = function_get_element_mass(line_Si29[1]) * scale_factor
            Si30_mass = function_get_element_mass(line_Si30[1]) * scale_factor
            Si_mass = Si28_mass + Si29_mass + Si30_mass
            Si28_num = Si28_mass / element_weight_table.function_element_weight("Si28")
            Si29_num = Si29_mass / element_weight_table.function_element_weight("Si29")
            Si30_num = Si30_mass / element_weight_table.function_element_weight("Si30")
            Si_num = Si28_num + Si29_num + Si30_num
            S32_mass = function_get_element_mass(line_S32[1]) * scale_factor
            S33_mass = function_get_element_mass(line_S33[1]) * scale_factor
            S34_mass = function_get_element_mass(line_S34[1]) * scale_factor
            S36_mass = function_get_element_mass(line_S36[1]) * scale_factor
            S_mass = S32_mass + S33_mass + S34_mass + S36_mass
            S32_num = S32_mass / element_weight_table.function_element_weight("S32")
            S33_num = S33_mass / element_weight_table.function_element_weight("S33")
            S34_num = S34_mass / element_weight_table.function_element_weight("S34")
            S36_num = S36_mass / element_weight_table.function_element_weight("S36")
            S_num = S32_num + S33_num + S34_num + S36_num
            Ar36_mass = function_get_element_mass(line_Ar36[1]) * scale_factor
            Ar38_mass = function_get_element_mass(line_Ar38[1]) * scale_factor
            Ar40_mass = function_get_element_mass(line_Ar40[1]) * scale_factor
            Ar_mass = Ar36_mass + Ar38_mass + Ar40_mass
            Ar36_num = Ar36_mass / element_weight_table.function_element_weight("Ar36")
            Ar38_num = Ar38_mass / element_weight_table.function_element_weight("Ar38")
            Ar40_num = Ar40_mass / element_weight_table.function_element_weight("Ar40")
            Ar_num = max(Ar36_num + Ar38_num + Ar40_num, 1e-30)
            if yield_table_name == "Limongi_R000" or "Limongi_R150" or "Limongi_R300" or "Nomoto_ZY_CCSN":
                if M_ini > 10:
                    Ca40_mass = function_get_element_mass(line_Ca40[1]) * scale_factor
                    Ca42_mass = function_get_element_mass(line_Ca42[1]) * scale_factor
                    Ca43_mass = function_get_element_mass(line_Ca43[1]) * scale_factor
                    Ca44_mass = function_get_element_mass(line_Ca44[1]) * scale_factor
                    Ca46_mass = function_get_element_mass(line_Ca46[1]) * scale_factor
                    Ca48_mass = function_get_element_mass(line_Ca48[1]) * scale_factor
                    Ca_mass = Ca40_mass + Ca42_mass + Ca43_mass + Ca44_mass + Ca46_mass + Ca48_mass
                    Ca40_num = Ca40_mass / element_weight_table.function_element_weight("Ca40")
                    Ca42_num = Ca42_mass / element_weight_table.function_element_weight("Ca42")
                    Ca43_num = Ca43_mass / element_weight_table.function_element_weight("Ca43")
                    Ca44_num = Ca44_mass / element_weight_table.function_element_weight("Ca44")
                    Ca46_num = Ca46_mass / element_weight_table.function_element_weight("Ca46")
                    Ca48_num = Ca48_mass / element_weight_table.function_element_weight("Ca48")
                    Ca_num = Ca40_num + Ca42_num + Ca43_num + Ca44_num + Ca46_num + Ca48_num
                # Ca yield is 1.000E-30 in Karakas10 (adopted by Limongi and Nomoto) for low mass stars.
                elif Z_ini < 0.0004:
                    Ca_mass = function_get_element_mass(line_Fe56[1]) / 56 * 40 * 10 ** (
                                solar_Ca - solar_Fe + 0.052) * scale_factor
                    Ca_num = Ca_mass / element_weight_table.function_element_weight("Ca40")
                else:
                    Ca_mass = function_get_element_mass(line_Fe56[1]) / 56 * 40 * 10 ** (
                                solar_Ca - solar_Fe + 0.002) * scale_factor
                    Ca_num = Ca_mass / element_weight_table.function_element_weight("Ca40")
            else:
                Ca40_mass = function_get_element_mass(line_Ca40[1]) * scale_factor
                Ca42_mass = function_get_element_mass(line_Ca42[1]) * scale_factor
                Ca43_mass = function_get_element_mass(line_Ca43[1]) * scale_factor
                Ca44_mass = function_get_element_mass(line_Ca44[1]) * scale_factor
                Ca46_mass = function_get_element_mass(line_Ca46[1]) * scale_factor
                Ca48_mass = function_get_element_mass(line_Ca48[1]) * scale_factor
                Ca_mass = Ca40_mass + Ca42_mass + Ca43_mass + Ca44_mass + Ca46_mass + Ca48_mass
                Ca40_num = Ca40_mass / element_weight_table.function_element_weight("Ca40")
                Ca42_num = Ca42_mass / element_weight_table.function_element_weight("Ca42")
                Ca43_num = Ca43_mass / element_weight_table.function_element_weight("Ca43")
                Ca44_num = Ca44_mass / element_weight_table.function_element_weight("Ca44")
                Ca46_num = Ca46_mass / element_weight_table.function_element_weight("Ca46")
                Ca48_num = Ca48_mass / element_weight_table.function_element_weight("Ca48")
                Ca_num = Ca40_num + Ca42_num + Ca43_num + Ca44_num + Ca46_num + Ca48_num
            Ti46_mass = function_get_element_mass(line_Ti46[1]) * scale_factor
            Ti47_mass = function_get_element_mass(line_Ti47[1]) * scale_factor
            Ti48_mass = function_get_element_mass(line_Ti48[1]) * scale_factor
            Ti49_mass = function_get_element_mass(line_Ti49[1]) * scale_factor
            Ti50_mass = function_get_element_mass(line_Ti50[1]) * scale_factor
            Ti_mass = Ti46_mass + Ti47_mass + Ti48_mass + Ti49_mass + Ti50_mass
            Ti46_num = Ti46_mass / element_weight_table.function_element_weight("Ti46")
            Ti47_num = Ti47_mass / element_weight_table.function_element_weight("Ti47")
            Ti48_num = Ti48_mass / element_weight_table.function_element_weight("Ti48")
            Ti49_num = Ti49_mass / element_weight_table.function_element_weight("Ti49")
            Ti50_num = Ti50_mass / element_weight_table.function_element_weight("Ti50")
            Ti_num = max(Ti46_num + Ti47_num + Ti48_num + Ti49_num + Ti50_num, 1e-30)
            Cr50_mass = function_get_element_mass(line_Cr50[1]) * scale_factor
            Cr52_mass = function_get_element_mass(line_Cr52[1]) * scale_factor
            Cr53_mass = function_get_element_mass(line_Cr53[1]) * scale_factor
            Cr54_mass = function_get_element_mass(line_Cr54[1]) * scale_factor
            Cr_mass = Cr50_mass + Cr52_mass + Cr53_mass + Cr54_mass
            Cr50_num = Cr50_mass / element_weight_table.function_element_weight("Cr50")
            Cr52_num = Cr52_mass / element_weight_table.function_element_weight("Cr52")
            Cr53_num = Cr53_mass / element_weight_table.function_element_weight("Cr53")
            Cr54_num = Cr54_mass / element_weight_table.function_element_weight("Cr54")
            Cr_num = max(Cr50_num + Cr52_num + Cr53_num + Cr54_num, 1e-30)
            Mn_mass = function_get_element_mass(line_Mn[1]) * scale_factor
            Mn_num = max(Mn_mass / element_weight_table.function_element_weight("Mn"), 1e-30)
            Ni58_mass = function_get_element_mass(line_Ni58[1]) * scale_factor
            Ni60_mass = function_get_element_mass(line_Ni60[1]) * scale_factor
            Ni61_mass = function_get_element_mass(line_Ni61[1]) * scale_factor
            Ni62_mass = function_get_element_mass(line_Ni62[1]) * scale_factor
            Ni64_mass = function_get_element_mass(line_Ni64[1]) * scale_factor
            Ni_mass = Ni58_mass + Ni60_mass + Ni61_mass + Ni62_mass + Ni64_mass
            Ni58_num = Ni58_mass / element_weight_table.function_element_weight("Ni58")
            Ni60_num = Ni60_mass / element_weight_table.function_element_weight("Ni60")
            Ni61_num = Ni61_mass / element_weight_table.function_element_weight("Ni61")
            Ni62_num = Ni62_mass / element_weight_table.function_element_weight("Ni62")
            Ni64_num = Ni64_mass / element_weight_table.function_element_weight("Ni64")
            Ni_num = Ni58_num + Ni60_num + Ni61_num + Ni62_num + Ni64_num
            if Si_mass == 0:
                Si_mass = 1e-30
            if S_mass == 0:
                S_mass = 1e-30
            if Ar_mass == 0:
                Ar_mass = 1e-30
            if Ca_mass == 0:
                Ca_mass = 1e-30
            if Cr_mass == 0:
                Cr_mass = 1e-30
            if Mn_mass == 0:
                Mn_mass = 1e-30
            if Ni_mass == 0:
                Ni_mass = 1e-30
            if Ti_mass == 0:
                Ti_mass = 1e-30
            Fe54_mass = function_get_element_mass(line_Fe54[1]) * scale_factor
            Fe56_mass = function_get_element_mass(line_Fe56[1]) * scale_factor
            Fe57_mass = function_get_element_mass(line_Fe57[1]) * scale_factor
            Fe58_mass = function_get_element_mass(line_Fe58[1]) * scale_factor
            Fe_mass = Fe54_mass + Fe56_mass + Fe57_mass + Fe58_mass
            Fe54_num = Fe54_mass / element_weight_table.function_element_weight("Fe54")
            Fe56_num = Fe56_mass / element_weight_table.function_element_weight("Fe56")
            Fe57_num = Fe57_mass / element_weight_table.function_element_weight("Fe57")
            Fe58_num = Fe58_mass / element_weight_table.function_element_weight("Fe58")
            Fe_num = Fe54_num + Fe56_num + Fe57_num + Fe58_num
            O_over_Mg = math.log(O_num / Mg_num, 10) - solar_O + solar_Mg
            C13_over_O17 = math.log(C13_num / O17_num, 10)
            C_over_O = math.log(C_num / O_num, 10)
            C13_over_C = math.log(C13_num / C_num, 10)
            O17_over_O = math.log(O17_num / O_num, 10)
            Al_over_Mg = math.log(Al_num / Mg_num, 10) - solar_Al + solar_Mg
            Al_over_O = math.log(Al_num / O_num, 10) - solar_Al + solar_O
            Ne_over_H = math.log(Ne_num / H_num, 10) - solar_Ne + solar_H
            Na_over_H = math.log(Na_num / H_num, 10) - solar_Na + solar_H
            Mg_over_H = math.log(Mg_num / H_num, 10) - solar_Mg + solar_H
            Al_over_H = math.log(Al_num / H_num, 10) - solar_Al + solar_H
            Si_over_H = math.log(Si_num / H_num, 10) - solar_Si + solar_H
            S_over_H = math.log(S_num / H_num, 10) - solar_S + solar_H
            Ar_over_H = math.log(Ar_num / H_num, 10) - solar_Ar + solar_H
            Ca_over_H = math.log(Ca_num / H_num, 10) - solar_Ca + solar_H
            Ti_over_H = math.log(Ti_num / H_num, 10) - solar_Ti + solar_H
            Cr_over_H = math.log(Cr_num / H_num, 10) - solar_Cr + solar_H
            Mn_over_H = math.log(Mn_num / H_num, 10) - solar_Mn + solar_H
            Fe_over_H = math.log(Fe_num / H_num, 10) - solar_Fe + solar_H
            Ni_over_H = math.log(Ni_num / H_num, 10) - solar_Ni + solar_H
            C_over_H = math.log(C_num / H_num, 10) - solar_C + solar_H
            N_over_H = math.log(N_num / H_num, 10) - solar_N + solar_H
            O_over_H = math.log(O_num / H_num, 10) - solar_O + solar_H
            C13_over_H = math.log(C13_num / H_num, 10) - solar_C13 + solar_H
            O17_over_H = math.log(O17_num / H_num, 10) - solar_O17 + solar_H
            O18_over_H = math.log(O18_num / H_num, 10) - solar_O18 + solar_H
            Mg_over_Fe = math.log(Mg_num / Fe_num, 10) - solar_Mg + solar_Fe
            Al_over_Fe = math.log(Al_num / Fe_num, 10) - solar_Al + solar_Fe
            Si_over_Fe = math.log(Si_num / Fe_num, 10) - solar_Si + solar_Fe
            Ca_over_Fe = math.log(Ca_num / Fe_num, 10) - solar_Ca + solar_Fe
            Ti_over_Fe = math.log(Ti_num / Fe_num, 10) - solar_Ti + solar_Fe
            Cr_over_Fe = math.log(Cr_num / Fe_num, 10) - solar_Cr + solar_Fe
            Mn_over_Fe = math.log(Mn_num / Fe_num, 10) - solar_Mn + solar_Fe
            Ni_over_Fe = math.log(Ni_num / Fe_num, 10) - solar_Ni + solar_Fe
            O_over_Fe = math.log(O_num / Fe_num, 10) - solar_O + solar_Fe
            Metal_mass = round((ejecta_mass - H_mass - He_mass), 5)  ####################
            if Metal_mass < 0 or Metal_mass == 0:
                print("Warning: Metal_mass=", Metal_mass, "<0")
                print("check stellar yield table with metallicity and mass being:", Z_ini, "&", M_ini)
                Metal_mass = 1e-6
            Z_over_H = math.log(Metal_mass / H_mass, 10) - math.log(0.0134 / 0.7381, 10)
            if len(Z_list) == 0:
                Z_list.append(Z_ini)
                Z_n = 0
                M_list.append([])
                eject_mass_list.append([])
                lifetime_list.append([])
                Mfinal_list.append([])
                H_eject_mass_list.append([])
                He_eject_mass_list.append([])
                C_eject_mass_list.append([])
                C13_eject_mass_list.append([])
                N_eject_mass_list.append([])
                O_eject_mass_list.append([])
                O17_eject_mass_list.append([])
                O18_eject_mass_list.append([])
                Ne_eject_mass_list.append([])
                Na_eject_mass_list.append([])
                Mg_eject_mass_list.append([])
                Al_eject_mass_list.append([])
                Si_eject_mass_list.append([])
                S_eject_mass_list.append([])
                Ar_eject_mass_list.append([])
                Ca_eject_mass_list.append([])
                Ti_eject_mass_list.append([])
                Cr_eject_mass_list.append([])
                Mn_eject_mass_list.append([])
                Fe_eject_mass_list.append([])
                Ni_eject_mass_list.append([])
                Metal_eject_mass_list.append([])
                Z_over_H_list.append([])
                O_over_Mg_list.append([])
                C13_over_O17_list.append([])
                C_over_O_list.append([])
                C13_over_C_list.append([])
                O17_over_O_list.append([])
                Al_over_Mg_list.append([])
                Al_over_O_list.append([])
                O_over_Fe_list.append([])
                Mg_over_Fe_list.append([])
                Al_over_Fe_list.append([])
                Si_over_Fe_list.append([])
                Ca_over_Fe_list.append([])
                Ti_over_Fe_list.append([])
                Cr_over_Fe_list.append([])
                Mn_over_Fe_list.append([])
                Ni_over_Fe_list.append([])
                C_over_H_list.append([])
                N_over_H_list.append([])
                O_over_H_list.append([])
                C13_over_H_list.append([])
                O17_over_H_list.append([])
                O18_over_H_list.append([])
                Ne_over_H_list.append([])
                Na_over_H_list.append([])
                Mg_over_H_list.append([])
                Al_over_H_list.append([])
                Si_over_H_list.append([])
                S_over_H_list.append([])
                Ar_over_H_list.append([])
                Ca_over_H_list.append([])
                Ti_over_H_list.append([])
                Cr_over_H_list.append([])
                Mn_over_H_list.append([])
                Fe_over_H_list.append([])
                Ni_over_H_list.append([])
            if Z_ini != Z_list[-1]:
                Z_list.append(Z_ini)
                Z_n += 1
                M_list.append([])
                eject_mass_list.append([])
                lifetime_list.append([])
                Mfinal_list.append([])
                H_eject_mass_list.append([])
                He_eject_mass_list.append([])
                C_eject_mass_list.append([])
                C13_eject_mass_list.append([])
                N_eject_mass_list.append([])
                O_eject_mass_list.append([])
                O17_eject_mass_list.append([])
                O18_eject_mass_list.append([])
                Ne_eject_mass_list.append([])
                Na_eject_mass_list.append([])
                Mg_eject_mass_list.append([])
                Al_eject_mass_list.append([])
                Si_eject_mass_list.append([])
                S_eject_mass_list.append([])
                Ar_eject_mass_list.append([])
                Ca_eject_mass_list.append([])
                Ti_eject_mass_list.append([])
                Cr_eject_mass_list.append([])
                Mn_eject_mass_list.append([])
                Fe_eject_mass_list.append([])
                Ni_eject_mass_list.append([])
                Metal_eject_mass_list.append([])
                O_over_Mg_list.append([])
                C13_over_O17_list.append([])
                C_over_O_list.append([])
                C13_over_C_list.append([])
                O17_over_O_list.append([])
                Al_over_Mg_list.append([])
                Al_over_O_list.append([])
                O_over_Fe_list.append([])
                Mg_over_Fe_list.append([])
                Al_over_Fe_list.append([])
                Si_over_Fe_list.append([])
                Ca_over_Fe_list.append([])
                Ti_over_Fe_list.append([])
                Cr_over_Fe_list.append([])
                Mn_over_Fe_list.append([])
                Ni_over_Fe_list.append([])
                C_over_H_list.append([])
                N_over_H_list.append([])
                O_over_H_list.append([])
                C13_over_H_list.append([])
                O17_over_H_list.append([])
                O18_over_H_list.append([])
                Ne_over_H_list.append([])
                Na_over_H_list.append([])
                Mg_over_H_list.append([])
                Al_over_H_list.append([])
                Si_over_H_list.append([])
                S_over_H_list.append([])
                Ar_over_H_list.append([])
                Ca_over_H_list.append([])
                Ti_over_H_list.append([])
                Cr_over_H_list.append([])
                Mn_over_H_list.append([])
                Fe_over_H_list.append([])
                Ni_over_H_list.append([])
                Z_over_H_list.append([])
            M_list[Z_n].append(M_ini)
            eject_mass_list[Z_n].append(ejecta_mass)
            lifetime_list[Z_n].append(lifetime)
            Mfinal_list[Z_n].append(Mfinal)
            H_eject_mass_list[Z_n].append(H_mass)
            He_eject_mass_list[Z_n].append(He_mass)
            C_eject_mass_list[Z_n].append(C_mass)
            C13_eject_mass_list[Z_n].append(C13_mass)
            N_eject_mass_list[Z_n].append(N_mass)
            O_eject_mass_list[Z_n].append(O_mass)
            O17_eject_mass_list[Z_n].append(O17_mass)
            O18_eject_mass_list[Z_n].append(O18_mass)
            Ne_eject_mass_list[Z_n].append(Ne_mass)
            Na_eject_mass_list[Z_n].append(Na_mass)
            Mg_eject_mass_list[Z_n].append(Mg_mass)
            Al_eject_mass_list[Z_n].append(Al_mass)
            Si_eject_mass_list[Z_n].append(Si_mass)
            S_eject_mass_list[Z_n].append(S_mass)
            Ar_eject_mass_list[Z_n].append(Ar_mass)
            Ca_eject_mass_list[Z_n].append(Ca_mass)
            Ti_eject_mass_list[Z_n].append(Ti_mass)
            Cr_eject_mass_list[Z_n].append(Cr_mass)
            Mn_eject_mass_list[Z_n].append(Mn_mass)
            Fe_eject_mass_list[Z_n].append(Fe_mass)
            Ni_eject_mass_list[Z_n].append(Ni_mass)
            Metal_eject_mass_list[Z_n].append(Metal_mass)
            O_over_Mg_list[Z_n].append(O_over_Mg)
            C13_over_O17_list[Z_n].append(C13_over_O17)
            C_over_O_list[Z_n].append(C_over_O)
            C13_over_C_list[Z_n].append(C13_over_C)
            O17_over_O_list[Z_n].append(O17_over_O)
            Al_over_Mg_list[Z_n].append(Al_over_Mg)
            Al_over_O_list[Z_n].append(Al_over_O)
            Mg_over_Fe_list[Z_n].append(Mg_over_Fe)
            Al_over_Fe_list[Z_n].append(Al_over_Fe)
            Si_over_Fe_list[Z_n].append(Si_over_Fe)
            Ca_over_Fe_list[Z_n].append(Ca_over_Fe)
            Ti_over_Fe_list[Z_n].append(Ti_over_Fe)
            Cr_over_Fe_list[Z_n].append(Cr_over_Fe)
            Mn_over_Fe_list[Z_n].append(Mn_over_Fe)
            Ni_over_Fe_list[Z_n].append(Ni_over_Fe)
            Ne_over_H_list[Z_n].append(Ne_over_H)
            Na_over_H_list[Z_n].append(Na_over_H)
            Mg_over_H_list[Z_n].append(Mg_over_H)
            Al_over_H_list[Z_n].append(Al_over_H)
            Si_over_H_list[Z_n].append(Si_over_H)
            S_over_H_list[Z_n].append(S_over_H)
            Ar_over_H_list[Z_n].append(Ar_over_H)
            Ca_over_H_list[Z_n].append(Ca_over_H)
            Ti_over_H_list[Z_n].append(Ti_over_H)
            C_over_H_list[Z_n].append(C_over_H)
            N_over_H_list[Z_n].append(N_over_H)
            O_over_H_list[Z_n].append(O_over_H)
            C13_over_H_list[Z_n].append(C13_over_H)
            O17_over_H_list[Z_n].append(O17_over_H)
            O18_over_H_list[Z_n].append(O18_over_H)
            Z_over_H_list[Z_n].append(Z_over_H)
            Fe_over_H_list[Z_n].append(Fe_over_H)
            Cr_over_H_list[Z_n].append(Cr_over_H)
            Mn_over_H_list[Z_n].append(Mn_over_H)
            Ni_over_H_list[Z_n].append(Ni_over_H)
            O_over_Fe_list[Z_n].append(O_over_Fe)
        (i) = (i + 1)
    # print(Z_list)
    # print(O_eject_mass_list[4])
    # print(M_list)
    # print(eject_mass_list)
    # print(H_eject_mass_list)
    # print(Fe_eject_mass_list)

    ###########################
    ### write data to files ###
    ###########################

    write_data()

    return


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


def function_get_mass_grid(yield_table_name):  # read in a grid from 0.08 to 150 Msun
    if yield_table_name == "Kobayashi06":
        file_yield = open('yield_tables/agb_and_massive_stars_Kobayashi06_marigo01_gce_totalyields.txt', 'r')
        # Use net yields of Kobayashi C., Umeda H., Nomoto K., Tominaga N., Ohkubo T., 2006, ApJ, 653, 1145
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "WW95":
        file_yield = open('yield_tables/massive_stars_WW95_totalyields.txt', 'r')
        # Use net yields of Woosley S. E., Weaver T. A., 1995, ApJS, 101, 181 (WW95)
        # Use WW95 model B which has the highest [Mg/Fe].
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "marigo01":
        file_yield = open('yield_tables/agb_marigo01_totalyields.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Karakas10":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_K06_0.0HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_R000":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_LC18_R000.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_R":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_LC18_R_mix.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_M000":
        file_yield = open('yield_tables/agb_and_massive_stars_K10C20_lc18_r0_M.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_M150":
        file_yield = open('yield_tables/agb_and_massive_stars_K10C20_lc18_r150_M.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_M300":
        file_yield = open('yield_tables/agb_and_massive_stars_K10C20_lc18_r300_M.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_R150":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_LC18_R150.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Limongi_R300":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_LC18_R300.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_N13_0_0_HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_HNe":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_N13_1_0_HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_HNe_05":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_N13_0_5_HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "nugrid_FRUITY":
        file_yield = open('yield_tables/agb_and_massive_stars_nugrid_FRUITY.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "nugrid_N13":
        file_yield = open('yield_tables/agb_and_massive_stars_nugrid_N13.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "nugrid_K06":
        file_yield = open('yield_tables/agb_and_massive_stars_nugrid_K06.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "nugrid_MESAonly_ye":
        file_yield = open('yield_tables/agb_and_massive_stars_nugrid_MESAonly_ye.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_ZY_hypernova":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_2.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_ZY_CCSN":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_4.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_ZY_CCSN_popIII":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_4_correctPopIII.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_ZY_CCSN_ECap":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_5.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "K10_N13_HegerPopIII":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_HegerPopIII.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "Nomoto_ZY_6":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_6.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "C15_N13_HNe10":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_N13_1_0_HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "K10_K06_HNe10":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_K06_1.0HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "C15_N13_HNe00":
        file_yield = open('yield_tables/agb_and_massive_stars_C15_N13_0_0_HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "K10_K06_HNe00":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_K06_0.0HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "K10_N13_HNe00":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_N13_0.0HNe.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "popIII_heger10":
        file_yield = open('yield_tables/popIII_heger10.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    elif yield_table_name == "popIII_N13":
        file_yield = open('yield_tables/popIII_N13.txt', 'r')
        data = file_yield.readlines()
        file_yield.close()
    pattern = r"\(M=([\d.]+)"
    m_values = []
    i = 0
    while i < len(data):
        line_i = str.split(data[i])
        if line_i[1] == 'Table:':  # Here select the lines being: ['H', 'Table:', '(M=xxx,Z=xxx)']
            match = re.search(pattern, line_i[2])
            if match:
                m_values.append(float(match.group(1)))
        (i) = (i + 1)
    mass_grid = [0.08] + sorted(set(m_values)) + [150]
    return mass_grid


def time_list(list, number, mass):
    new_list = []
    for i in range(len(list)):
        new_list.append(list[i] * number)
    return new_list


def time_list_C(list, number, mass):
    new_list = []
    for i in range(len(list)):
        if mass[i] > 25 or 3 < mass[i] < 8:
            new_list.append(list[i] * (number ** 1))
        else:
            new_list.append(list[i])
    return new_list


def time_list_N(list, number, mass):
    new_list = []
    for i in range(len(list)):
        if mass[i] > 8 or mass[i] < 3:
            new_list.append(list[i] * number)
        else:
            new_list.append(list[i])
    return new_list


def time_list_O(list, number, mass):
    new_list = []
    for i in range(len(list)):
        if mass[i] < 9 or mass[i] > 39:
            new_list.append(list[i] * number)
        else:
            new_list.append(list[i])
    return new_list


def time_list_Fe(list, number, mass):
    new_list = []
    for i in range(len(list)):
        if mass[i] < 9 or mass[i] > 19:
            new_list.append(list[i] * number)
        else:
            new_list.append(list[i])
    return new_list


def time_list_alpha(list, number, mass):
    new_list = []
    for i in range(len(list)):
        if mass[i] < 9:
            new_list.append(list[i] * number)
        else:
            new_list.append(list[i])
    return new_list


def write_data():
    global M_list, Z_list, eject_mass_list, lifetime_list, Mfinal_list, H_eject_mass_list, He_eject_mass_list, C_eject_mass_list, C13_eject_mass_list, \
        N_eject_mass_list, O_eject_mass_list, O17_eject_mass_list, O18_eject_mass_list, Ne_eject_mass_list, Na_eject_mass_list, Mg_eject_mass_list, Al_eject_mass_list, Si_eject_mass_list, \
        S_eject_mass_list, Ar_eject_mass_list, Ca_eject_mass_list, Ti_eject_mass_list, Cr_eject_mass_list, Mn_eject_mass_list, Fe_eject_mass_list, Ni_eject_mass_list, Metal_eject_mass_list

    if make_dir == True:
        os.makedirs('yield_tables_fractional/rearranged___/setllar_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_lifetime_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Mfinal_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_H_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_He_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_C_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_C13_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_N_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_O_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_O17_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_O18_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Ne_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Na_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Mg_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Al_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Si_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_S_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Ar_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Ca_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Ti_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Cr_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Mn_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Fe_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Ni_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Y_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Ba_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Ce_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        # os.makedirs('yield_tables_fractional/rearranged___/setllar_Eu_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
        os.makedirs('yield_tables_fractional/rearranged___/setllar_Metal_eject_mass_from_{}'.format(yield_table_name),
                    exist_ok=True)
    mass_grid = function_get_mass_grid(yield_table_name)
    print("mass_grid:", mass_grid)
    # splitted_mass_grid = lindexsplit(mass_grid, 153)

    # print(Z_list)
    # Z_extrapolation_list = [Z_list[0]*4, Z_list[0]*2] + Z_list
    Z_extrapolation_list = Z_list
    # Z_extrapolation_list = [Z_list[0] * 2] + Z_list + [Z_list[-1] / 10]
    # Z_extrapolation_list = list(reversed([Z_list[0]*2] + Z_list + [Z_list[-1]/10]))
    # Z_extrapolation_list = [Z_list[0]*4, Z_list[0]*2] + Z_list + [Z_list[-1]/10, Z_list[-1]/1E2, Z_list[-1]/1E3]
    # Z_extrapolation_list = [Z_list[0]*4, Z_list[0]*2] + Z_list + [Z_list[-1]/10, Z_list[-1]/1E2, Z_list[-1]/1E3, Z_list[-1]/1E4, Z_list[-1]/1E5]
    print("Z_list:", Z_extrapolation_list)


    for Z in range(len(Z_extrapolation_list)):
        metallicity = Z_extrapolation_list[Z]
        if metallicity > Z_list[0]:
            mass = M_list[0]
            eject_mass = eject_mass_list[0]
            lifetime = lifetime_list[0]
            Mfinal = Mfinal_list[0]
            H_eject_mass = H_eject_mass_list[0]
            He_eject_mass = He_eject_mass_list[0]
            C_eject_mass = time_list(C_eject_mass_list[0], 1, mass)
            C13_eject_mass = time_list(C13_eject_mass_list[0], metallicity / Z_list[0], mass)
            N_eject_mass = time_list_N(N_eject_mass_list[0], metallicity / Z_list[0], mass)
            O_eject_mass = time_list_alpha(O_eject_mass_list[0], metallicity / Z_list[0], mass)
            O17_eject_mass = time_list(O17_eject_mass_list[0], (metallicity / Z_list[0]) ** 2, mass)
            O18_eject_mass = time_list(O18_eject_mass_list[0], (metallicity / Z_list[0]) ** 2, mass)
            Ne_eject_mass = time_list_Fe(Ne_eject_mass_list[0], metallicity / Z_list[0], mass)
            Na_eject_mass = time_list_Fe(Na_eject_mass_list[0], metallicity / Z_list[0], mass)
            Mg_eject_mass = time_list_alpha(Mg_eject_mass_list[0], metallicity / Z_list[0], mass)
            Al_eject_mass = time_list_Fe(Al_eject_mass_list[0], metallicity / Z_list[0], mass)
            Si_eject_mass = time_list_alpha(Si_eject_mass_list[0], metallicity / Z_list[0], mass)
            S_eject_mass = time_list_alpha(S_eject_mass_list[0], metallicity / Z_list[0], mass)
            Ar_eject_mass = time_list_alpha(Ar_eject_mass_list[0], metallicity / Z_list[0], mass)
            Ca_eject_mass = time_list_alpha(Ca_eject_mass_list[0], metallicity / Z_list[0], mass)
            Ti_eject_mass = time_list(Ti_eject_mass_list[0], 1, mass)
            Fe_eject_mass = time_list_Fe(Fe_eject_mass_list[0], metallicity / Z_list[0], mass)
            Cr_eject_mass = time_list(Cr_eject_mass_list[0], 1, mass)
            Mn_eject_mass = time_list(Mn_eject_mass_list[0], 1, mass)
            Ni_eject_mass = time_list(Ni_eject_mass_list[0], metallicity / Z_list[0], mass)
            Metal_eject_mass = time_list(Metal_eject_mass_list[0], metallicity / Z_list[0], mass)
        elif metallicity < Z_list[-1]:
            mass = M_list[-1]
            eject_mass = eject_mass_list[-1]
            lifetime = lifetime_list[-1]
            Mfinal = Mfinal_list[-1]
            H_eject_mass = H_eject_mass_list[-1]
            He_eject_mass = He_eject_mass_list[-1]
            C_eject_mass = time_list(C_eject_mass_list[-1], 1, mass)
            C13_eject_mass = time_list(C13_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            N_eject_mass = time_list_N(N_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            O_eject_mass = time_list_alpha(O_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            O17_eject_mass = time_list(O17_eject_mass_list[-1], (metallicity / Z_list[-1]) ** 2, mass)
            O18_eject_mass = time_list(O18_eject_mass_list[-1], (metallicity / Z_list[-1]) ** 2, mass)
            Ne_eject_mass = time_list_Fe(Ne_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            Na_eject_mass = time_list_Fe(Na_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            Mg_eject_mass = time_list_alpha(Mg_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            Al_eject_mass = time_list_Fe(Al_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            Si_eject_mass = time_list_alpha(Si_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            S_eject_mass = time_list_alpha(S_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            Ar_eject_mass = time_list_alpha(Ar_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            Ca_eject_mass = time_list_alpha(Ca_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            Ti_eject_mass = time_list(Ti_eject_mass_list[-1], 1, mass)
            Fe_eject_mass = time_list_Fe(Fe_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            Cr_eject_mass = time_list(Cr_eject_mass_list[-1], 1, mass)
            Mn_eject_mass = time_list(Mn_eject_mass_list[-1], 1, mass)
            Ni_eject_mass = time_list(Ni_eject_mass_list[-1], metallicity / Z_list[-1], mass)
            Metal_eject_mass = time_list(Metal_eject_mass_list[-1], metallicity / Z_list[-1], mass)
        else:
            mass = M_list[Z]
            eject_mass = eject_mass_list[Z]
            lifetime = lifetime_list[Z]
            Mfinal = Mfinal_list[Z]
            H_eject_mass = H_eject_mass_list[Z]
            He_eject_mass = He_eject_mass_list[Z]
            C_eject_mass = C_eject_mass_list[Z]
            C13_eject_mass = C13_eject_mass_list[Z]
            N_eject_mass = N_eject_mass_list[Z]
            O_eject_mass = O_eject_mass_list[Z]
            O17_eject_mass = O17_eject_mass_list[Z]
            O18_eject_mass = O18_eject_mass_list[Z]
            Ne_eject_mass = Ne_eject_mass_list[Z]
            Na_eject_mass = Na_eject_mass_list[Z]
            Mg_eject_mass = Mg_eject_mass_list[Z]
            Al_eject_mass = Al_eject_mass_list[Z]
            Si_eject_mass = Si_eject_mass_list[Z]
            S_eject_mass = S_eject_mass_list[Z]
            Ar_eject_mass = Ar_eject_mass_list[Z]
            Ca_eject_mass = Ca_eject_mass_list[Z]
            Ti_eject_mass = Ti_eject_mass_list[Z]
            Fe_eject_mass = Fe_eject_mass_list[Z]
            Cr_eject_mass = Cr_eject_mass_list[Z]
            Mn_eject_mass = Mn_eject_mass_list[Z]
            Ni_eject_mass = Ni_eject_mass_list[Z]
            Metal_eject_mass = Metal_eject_mass_list[Z]

        if (yield_table_name == "Limongi_M000" or yield_table_name == "Limongi_M150" or yield_table_name == "Limongi_M300") and (metallicity < Z_list[-1] or metallicity == Z_list[-1]):
            for i in range(len(mass)):
                if mass[i] < 8.1:
                    Mfinal[i] = Mfinal_list[-2][i]

        ### Interpolate the metal yield ###

        # # portinari98 or marigo01:
        # eject_mass = np.interp(mass_grid, mass, eject_mass).tolist()
        # H_eject_mass = np.interp(mass_grid, mass, H_eject_mass).tolist()
        # He_eject_mass = np.interp(mass_grid, mass, He_eject_mass).tolist()
        # C_eject_mass = np.interp(mass_grid, mass, C_eject_mass).tolist()
        # N_eject_mass = np.interp(mass_grid, mass, N_eject_mass).tolist()
        # O_eject_mass = np.interp(mass_grid, mass, O_eject_mass).tolist()
        # Ne_eject_mass = np.interp(mass_grid, mass, Ne_eject_mass).tolist()
        # Na_eject_mass = np.interp(mass_grid, mass, Na_eject_mass).tolist()
        # Mg_eject_mass = np.interp(mass_grid, mass, Mg_eject_mass).tolist()
        # Si_eject_mass = np.interp(mass_grid, mass, Si_eject_mass).tolist()
        # S_eject_mass = np.interp(mass_grid, mass, S_eject_mass).tolist()
        # Ca_eject_mass = np.interp(mass_grid, mass, Ca_eject_mass).tolist()
        # Metal_eject_mass = np.interp(mass_grid, mass, Metal_eject_mass).tolist()
        # Fe_eject_mass = np.interp(mass_grid, mass, Fe_eject_mass).tolist()

        # # WW95
        # eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # H_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # He_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # C_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # N_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # O_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # Ne_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # Na_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # Mg_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # Si_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # S_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # Ca_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # Metal_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        # Fe_eject_mass_low = np.interp(splitted_mass_grid[0], [1], [0]).tolist()
        #
        # eject_mass_high = np.interp(splitted_mass_grid[1], mass, eject_mass).tolist()
        # H_eject_mass_high = np.interp(splitted_mass_grid[1], mass, H_eject_mass).tolist()
        # He_eject_mass_high = np.interp(splitted_mass_grid[1], mass, He_eject_mass).tolist()
        # C_eject_mass_high = np.interp(splitted_mass_grid[1], mass, C_eject_mass).tolist()
        # N_eject_mass_high = np.interp(splitted_mass_grid[1], mass, N_eject_mass).tolist()
        # O_eject_mass_high = np.interp(splitted_mass_grid[1], mass, O_eject_mass).tolist()
        # Ne_eject_mass_high = np.interp(splitted_mass_grid[1], mass, Ne_eject_mass).tolist()
        # Na_eject_mass_high = np.interp(splitted_mass_grid[1], mass, Na_eject_mass).tolist()
        # Mg_eject_mass_high = np.interp(splitted_mass_grid[1], mass, Mg_eject_mass).tolist()
        # Si_eject_mass_high = np.interp(splitted_mass_grid[1], mass, Si_eject_mass).tolist()
        # S_eject_mass_high = np.interp(splitted_mass_grid[1], mass, S_eject_mass).tolist()
        # Ca_eject_mass_high = np.interp(splitted_mass_grid[1], mass,Ca_eject_mass).tolist()
        # Metal_eject_mass_high = np.interp(splitted_mass_grid[1], mass, Metal_eject_mass).tolist()
        # Fe_eject_mass_high = np.interp(splitted_mass_grid[1], mass, Fe_eject_mass).tolist()
        #
        # eject_mass = eject_mass_low + eject_mass_high
        # H_eject_mass = H_eject_mass_low + H_eject_mass_high
        # He_eject_mass = He_eject_mass_low + He_eject_mass_high
        # C_eject_mass = C_eject_mass_low + C_eject_mass_high
        # N_eject_mass = N_eject_mass_low + N_eject_mass_high
        # O_eject_mass = O_eject_mass_low + O_eject_mass_high
        # Ne_eject_mass = Ne_eject_mass_low + Ne_eject_mass_high
        # Na_eject_mass = Na_eject_mass_low + Na_eject_mass_high
        # Mg_eject_mass = Mg_eject_mass_low + Mg_eject_mass_high
        # Si_eject_mass = Si_eject_mass_low + Si_eject_mass_high
        # S_eject_mass = S_eject_mass_low + S_eject_mass_high
        # Ca_eject_mass = Ca_eject_mass_low + Ca_eject_mass_high
        # Metal_eject_mass = Metal_eject_mass_low + Metal_eject_mass_high
        # Fe_eject_mass = Fe_eject_mass_low + Fe_eject_mass_high

        # ## Limongi 2018
        # for i in range(len(eject_mass)):
        #     if M_list[0][i] == 8:
        #         eject_mass[i] = eject_mass[i - 1] * M_list[0][i] / 6
        #         H_eject_mass[i] = H_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         He_eject_mass[i] = He_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         C_eject_mass[i] = C_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         C13_eject_mass[i] = C13_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         N_eject_mass[i] = N_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         O_eject_mass[i] = O_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         O17_eject_mass[i] = O17_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         O18_eject_mass[i] = O18_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Ne_eject_mass[i] = Ne_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Na_eject_mass[i] = Na_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Mg_eject_mass[i] = Mg_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Al_eject_mass[i] = Al_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Si_eject_mass[i] = Si_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         S_eject_mass[i] = S_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Ar_eject_mass[i] = Ar_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Ca_eject_mass[i] = Ca_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Ti_eject_mass[i] = Ti_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Cr_eject_mass[i] = Cr_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Mn_eject_mass[i] = Mn_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Ni_eject_mass[i] = Ni_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Metal_eject_mass[i] = Metal_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]
        #         Fe_eject_mass[i] = Fe_eject_mass[i] * eject_mass[i] / eject_mass[i - 1]

        # The interpretation of the tables are done in logarithmic scale, below first change to log, then change back.
        eject_mass__ = []
        lifetime__ = []
        Mfinal__ = []
        H_eject_mass__ = []
        He_eject_mass__ = []
        C_eject_mass__ = []
        C13_eject_mass__ = []
        N_eject_mass__ = []
        O_eject_mass__ = []
        O17_eject_mass__ = []
        O18_eject_mass__ = []
        Ne_eject_mass__ = []
        Na_eject_mass__ = []
        Mg_eject_mass__ = []
        Al_eject_mass__ = []
        Si_eject_mass__ = []
        S_eject_mass__ = []
        Ar_eject_mass__ = []
        Ca_eject_mass__ = []
        Ti_eject_mass__ = []
        Cr_eject_mass__ = []
        Mn_eject_mass__ = []
        Ni_eject_mass__ = []
        Metal_eject_mass__ = []
        Fe_eject_mass__ = []
        for i in range(len(eject_mass)):
            eject_mass__.append(math.log(eject_mass[i], 10))
            lifetime__.append(math.log(lifetime[i], 10))
            Mfinal__.append(math.log(Mfinal[i], 10))
            H_eject_mass__.append(math.log(H_eject_mass[i], 10))
            He_eject_mass__.append(math.log(He_eject_mass[i], 10))
            C_eject_mass__.append(math.log(C_eject_mass[i], 10))
            C13_eject_mass__.append(math.log(C13_eject_mass[i], 10))
            N_eject_mass__.append(math.log(N_eject_mass[i], 10))
            O_eject_mass__.append(math.log(O_eject_mass[i], 10))
            O17_eject_mass__.append(math.log(O17_eject_mass[i], 10))
            O18_eject_mass__.append(math.log(O18_eject_mass[i], 10))
            Ne_eject_mass__.append(math.log(Ne_eject_mass[i], 10))
            Na_eject_mass__.append(math.log(Na_eject_mass[i], 10))
            Mg_eject_mass__.append(math.log(Mg_eject_mass[i], 10))
            Al_eject_mass__.append(math.log(Al_eject_mass[i], 10))
            Si_eject_mass__.append(math.log(Si_eject_mass[i], 10))
            S_eject_mass__.append(math.log(S_eject_mass[i], 10))
            Ar_eject_mass__.append(math.log(Ar_eject_mass[i], 10))
            Ca_eject_mass__.append(math.log(Ca_eject_mass[i], 10))
            Ti_eject_mass__.append(math.log(Ti_eject_mass[i], 10))
            Cr_eject_mass__.append(math.log(Cr_eject_mass[i], 10))
            Mn_eject_mass__.append(math.log(Mn_eject_mass[i], 10))
            Ni_eject_mass__.append(math.log(Ni_eject_mass[i], 10))
            Metal_eject_mass__.append(math.log(Metal_eject_mass[i], 10))
            Fe_eject_mass__.append(math.log(Fe_eject_mass[i], 10))

        eject_mass = np.interp(mass_grid, mass, eject_mass__).tolist()
        lifetime = np.interp(mass_grid, mass, lifetime__).tolist()
        Mfinal = np.interp(mass_grid, mass, Mfinal__).tolist()
        H_eject_mass = np.interp(mass_grid, mass, H_eject_mass__).tolist()
        He_eject_mass = np.interp(mass_grid, mass, He_eject_mass__).tolist()
        C_eject_mass = np.interp(mass_grid, mass, C_eject_mass__).tolist()
        C13_eject_mass = np.interp(mass_grid, mass, C13_eject_mass__).tolist()
        N_eject_mass = np.interp(mass_grid, mass, N_eject_mass__).tolist()
        O_eject_mass = np.interp(mass_grid, mass, O_eject_mass__).tolist()
        O17_eject_mass = np.interp(mass_grid, mass, O17_eject_mass__).tolist()
        O18_eject_mass = np.interp(mass_grid, mass, O18_eject_mass__).tolist()
        Ne_eject_mass = np.interp(mass_grid, mass, Ne_eject_mass__).tolist()
        Na_eject_mass = np.interp(mass_grid, mass, Na_eject_mass__).tolist()
        Mg_eject_mass = np.interp(mass_grid, mass, Mg_eject_mass__).tolist()
        Al_eject_mass = np.interp(mass_grid, mass, Al_eject_mass__).tolist()
        Si_eject_mass = np.interp(mass_grid, mass, Si_eject_mass__).tolist()
        S_eject_mass = np.interp(mass_grid, mass, S_eject_mass__).tolist()
        Ar_eject_mass = np.interp(mass_grid, mass, Ar_eject_mass__).tolist()
        Ca_eject_mass = np.interp(mass_grid, mass, Ca_eject_mass__).tolist()
        Ti_eject_mass = np.interp(mass_grid, mass, Ti_eject_mass__).tolist()
        Cr_eject_mass = np.interp(mass_grid, mass, Cr_eject_mass__).tolist()
        Mn_eject_mass = np.interp(mass_grid, mass, Mn_eject_mass__).tolist()
        Ni_eject_mass = np.interp(mass_grid, mass, Ni_eject_mass__).tolist()
        Metal_eject_mass = np.interp(mass_grid, mass, Metal_eject_mass__).tolist()
        Fe_eject_mass = np.interp(mass_grid, mass, Fe_eject_mass__).tolist()

        for i in range(len(eject_mass)):
            eject_mass[i] = 10 ** eject_mass[i]
            lifetime[i] = 10 ** lifetime[i]
            Mfinal[i] = 10 ** Mfinal[i]
            H_eject_mass[i] = 10 ** H_eject_mass[i]
            He_eject_mass[i] = 10 ** He_eject_mass[i]
            C_eject_mass[i] = 10 ** C_eject_mass[i]
            C13_eject_mass[i] = 10 ** C13_eject_mass[i]
            N_eject_mass[i] = 10 ** N_eject_mass[i]
            O_eject_mass[i] = 10 ** O_eject_mass[i]
            O17_eject_mass[i] = 10 ** O17_eject_mass[i]
            O18_eject_mass[i] = 10 ** O18_eject_mass[i]
            Ne_eject_mass[i] = 10 ** Ne_eject_mass[i]
            Na_eject_mass[i] = 10 ** Na_eject_mass[i]
            Mg_eject_mass[i] = 10 ** Mg_eject_mass[i]
            Al_eject_mass[i] = 10 ** Al_eject_mass[i]
            Si_eject_mass[i] = 10 ** Si_eject_mass[i]
            S_eject_mass[i] = 10 ** S_eject_mass[i]
            Ar_eject_mass[i] = 10 ** Ar_eject_mass[i]
            Ca_eject_mass[i] = 10 ** Ca_eject_mass[i]
            Ti_eject_mass[i] = 10 ** Ti_eject_mass[i]
            Cr_eject_mass[i] = 10 ** Cr_eject_mass[i]
            Mn_eject_mass[i] = 10 ** Mn_eject_mass[i]
            Ni_eject_mass[i] = 10 ** Ni_eject_mass[i]
            Metal_eject_mass[i] = 10 ** Metal_eject_mass[i]
            Fe_eject_mass[i] = 10 ** Fe_eject_mass[i]

        for i in range(len(mass_grid)):
            ####################### Extrapolation for low-mass stars assuming the same yield #######################
            if mass_grid[i] < mass[0]:
                eject_mass[i] = eject_mass[i] * mass_grid[i] / mass[0]
                Mfinal[i] = Mfinal[i] * mass_grid[i] / mass[0]
                H_eject_mass[i] = H_eject_mass[i] * mass_grid[i] / mass[0]
                He_eject_mass[i] = He_eject_mass[i] * mass_grid[i] / mass[0]
                C_eject_mass[i] = C_eject_mass[i] * mass_grid[i] / mass[0]
                C13_eject_mass[i] = C13_eject_mass[i] * mass_grid[i] / mass[0]
                N_eject_mass[i] = N_eject_mass[i] * mass_grid[i] / mass[0]
                O_eject_mass[i] = O_eject_mass[i] * mass_grid[i] / mass[0]
                O17_eject_mass[i] = O17_eject_mass[i] * mass_grid[i] / mass[0]
                O18_eject_mass[i] = O18_eject_mass[i] * mass_grid[i] / mass[0]
                Ne_eject_mass[i] = Ne_eject_mass[i] * mass_grid[i] / mass[0]
                Na_eject_mass[i] = Na_eject_mass[i] * mass_grid[i] / mass[0]
                Mg_eject_mass[i] = Mg_eject_mass[i] * mass_grid[i] / mass[0]
                Al_eject_mass[i] = Al_eject_mass[i] * mass_grid[i] / mass[0]
                Si_eject_mass[i] = Si_eject_mass[i] * mass_grid[i] / mass[0]
                S_eject_mass[i] = S_eject_mass[i] * mass_grid[i] / mass[0]
                Ar_eject_mass[i] = Ar_eject_mass[i] * mass_grid[i] / mass[0]
                Ca_eject_mass[i] = Ca_eject_mass[i] * mass_grid[i] / mass[0]
                Ti_eject_mass[i] = Ti_eject_mass[i] * mass_grid[i] / mass[0]
                Cr_eject_mass[i] = Cr_eject_mass[i] * mass_grid[i] / mass[0]
                Mn_eject_mass[i] = Mn_eject_mass[i] * mass_grid[i] / mass[0]
                Ni_eject_mass[i] = Ni_eject_mass[i] * mass_grid[i] / mass[0]
                Metal_eject_mass[i] = Metal_eject_mass[i] * mass_grid[i] / mass[0]
                Fe_eject_mass[i] = Fe_eject_mass[i] * mass_grid[i] / mass[0]

            ### No extrapolation for massive stars, assuming the same yield in mass for metals,
            ### but the same yield in faction for H, He, and remnant.
            elif mass_grid[i] > mass[-1]:
                eject_mass[i] = eject_mass[i] * mass_grid[i] / mass[-1]
                Mfinal[i] = Mfinal[i] * mass_grid[i] / mass[-1]
                H_eject_mass[i] = H_eject_mass[i] * mass_grid[i] / mass[-1]
                He_eject_mass[i] = He_eject_mass[i] * mass_grid[i] / mass[-1]
            ### Extrapolation for massive stars, assuming the same yield in faction for metals:
            # elif mass_grid[i] > mass[-1]:
            #     eject_mass[i] = eject_mass[i] * mass_grid[i] / mass[-1]
            #     Mfinal[i] = Mfinal[i] * mass_grid[i] / mass[-1]
            #     H_eject_mass[i] = H_eject_mass[i] * mass_grid[i] / mass[-1]
            #     He_eject_mass[i] = He_eject_mass[i] * mass_grid[i] / mass[-1]
            #     C_eject_mass[i] = C_eject_mass[i] * mass_grid[i] / mass[-1]
            #     C13_eject_mass[i] = C13_eject_mass[i] * mass_grid[i] / mass[-1]
            #     N_eject_mass[i] = N_eject_mass[i] * mass_grid[i] / mass[-1]
            #     O_eject_mass[i] = O_eject_mass[i] * mass_grid[i] / mass[-1]
            #     O17_eject_mass[i] = O17_eject_mass[i] * mass_grid[i] / mass[-1]
            #     O18_eject_mass[i] = O18_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Ne_eject_mass[i] = Ne_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Na_eject_mass[i] = Na_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Mg_eject_mass[i] = Mg_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Al_eject_mass[i] = Al_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Si_eject_mass[i] = Si_eject_mass[i] * mass_grid[i] / mass[-1]
            #     S_eject_mass[i] = S_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Ar_eject_mass[i] = Ar_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Ca_eject_mass[i] = Ca_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Ti_eject_mass[i] = Ti_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Cr_eject_mass[i] = Cr_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Mn_eject_mass[i] = Mn_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Ni_eject_mass[i] = Ni_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Metal_eject_mass[i] = Metal_eject_mass[i] * mass_grid[i] / mass[-1]
            #     Fe_eject_mass[i] = Fe_eject_mass[i] * mass_grid[i] / mass[-1]

        # write file eject_mass
        out_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_eject_mass += '{} '.format(mass_grid[n])
        out_eject_mass += '\n# eject_mass\n'
        for n in range(len(eject_mass)):
            out_eject_mass += '{} '.format(eject_mass[n]/mass_grid[n])
        name_eject_mass = 'yield_tables_fractional/rearranged___/setllar_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name,
                                                                                                     yield_table_name,
                                                                                                     metallicity)
        file_eject_mass = open(name_eject_mass, 'w')
        file_eject_mass.write(out_eject_mass)
        file_eject_mass.close()

        # # write file lifetime. This is not written as we use the smooth spline fitted lifetime of Portinari98
        # out_lifetime = '# metallicity\n{}\n# mass\n'.format(metallicity)
        # for n in range(len(mass_grid)):
        #     out_lifetime += '{} '.format(mass_grid[n])
        # out_lifetime += '\n# lifetime\n'
        # for n in range(len(lifetime)):
        #     out_lifetime += '{} '.format(lifetime[n])
        # name_lifetime = 'yield_tables_fractional/rearranged___/setllar_lifetime_from_{}/{}_Z={}.txt'.format(yield_table_name,
        #                                                                                              yield_table_name,
        #                                                                                              metallicity)
        # file_lifetime = open(name_lifetime, 'w')
        # file_lifetime.write(out_lifetime)
        # file_lifetime.close()

        # write file Mfinal
        out_Mfinal = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Mfinal += '{} '.format(mass_grid[n])
        out_Mfinal += '\n# Mfinal\n'
        for n in range(len(Mfinal)):
            out_Mfinal += '{} '.format(Mfinal[n]/mass_grid[n])
        name_Mfinal = 'yield_tables_fractional/rearranged___/setllar_Mfinal_from_{}/{}_Z={}.txt'.format(yield_table_name,
                                                                                                     yield_table_name,
                                                                                                     metallicity)
        file_Mfinal = open(name_Mfinal, 'w')
        file_Mfinal.write(out_Mfinal)
        file_Mfinal.close()

        # write file H_eject_mass
        out_H_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_H_eject_mass += '{} '.format(mass_grid[n])
        out_H_eject_mass += '\n# H_eject_mass\n'
        for n in range(len(H_eject_mass)):
            out_H_eject_mass += '{} '.format(H_eject_mass[n]/mass_grid[n])
        # out_H_eject_mass += '\n# p30\n'
        # p30 = np.polyfit(mass_grid, H_eject_mass, 30)
        # for n in range(len(p30)):
        #     out_H_eject_mass += '{} '.format(p30[n])
        name_H_eject_mass = 'yield_tables_fractional/rearranged___/setllar_H_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_H_eject_mass = open(name_H_eject_mass, 'w')
        file_H_eject_mass.write(out_H_eject_mass)
        file_H_eject_mass.close()

        # write file He_eject_mass
        out_He_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_He_eject_mass += '{} '.format(mass_grid[n])
        out_He_eject_mass += '\n# He_eject_mass\n'
        for n in range(len(He_eject_mass)):
            out_He_eject_mass += '{} '.format(He_eject_mass[n]/mass_grid[n])
        # out_H_eject_mass += '\n# p30\n'
        # p30 = np.polyfit(mass_grid, He_eject_mass, 30)
        # for n in range(len(p30)):
        #     out_H_eject_mass += '{} '.format(p30[n])
        name_He_eject_mass = 'yield_tables_fractional/rearranged___/setllar_He_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_He_eject_mass = open(name_He_eject_mass, 'w')
        file_He_eject_mass.write(out_He_eject_mass)
        file_He_eject_mass.close()

        # write file C_eject_mass
        out_C_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_C_eject_mass += '{} '.format(mass_grid[n])
        out_C_eject_mass += '\n# C_eject_mass\n'
        for n in range(len(C_eject_mass)):
            out_C_eject_mass += '{} '.format(C_eject_mass[n]/mass_grid[n])
        name_C_eject_mass = 'yield_tables_fractional/rearranged___/setllar_C_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_C_eject_mass = open(name_C_eject_mass, 'w')
        file_C_eject_mass.write(out_C_eject_mass)
        file_C_eject_mass.close()

        # write file C13_eject_mass
        out_C13_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_C13_eject_mass += '{} '.format(mass_grid[n])
        out_C13_eject_mass += '\n# C13_eject_mass\n'
        for n in range(len(C13_eject_mass)):
            out_C13_eject_mass += '{} '.format(C13_eject_mass[n]/mass_grid[n])
        name_C13_eject_mass = 'yield_tables_fractional/rearranged___/setllar_C13_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_C13_eject_mass = open(name_C13_eject_mass, 'w')
        file_C13_eject_mass.write(out_C13_eject_mass)
        file_C13_eject_mass.close()

        # write file N_eject_mass
        out_N_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_N_eject_mass += '{} '.format(mass_grid[n])
        out_N_eject_mass += '\n# N_eject_mass\n'
        for n in range(len(N_eject_mass)):
            out_N_eject_mass += '{} '.format(N_eject_mass[n]/mass_grid[n])
        name_N_eject_mass = 'yield_tables_fractional/rearranged___/setllar_N_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_N_eject_mass = open(name_N_eject_mass, 'w')
        file_N_eject_mass.write(out_N_eject_mass)
        file_N_eject_mass.close()

        # write file O_eject_mass
        out_O_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_O_eject_mass += '{} '.format(mass_grid[n])
        out_O_eject_mass += '\n# O_eject_mass\n'
        for n in range(len(O_eject_mass)):
            out_O_eject_mass += '{} '.format(O_eject_mass[n]/mass_grid[n])
        name_O_eject_mass = 'yield_tables_fractional/rearranged___/setllar_O_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_O_eject_mass = open(name_O_eject_mass, 'w')
        file_O_eject_mass.write(out_O_eject_mass)
        file_O_eject_mass.close()

        # write file O17_eject_mass
        out_O17_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_O17_eject_mass += '{} '.format(mass_grid[n])
        out_O17_eject_mass += '\n# O17_eject_mass\n'
        for n in range(len(O17_eject_mass)):
            out_O17_eject_mass += '{} '.format(O17_eject_mass[n]/mass_grid[n])
        name_O17_eject_mass = 'yield_tables_fractional/rearranged___/setllar_O17_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_O17_eject_mass = open(name_O17_eject_mass, 'w')
        file_O17_eject_mass.write(out_O17_eject_mass)
        file_O17_eject_mass.close()

        # write file O18_eject_mass
        out_O18_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_O18_eject_mass += '{} '.format(mass_grid[n])
        out_O18_eject_mass += '\n# O18_eject_mass\n'
        for n in range(len(O18_eject_mass)):
            out_O18_eject_mass += '{} '.format(O18_eject_mass[n]/mass_grid[n])
        name_O18_eject_mass = 'yield_tables_fractional/rearranged___/setllar_O18_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_O18_eject_mass = open(name_O18_eject_mass, 'w')
        file_O18_eject_mass.write(out_O18_eject_mass)
        file_O18_eject_mass.close()

        # write file Ne_eject_mass
        out_Ne_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Ne_eject_mass += '{} '.format(mass_grid[n])
        out_Ne_eject_mass += '\n# Ne_eject_mass\n'
        for n in range(len(Ne_eject_mass)):
            out_Ne_eject_mass += '{} '.format(Ne_eject_mass[n]/mass_grid[n])
        name_Ne_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Ne_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Ne_eject_mass = open(name_Ne_eject_mass, 'w')
        file_Ne_eject_mass.write(out_Ne_eject_mass)
        file_Ne_eject_mass.close()

        # write file Na_eject_mass
        out_Na_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Na_eject_mass += '{} '.format(mass_grid[n])
        out_Na_eject_mass += '\n# Na_eject_mass\n'
        for n in range(len(Na_eject_mass)):
            out_Na_eject_mass += '{} '.format(Na_eject_mass[n]/mass_grid[n])
        name_Na_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Na_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Na_eject_mass = open(name_Na_eject_mass, 'w')
        file_Na_eject_mass.write(out_Na_eject_mass)
        file_Na_eject_mass.close()

        # write file Mg_eject_mass
        out_Mg_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Mg_eject_mass += '{} '.format(mass_grid[n])
        out_Mg_eject_mass += '\n# Mg_eject_mass\n'
        for n in range(len(Mg_eject_mass)):
            out_Mg_eject_mass += '{} '.format(Mg_eject_mass[n]/mass_grid[n])
        name_Mg_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Mg_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Mg_eject_mass = open(name_Mg_eject_mass, 'w')
        file_Mg_eject_mass.write(out_Mg_eject_mass)
        file_Mg_eject_mass.close()

        # write file Al_eject_mass
        out_Al_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Al_eject_mass += '{} '.format(mass_grid[n])
        out_Al_eject_mass += '\n# Al_eject_mass\n'
        for n in range(len(Al_eject_mass)):
            out_Al_eject_mass += '{} '.format(Al_eject_mass[n]/mass_grid[n])
        name_Al_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Al_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Al_eject_mass = open(name_Al_eject_mass, 'w')
        file_Al_eject_mass.write(out_Al_eject_mass)
        file_Al_eject_mass.close()

        # write file Si_eject_mass
        out_Si_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Si_eject_mass += '{} '.format(mass_grid[n])
        out_Si_eject_mass += '\n# Si_eject_mass\n'
        for n in range(len(Si_eject_mass)):
            out_Si_eject_mass += '{} '.format(Si_eject_mass[n]/mass_grid[n])
        name_Si_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Si_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Si_eject_mass = open(name_Si_eject_mass, 'w')
        file_Si_eject_mass.write(out_Si_eject_mass)
        file_Si_eject_mass.close()

        # write file S_eject_mass
        out_S_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_S_eject_mass += '{} '.format(mass_grid[n])
        out_S_eject_mass += '\n# S_eject_mass\n'
        for n in range(len(S_eject_mass)):
            out_S_eject_mass += '{} '.format(S_eject_mass[n]/mass_grid[n])
        name_S_eject_mass = 'yield_tables_fractional/rearranged___/setllar_S_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_S_eject_mass = open(name_S_eject_mass, 'w')
        file_S_eject_mass.write(out_S_eject_mass)
        file_S_eject_mass.close()

        # write file Ar_eject_mass
        out_Ar_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Ar_eject_mass += '{} '.format(mass_grid[n])
        out_Ar_eject_mass += '\n# Ar_eject_mass\n'
        for n in range(len(Ar_eject_mass)):
            out_Ar_eject_mass += '{} '.format(Ar_eject_mass[n]/mass_grid[n])
        name_Ar_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Ar_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Ar_eject_mass = open(name_Ar_eject_mass, 'w')
        file_Ar_eject_mass.write(out_Ar_eject_mass)
        file_Ar_eject_mass.close()

        # write file Ca_eject_mass
        out_Ca_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Ca_eject_mass += '{} '.format(mass_grid[n])
        out_Ca_eject_mass += '\n# Ca_eject_mass\n'
        for n in range(len(Ca_eject_mass)):
            out_Ca_eject_mass += '{} '.format(Ca_eject_mass[n]/mass_grid[n])
        name_Ca_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Ca_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Ca_eject_mass = open(name_Ca_eject_mass, 'w')
        file_Ca_eject_mass.write(out_Ca_eject_mass)
        file_Ca_eject_mass.close()

        # write file Ti_eject_mass
        out_Ti_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Ti_eject_mass += '{} '.format(mass_grid[n])
        out_Ti_eject_mass += '\n# Ti_eject_mass\n'
        for n in range(len(Ti_eject_mass)):
            out_Ti_eject_mass += '{} '.format(Ti_eject_mass[n]/mass_grid[n])
        name_Ti_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Ti_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Ti_eject_mass = open(name_Ti_eject_mass, 'w')
        file_Ti_eject_mass.write(out_Ti_eject_mass)
        file_Ti_eject_mass.close()

        # write file Cr_eject_mass
        out_Cr_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Cr_eject_mass += '{} '.format(mass_grid[n])
        out_Cr_eject_mass += '\n# Cr_eject_mass\n'
        for n in range(len(Cr_eject_mass)):
            out_Cr_eject_mass += '{} '.format(Cr_eject_mass[n]/mass_grid[n])
        name_Cr_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Cr_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Cr_eject_mass = open(name_Cr_eject_mass, 'w')
        file_Cr_eject_mass.write(out_Cr_eject_mass)
        file_Cr_eject_mass.close()

        # write file Mn_eject_mass
        out_Mn_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Mn_eject_mass += '{} '.format(mass_grid[n])
        out_Mn_eject_mass += '\n# Mn_eject_mass\n'
        for n in range(len(Mn_eject_mass)):
            out_Mn_eject_mass += '{} '.format(Mn_eject_mass[n]/mass_grid[n])
        name_Mn_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Mn_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Mn_eject_mass = open(name_Mn_eject_mass, 'w')
        file_Mn_eject_mass.write(out_Mn_eject_mass)
        file_Mn_eject_mass.close()

        # write file Ni_eject_mass
        out_Ni_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Ni_eject_mass += '{} '.format(mass_grid[n])
        out_Ni_eject_mass += '\n# Ni_eject_mass\n'
        for n in range(len(Ni_eject_mass)):
            out_Ni_eject_mass += '{} '.format(Ni_eject_mass[n]/mass_grid[n])
        name_Ni_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Ni_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Ni_eject_mass = open(name_Ni_eject_mass, 'w')
        file_Ni_eject_mass.write(out_Ni_eject_mass)
        file_Ni_eject_mass.close()

        # write file Fe_eject_mass
        out_Fe_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Fe_eject_mass += '{} '.format(mass_grid[n])
        out_Fe_eject_mass += '\n# Fe_eject_mass\n'
        for n in range(len(Fe_eject_mass)):
            out_Fe_eject_mass += '{} '.format(Fe_eject_mass[n]/mass_grid[n])
        name_Fe_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Fe_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Fe_eject_mass = open(name_Fe_eject_mass, 'w')
        file_Fe_eject_mass.write(out_Fe_eject_mass)
        file_Fe_eject_mass.close()

        # write file Metal_eject_mass
        out_Metal_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Metal_eject_mass += '{} '.format(mass_grid[n])
        out_Metal_eject_mass += '\n# Metal_eject_mass\n'
        for n in range(len(Metal_eject_mass)):
            out_Metal_eject_mass += '{} '.format(Metal_eject_mass[n]/mass_grid[n])
        name_Metal_eject_mass = 'yield_tables_fractional/rearranged___/setllar_Metal_eject_mass_from_{}/{}_Z={}.txt'.format(
            yield_table_name, yield_table_name, metallicity)
        file_Metal_eject_mass = open(name_Metal_eject_mass, 'w')
        file_Metal_eject_mass.write(out_Metal_eject_mass)
        file_Metal_eject_mass.close()


    print('yield_tables_fractional/rearranged___/setllar_...eject_mass_from_{}/{}_Z=....txt saved'.format(yield_table_name,
                                                                                               yield_table_name))
    return


def function_get_Mfinal_and_Lifetime(input_string):
    i_end = len(input_string)
    i = 0
    in_str = ''
    while i < i_end:
        in_str += input_string[i]
        (i) = (i + 1)
    output = float(in_str)
    return output


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
            while j < 300:
                line_j = str.split(data[j])
                if line_j[0] == '&' + element:
                    end = j
                    element_relative_line_number = j - i
                    break
                (j) = (j + 1)
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
            i_M_start = i + 2
        if M_Z_string[i] == ',':
            i_M_end = i
            i_Z_start = i + 3
        if M_Z_string[i] == ')':
            i_Z_end = i
        (i) = (i + 1)
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
    global O_over_Mg_list, C13_over_O17_list, C_over_O_list, C13_over_C_list, O17_over_O_list, Al_over_Mg_list, \
        Al_over_O_list, Mg_over_Fe_list, Al_over_Fe_list, Ne_over_H_list, Na_over_H_list, Mg_over_H_list, Fe_over_H_list, Fe_eject_mass_list, \
        C_over_H_list, N_over_H_list, O_over_H_list, C13_over_H_list, O17_over_H_list, O18_over_H_list, Si_over_H_list, S_over_H_list, Ar_over_H_list, Ca_over_H_list, \
        Cr_over_H_list, Mn_over_H_list, Ni_over_H_list, Ti_over_H_list, Z_over_H_list, O_over_Fe_list, Si_over_Fe_list, \
        Ca_over_Fe_list, Cr_over_Fe_list, Mn_over_Fe_list, Ni_over_Fe_list, Ti_over_Fe_list, M_list, Z_list, Al_over_H_list
    j = 0
    while j < len(M_list):
        i = 0
        while i < len(M_list[j]):
            M_list[j][i] = math.log(M_list[j][i], 10)
            (i) = (i + 1)
        (j) = (j + 1)

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(30, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    i = 0
    while i < len(M_list):
        plt.plot([10 ** x__ for x__ in M_list[i]], Fe_eject_mass_list[i], color=colors[i],
                 label='Z={}'.format(Z_list[i]))
        # spline = UnivariateSpline(M_list[i], Fe_eject_mass_list[i], s=0, k=5)
        # xs = np.linspace(M_list[i][-1], M_list[i][0], 1000)
        # ys = spline(xs)
        # plt.plot(xs, ys, ls="dashed")

        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'Fe', yield_table_name, yield_table_name, Z_list[i]), 'r')
        data = file_M_eject.readlines()
        stellar_mass = [float(x) for x in data[3].split()]
        M_eject_txt = data[5]
        file_M_eject.close()
        M_eject_table = [float(x) for x in M_eject_txt.split()]
        plt.scatter(stellar_mass, M_eject_table)
        # spline = UnivariateSpline(mass_grid, M_eject_table, s=0, k=1)
        # xs = np.linspace(mass_grid[0], mass_grid[-1], 1000)
        # ys = spline(xs)
        # plt.plot(xs, ys, ls="dashed")
        (i) = (i + 1)
    plt.xscale("log")
    # plt.yscale("log")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'Fe eject')
    plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(31, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    i = 0
    while i < len(M_list):
        plt.plot([10 ** x__ for x__ in M_list[i]], O_eject_mass_list[i], color=colors[i],
                 label='Z={}'.format(Z_list[i]))
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'O', yield_table_name, yield_table_name, Z_list[i]), 'r')
        data = file_M_eject.readlines()
        # stellar_mass = [float(x) for x in data[3].split()]
        M_eject_txt = data[5]
        file_M_eject.close()
        mass_grid = function_get_mass_grid(yield_table_name)
        M_eject_table = [float(x) for x in M_eject_txt.split()]
        plt.scatter(mass_grid, M_eject_table)
        (i) = (i + 1)
    plt.xscale("log")
    # plt.yscale("log")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'O eject')
    plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(32, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    i = 0
    while i < len(M_list):
        plt.plot([10 ** x__ for x__ in M_list[i]], C_eject_mass_list[i], color=colors[i],
                 label='Z={}'.format(Z_list[i]))
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'C', yield_table_name, yield_table_name, Z_list[i]), 'r')
        data = file_M_eject.readlines()
        # stellar_mass = [float(x) for x in data[3].split()]
        M_eject_txt = data[5]
        file_M_eject.close()
        mass_grid = function_get_mass_grid(yield_table_name)
        M_eject_table = [float(x) for x in M_eject_txt.split()]
        plt.scatter(mass_grid, M_eject_table)
        (i) = (i + 1)
    plt.xscale("log")
    # plt.yscale("log")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'C eject')
    plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(33, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    i = 0
    while i < len(M_list):
        plt.plot([10 ** x__ for x__ in M_list[i]], N_eject_mass_list[i], color=colors[i],
                 label='Z={}'.format(Z_list[i]))
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'N', yield_table_name, yield_table_name, Z_list[i]), 'r')
        data = file_M_eject.readlines()
        # stellar_mass = [float(x) for x in data[3].split()]
        M_eject_txt = data[5]
        file_M_eject.close()
        mass_grid = function_get_mass_grid(yield_table_name)
        M_eject_table = [float(x) for x in M_eject_txt.split()]
        plt.scatter(mass_grid, M_eject_table)
        (i) = (i + 1)
    plt.xscale("log")
    # plt.yscale("log")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'N eject')
    plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(34, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    i = 0
    while i < len(M_list):
        plt.plot([10 ** x__ for x__ in M_list[i]], Mg_eject_mass_list[i], color=colors[i],
                 label='Z={}'.format(Z_list[i]))
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'Mg', yield_table_name, yield_table_name, Z_list[i]), 'r')
        data = file_M_eject.readlines()
        # stellar_mass = [float(x) for x in data[3].split()]
        M_eject_txt = data[5]
        file_M_eject.close()
        mass_grid = function_get_mass_grid(yield_table_name)
        M_eject_table = [float(x) for x in M_eject_txt.split()]
        plt.scatter(mass_grid, M_eject_table)
        (i) = (i + 1)
    plt.xscale("log")
    # plt.yscale("log")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'Mg eject')
    plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(35, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    i = 0
    while i < len(Z_list):
        # plt.plot([10 ** x__ for x__ in M_list[i]], Ca_eject_mass_list[i], color=colors[i],
        #          label='Z={}'.format(Z_list[i]))
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'Ca', yield_table_name, yield_table_name, Z_list[i]), 'r')
        Ca_data = file_M_eject.readlines()
        file_M_eject.close()
        Ca_eject_txt = Ca_data[5]
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'Fe', yield_table_name, yield_table_name, Z_list[i]), 'r')
        Fe_data = file_M_eject.readlines()
        file_M_eject.close()
        Fe_eject_txt = Fe_data[5]
        # stellar_mass = [float(x) for x in data[3].split()]
        mass_grid = function_get_mass_grid(yield_table_name)
        Ca_eject_table = [float(x) for x in Ca_eject_txt.split()]
        Fe_eject_table = [float(x) for x in Fe_eject_txt.split()]
        plt.plot(mass_grid, [math.log(Ca_eject_table[j]/element_weight_table.function_element_weight("Ca")/Fe_eject_table[j]*element_weight_table.function_element_weight("Fe"), 10) - solar_Ca + solar_Fe for j in range(len(Ca_eject_table))], label='Z={}'.format(Z_list[i]), color=colors[i])
        # plt.plot([10**x for x in M_list[i]], Ca_over_Fe_list[i], c='k', lw=4)
        if Z_list[i] == 0.0003:
            plt.plot(mass_grid, [math.log(
                Ca_eject_table[j] / element_weight_table.function_element_weight("Ca") / Fe_eject_table[
                    j] * element_weight_table.function_element_weight("Fe"), 10) - solar_Ca + solar_Fe for j in
                                 range(len(Ca_eject_table))], label='Z={}'.format(Z_list[i]),
                     color="b", lw=2)
        if Z_list[i] == 0.0001:
            plt.plot(mass_grid, [math.log(
                Ca_eject_table[j] / element_weight_table.function_element_weight("Ca") / Fe_eject_table[
                    j] * element_weight_table.function_element_weight("Fe"), 10) - solar_Ca + solar_Fe for j in
                                 range(len(Ca_eject_table))], label='Z={}'.format(Z_list[i]),
                     color="r", lw=2)
        (i) = (i + 1)
    plt.xscale("log")
    # plt.yscale("log")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'Stellar mass [M$_\odot$]')
    plt.ylabel(r'[Ca/Fe]')
    plt.tight_layout()



    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(36, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    i = 0
    while i < len(Z_list):
        # plt.plot([10 ** x__ for x__ in M_list[i]], Si_eject_mass_list[i], color=colors[i],
        #          label='Z={}'.format(Z_list[i]))
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'Si', yield_table_name, yield_table_name, Z_list[i]), 'r')
        Si_data = file_M_eject.readlines()
        file_M_eject.close()
        Si_eject_txt = Si_data[5]
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'Fe', yield_table_name, yield_table_name, Z_list[i]), 'r')
        Fe_data = file_M_eject.readlines()
        file_M_eject.close()
        Fe_eject_txt = Fe_data[5]
        # stellar_mass = [float(x) for x in data[3].split()]
        mass_grid = function_get_mass_grid(yield_table_name)
        Si_eject_table = [float(x) for x in Si_eject_txt.split()]
        Fe_eject_table = [float(x) for x in Fe_eject_txt.split()]
        plt.plot(mass_grid, [math.log(Si_eject_table[j]/element_weight_table.function_element_weight("Si")/Fe_eject_table[j]*element_weight_table.function_element_weight("Fe"), 10) - solar_Si + solar_Fe for j in range(len(Si_eject_table))], label='Z={}'.format(Z_list[i]), color=colors[i])
        # plt.plot([10**x for x in M_list[i]], Si_over_Fe_list[i], c='k', lw=4)
        if Z_list[i] == 0.0003:
            plt.plot(mass_grid, [math.log(
                Si_eject_table[j] / element_weight_table.function_element_weight("Si") / Fe_eject_table[
                    j] * element_weight_table.function_element_weight("Fe"), 10) - solar_Si + solar_Fe for j in
                                 range(len(Si_eject_table))], label='Z={}'.format(Z_list[i]),
                     color="b", lw=2)
        if Z_list[i] == 0.0001:
            plt.plot(mass_grid, [math.log(
                Si_eject_table[j] / element_weight_table.function_element_weight("Si") / Fe_eject_table[
                    j] * element_weight_table.function_element_weight("Fe"), 10) - solar_Si + solar_Fe for j in
                                 range(len(Si_eject_table))], label='Z={}'.format(Z_list[i]),
                     color="r", lw=2)
        (i) = (i + 1)
    plt.xscale("log")
    # plt.yscale("log")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'Stellar mass [M$_\odot$]')
    plt.ylabel(r'[Si/Fe]')
    plt.tight_layout()

    # plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(36, figsize=(6, 5.25))
    # ax = fig.add_subplot(1, 1, 1)
    # i = 0
    # while i < len(M_list):
    #     plt.plot([10 ** x__ for x__ in M_list[i]], Si_eject_mass_list[i], color=colors[i],
    #              label='Z={}'.format(Z_list[i]))
    #     yield_path = 'yield_tables_fractional'
    #     file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
    #         'Si', yield_table_name, yield_table_name, Z_list[i]), 'r')
    #     data = file_M_eject.readlines()
    #     # stellar_mass = [float(x) for x in data[3].split()]
    #     M_eject_txt = data[5]
    #     file_M_eject.close()
    #     mass_grid = function_get_mass_grid(yield_table_name)
    #     M_eject_table = [float(x) for x in M_eject_txt.split()]
    #     plt.scatter(mass_grid, M_eject_table)
    #     (i) = (i + 1)
    # plt.xscale("log")
    # # plt.yscale("log")
    # plt.legend(prop={'size': 10}, loc='best')
    # plt.xlabel(r'log stellar mass [M$_\odot$]')
    # plt.ylabel(r'Si eject')
    # plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(37, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    i = 0
    while i < len(Z_list):
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'O', yield_table_name, yield_table_name, Z_list[i]), 'r')
        O_data = file_M_eject.readlines()
        file_M_eject.close()
        O_eject_txt = O_data[5]
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'Fe', yield_table_name, yield_table_name, Z_list[i]), 'r')
        Fe_data = file_M_eject.readlines()
        file_M_eject.close()
        Fe_eject_txt = Fe_data[5]
        # stellar_mass = [float(x) for x in data[3].split()]
        mass_grid = function_get_mass_grid(yield_table_name)
        O_eject_table = [float(x) for x in O_eject_txt.split()]
        Fe_eject_table = [float(x) for x in Fe_eject_txt.split()]
        plt.plot(mass_grid, [math.log(
            O_eject_table[j] / element_weight_table.function_element_weight("O") / Fe_eject_table[
                j] * element_weight_table.function_element_weight("Fe"), 10) - solar_O + solar_Fe for j in
                             range(len(O_eject_table))], label='Z={}'.format(Z_list[i]))
                             # range(len(Ca_eject_table))], label='Z={}'.format(Z_list[i]), color=colors[i])
        # if Z_list[i] == 0.0003:
        #     plt.plot(mass_grid, [math.log(
        #         Ca_eject_table[j] / element_weight_table.function_element_weight("O") / Fe_eject_table[
        #             j] * element_weight_table.function_element_weight("Fe"), 10) - solar_O + solar_Fe for j in
        #                          range(len(Ca_eject_table))], label='Z={}'.format(Z_list[i]),
        #              color="b", lw=2)
        (i) = (i + 1)
    plt.xscale("log")
    # plt.yscale("log")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'Stellar mass [M$_\odot$]')
    plt.ylabel(r'[O/Fe]')
    plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(38, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    i = 0
    while i < len(Z_list):
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'Al', yield_table_name, yield_table_name, Z_list[i]), 'r')
        Ca_data = file_M_eject.readlines()
        file_M_eject.close()
        Ca_eject_txt = Ca_data[5]
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'Fe', yield_table_name, yield_table_name, Z_list[i]), 'r')
        Fe_data = file_M_eject.readlines()
        file_M_eject.close()
        Fe_eject_txt = Fe_data[5]
        # stellar_mass = [float(x) for x in data[3].split()]
        mass_grid = function_get_mass_grid(yield_table_name)
        Ca_eject_table = [float(x) for x in Ca_eject_txt.split()]
        Fe_eject_table = [float(x) for x in Fe_eject_txt.split()]
        plt.plot(mass_grid, [math.log(
            Ca_eject_table[j] / element_weight_table.function_element_weight("Al") / Fe_eject_table[
                j] * element_weight_table.function_element_weight("Fe"), 10) - solar_Al + solar_Fe for j in
                             range(len(Ca_eject_table))], label='Z={}'.format(Z_list[i]), color=colors[i])
        (i) = (i + 1)
    plt.xscale("log")
    # plt.yscale("log")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'Stellar mass [M$_\odot$]')
    plt.ylabel(r'[Al/Fe]')
    plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(39, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    i = 0
    while i < len(Z_list):
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'O', yield_table_name, yield_table_name, Z_list[i]), 'r')
        O_data = file_M_eject.readlines()
        file_M_eject.close()
        O_eject_txt = O_data[5]
        yield_path = 'yield_tables_fractional'
        file_M_eject = open(yield_path + '/rearranged___/setllar_{}_eject_mass_from_{}/{}_Z={}.txt'.format(
            'N', yield_table_name, yield_table_name, Z_list[i]), 'r')
        N_data = file_M_eject.readlines()
        file_M_eject.close()
        N_eject_txt = N_data[5]
        # stellar_mass = [float(x) for x in data[3].split()]
        mass_grid = function_get_mass_grid(yield_table_name)
        O_eject_table = [float(x) for x in O_eject_txt.split()]
        N_eject_table = [float(x) for x in N_eject_txt.split()]
        plt.plot(mass_grid, [math.log(N_eject_table[j] / element_weight_table.function_element_weight("N") / O_eject_table[j] * element_weight_table.function_element_weight("O"), 10)
                                 for j in range(len(O_eject_table))], label='Z={}'.format(Z_list[i]), color=colors[i])
        # plt.plot(mass_grid, [math.log(N_eject_table[j] / element_weight_table.function_element_weight("N") / O_eject_table[j] * element_weight_table.function_element_weight("O"), 10) - solar_N + solar_O
        #                          for j in range(len(O_eject_table))], label='Z={}'.format(Z_list[i]), color=colors[i])
        (i) = (i + 1)
    # plt.xscale("log")
    # plt.yscale("log")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'Stellar mass [M$_\odot$]')
    plt.ylabel(r'N/O')
    plt.tight_layout()

    plt.show()
    return


if __name__ == '__main__':
    start_time = time.time()
    Z_list = []
    M_list = []
    eject_mass_list = []
    lifetime_list = []
    Mfinal_list = []
    H_eject_mass_list = []
    He_eject_mass_list = []
    C_eject_mass_list = []
    C13_eject_mass_list = []
    N_eject_mass_list = []
    O_eject_mass_list = []
    O17_eject_mass_list = []
    O18_eject_mass_list = []
    Ne_eject_mass_list = []
    Na_eject_mass_list = []
    Mg_eject_mass_list = []
    Al_eject_mass_list = []
    Si_eject_mass_list = []
    S_eject_mass_list = []
    Ar_eject_mass_list = []
    Ca_eject_mass_list = []
    Ti_eject_mass_list = []
    Cr_eject_mass_list = []
    Mn_eject_mass_list = []
    Fe_eject_mass_list = []
    Ni_eject_mass_list = []
    Metal_eject_mass_list = []
    O_over_Mg_list = []
    C13_over_O17_list = []
    C_over_O_list = []
    C13_over_C_list = []
    O17_over_O_list = []
    Al_over_Mg_list = []
    Al_over_O_list = []
    Ne_over_H_list = []
    Na_over_H_list = []
    Mg_over_H_list = []
    Al_over_H_list = []
    Si_over_H_list = []
    S_over_H_list = []
    Ar_over_H_list = []
    Ca_over_H_list = []
    Ti_over_H_list = []
    Cr_over_H_list = []
    Mn_over_H_list = []
    Fe_over_H_list = []
    Ni_over_H_list = []
    C_over_H_list = []
    N_over_H_list = []
    O_over_H_list = []
    C13_over_H_list = []
    O17_over_H_list = []
    O18_over_H_list = []
    Z_over_H_list = []
    O_over_Fe_list = []
    Mg_over_Fe_list = []
    Al_over_Fe_list = []
    Si_over_Fe_list = []
    Ca_over_Fe_list = []
    Ti_over_Fe_list = []
    Cr_over_Fe_list = []
    Mn_over_Fe_list = []
    Ni_over_Fe_list = []
    yield_table_name = "Limongi_M000"  # "K10_K06_HNe10" or "C15_N13_HNe10" or "WW95" or "portinari98" or "marigo01" or "Kobayashi06" or "Karakas10"
    # or "Nomoto" or "Nomoto_HNe" or "Nomoto_ZY_hypernova" or "Nomoto_ZY_CCSN" or "Nomoto_ZY_CCSN_popIII" or "K10_N13_HegerPopIII"
    # or "Limongi_R000" or "Limongi_R300" or "Limongi_R150" or "K10_N13_HNe00"
    # or "popIII_N13" or "popIII_heger10"
    # or "nugrid_N13" or "nugrid_FRUITY" or "nugrid_MESAonly_ye" or "nugrid_K06"
    make_dir = True  # mkdir for new yield tables.
    # make_dir = False
    function_read_file(yield_table_name)
    funtion_plot_yields()

    # yield_table_name = "WW95" or "portinari98" or "marigo01"
    # or marigo01+"Kobayashi06" or "Karakas10"+Kobayashi06 or Cristallo15+"Nomoto"13 or Karakas10+"Limongi_R000"18

    # plt.show()
    print(" - Run time: %s -" % round((time.time() - start_time), 2))