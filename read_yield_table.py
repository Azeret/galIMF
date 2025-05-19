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

SNIa_yield_table = "Iwamoto1999_WDD3" # "Iwamoto1999_W70"-0.39 #"Seitenzahl2013"-0.16 "Iwamoto1999"-0.24
C_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'C')
O_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'O')
Mg_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Mg')
Al_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Al')
Si_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Si')
S_mass_eject_SNIa = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'S')
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
solar_N = solar.function_solar_element_abundances("Asplund2009", "N")
solar_O = solar.function_solar_element_abundances("Asplund2009", "O")
solar_Ne = solar.function_solar_element_abundances("Asplund2009", "Ne")
solar_Mg = solar.function_solar_element_abundances("Asplund2009", "Mg")
solar_Al = solar.function_solar_element_abundances("Asplund2009", "Al")
solar_Si = solar.function_solar_element_abundances("Asplund2009", "Si")
solar_S = solar.function_solar_element_abundances("Asplund2009", "S")
solar_Ca = solar.function_solar_element_abundances("Asplund2009", "Ca")
solar_Ti = solar.function_solar_element_abundances("Asplund2009", "Ti")
solar_Cr = solar.function_solar_element_abundances("Asplund2009", "Cr")
solar_Mn = solar.function_solar_element_abundances("Asplund2009", "Mn")
solar_Fe = solar.function_solar_element_abundances("Asplund2009", "Fe")
solar_Ni = solar.function_solar_element_abundances("Asplund2009", "Ni")
solar_Na = solar.function_solar_element_abundances("Asplund2009", "Na")

def function_read_file(yield_table_name):

    ####################
    ### read in file ###
    ####################
    if yield_table_name == "portinari98":
        file_yield = open(
            'yield_tables/agb_and_massive_stars_portinari98_marigo01_gce_totalyields.txt', 'r')
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
    elif yield_table_name == "Limongi_R000" or yield_table_name == "Limongi" :
        file_yield = open('yield_tables/agb_and_massive_stars_K10_LC18_R000.txt', 'r')
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
    Na_relative_line_number = function_get_element_line_number(data, 'Na-23')
    Mg_relative_line_number = function_get_element_line_number(data, 'Mg-24')
    Al_relative_line_number = function_get_element_line_number(data, 'Al-27')
    Si_relative_line_number = function_get_element_line_number(data, 'Si-28')
    S_relative_line_number = function_get_element_line_number(data, 'S-32')
    Ca_relative_line_number = function_get_element_line_number(data, 'Ca-40')
    Ti_relative_line_number = function_get_element_line_number(data, 'Ti-48')
    Cr_relative_line_number = function_get_element_line_number(data, 'Cr-52')
    Mn_relative_line_number = function_get_element_line_number(data, 'Mn-55')
    Fe_relative_line_number = function_get_element_line_number(data, 'Fe-56')
    Ni_relative_line_number = function_get_element_line_number(data, 'Ni-58')
    #
    global M_list, Z_list, eject_mass_list, H_eject_mass_list, He_eject_mass_list, C_eject_mass_list, \
        N_eject_mass_list, O_eject_mass_list, Ne_eject_mass_list, Na_eject_mass_list, Mg_eject_mass_list, Al_eject_mass_list, Si_eject_mass_list, \
        S_eject_mass_list, Ca_eject_mass_list, Ti_eject_mass_list, Cr_eject_mass_list, Mn_eject_mass_list, Fe_eject_mass_list, Ni_eject_mass_list, Metal_eject_mass_list
    global O_over_Mg_list, Mg_over_Fe_list, Mg_over_H_list, Fe_over_H_list, O_over_H_list, Z_over_H_list, O_over_Fe_list
    #
    i = 0
    while i < len(data):
        line_i = str.split(data[i])
        if line_i[1] == 'Table:':  # Here select the lines being: ['H', 'Table:', '(M=xxx,Z=xxx)']
            line_H = str.split(data[i + H_relative_line_number])
            line_He = str.split(data[i + He_relative_line_number])
            line_C = str.split(data[i + C_relative_line_number])
            line_N = str.split(data[i + N_relative_line_number])
            line_O = str.split(data[i + O_relative_line_number])
            line_Ne = str.split(data[i + Ne_relative_line_number])
            line_Na = str.split(data[i + Na_relative_line_number])
            line_Mg = str.split(data[i + Mg_relative_line_number])
            line_Al = str.split(data[i + Al_relative_line_number])
            line_Si = str.split(data[i + Si_relative_line_number])
            line_S = str.split(data[i + S_relative_line_number])
            line_Ca = str.split(data[i + Ca_relative_line_number])
            line_Ti = str.split(data[i + Ti_relative_line_number])
            line_Cr = str.split(data[i + Cr_relative_line_number])
            line_Mn = str.split(data[i + Mn_relative_line_number])
            line_Fe = str.split(data[i + Fe_relative_line_number])
            line_Ni = str.split(data[i + Ni_relative_line_number])
            line_Lifetime = str.split(data[i + 1])
            Lifetime = function_get_Mfinal_and_Lifetime(line_Lifetime[2])
            line_Mfinal = str.split(data[i + 2])
            Mfinal = function_get_Mfinal_and_Lifetime(line_Mfinal[2])
            (Z_ini, M_ini) = function_get_Z_M(line_i[2])  # get the initial mass and metallicity of the star
            ejecta_mass = round((M_ini - Mfinal), 5)  ####################
            H_mass = function_get_element_mass(line_H[1])
            He_mass = function_get_element_mass(line_He[1])
            C_mass = function_get_element_mass(line_C[1])
            N_mass = function_get_element_mass(line_N[1])
            O_mass = function_get_element_mass(line_O[1])
            Ne_mass = function_get_element_mass(line_Ne[1])
            Na_mass = function_get_element_mass(line_Na[1])
            Mg_mass = function_get_element_mass(line_Mg[1])
            Al_mass = function_get_element_mass(line_Al[1])
            Si_mass = function_get_element_mass(line_Si[1])
            S_mass = function_get_element_mass(line_S[1])
            Ca_mass = function_get_element_mass(line_Ca[1])
            Ti_mass = function_get_element_mass(line_Ti[1])
            Cr_mass = function_get_element_mass(line_Cr[1])
            Mn_mass = function_get_element_mass(line_Mn[1])
            Fe_mass = function_get_element_mass(line_Fe[1])
            Ni_mass = function_get_element_mass(line_Ni[1])
            if Si_mass == 0:
                Si_mass = 1e-30
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
            H_num = H_mass / element_weight_table.function_element_weight("H")
            C_num = C_mass / element_weight_table.function_element_weight("C")
            N_num = N_mass / element_weight_table.function_element_weight("N")
            O_num = O_mass / element_weight_table.function_element_weight("O")
            Na_num = Na_mass / element_weight_table.function_element_weight("Na")
            Mg_num = Mg_mass / element_weight_table.function_element_weight("Mg")
            Al_num = Al_mass / element_weight_table.function_element_weight("Al")
            Si_num = Si_mass / element_weight_table.function_element_weight("Si")
            Ca_num = Ca_mass / element_weight_table.function_element_weight("Ca")
            Ti_num = Ti_mass / element_weight_table.function_element_weight("Ti")
            Cr_num = Cr_mass / element_weight_table.function_element_weight("Cr")
            Mn_num = Mn_mass / element_weight_table.function_element_weight("Mn")
            Fe_num = Fe_mass / element_weight_table.function_element_weight("Fe")
            Ni_num = Ni_mass / element_weight_table.function_element_weight("Ni")
            O_over_Mg = math.log(O_num/Mg_num, 10) - solar_O + solar_Mg
            Mg_over_H = math.log(Mg_num/H_num, 10) - solar_Mg + solar_H
            Fe_over_H = math.log(Fe_num/H_num, 10) - solar_Fe + solar_H
            O_over_H = math.log(O_num/H_num, 10) - solar_O + solar_H
            Mg_over_Fe = math.log(Mg_num/Fe_num, 10) - solar_Mg + solar_Fe
            O_over_Fe = math.log(O_num/Fe_num, 10) - solar_O + solar_Fe
            Metal_mass = round((ejecta_mass - H_mass - He_mass), 5)  ####################
            if Metal_mass<0 or Metal_mass==0:
                print("Warning: Metal_mass=", Metal_mass, "<0")
                print("check stellar yield table with metallicity and mass being:", Z_ini, "&", M_ini)
                Metal_mass = 1e-6
            Z_over_H = math.log(Metal_mass / H_mass, 10) - math.log(0.0134 / 0.7381, 10)
            if len(Z_list) == 0:
                Z_list.append(Z_ini)
                Z_n = 0
                M_list.append([])
                eject_mass_list.append([])
                H_eject_mass_list.append([])
                He_eject_mass_list.append([])
                C_eject_mass_list.append([])
                N_eject_mass_list.append([])
                O_eject_mass_list.append([])
                Ne_eject_mass_list.append([])
                Na_eject_mass_list.append([])
                Mg_eject_mass_list.append([])
                Al_eject_mass_list.append([])
                Si_eject_mass_list.append([])
                S_eject_mass_list.append([])
                Ca_eject_mass_list.append([])
                Ti_eject_mass_list.append([])
                Cr_eject_mass_list.append([])
                Mn_eject_mass_list.append([])
                Fe_eject_mass_list.append([])
                Ni_eject_mass_list.append([])
                Metal_eject_mass_list.append([])
                Z_over_H_list.append([])
                O_over_Mg_list.append([])
                Mg_over_Fe_list.append([])
                Mg_over_H_list.append([])
                Fe_over_H_list.append([])
                O_over_H_list.append([])
                O_over_Fe_list.append([])
            if Z_ini != Z_list[-1]:
                Z_list.append(Z_ini)
                Z_n += 1
                M_list.append([])
                eject_mass_list.append([])
                H_eject_mass_list.append([])
                He_eject_mass_list.append([])
                C_eject_mass_list.append([])
                N_eject_mass_list.append([])
                O_eject_mass_list.append([])
                Ne_eject_mass_list.append([])
                Na_eject_mass_list.append([])
                Mg_eject_mass_list.append([])
                Al_eject_mass_list.append([])
                Si_eject_mass_list.append([])
                S_eject_mass_list.append([])
                Ca_eject_mass_list.append([])
                Ti_eject_mass_list.append([])
                Cr_eject_mass_list.append([])
                Mn_eject_mass_list.append([])
                Fe_eject_mass_list.append([])
                Ni_eject_mass_list.append([])
                Metal_eject_mass_list.append([])
                O_over_Mg_list.append([])
                Mg_over_Fe_list.append([])
                Mg_over_H_list.append([])
                Fe_over_H_list.append([])
                O_over_H_list.append([])
                Z_over_H_list.append([])
                O_over_Fe_list.append([])
            M_list[Z_n].append(M_ini)
            eject_mass_list[Z_n].append(ejecta_mass)
            H_eject_mass_list[Z_n].append(H_mass)
            He_eject_mass_list[Z_n].append(He_mass)
            C_eject_mass_list[Z_n].append(C_mass)
            N_eject_mass_list[Z_n].append(N_mass)
            O_eject_mass_list[Z_n].append(O_mass)
            Ne_eject_mass_list[Z_n].append(Ne_mass)
            Na_eject_mass_list[Z_n].append(Na_mass)
            Mg_eject_mass_list[Z_n].append(Mg_mass)
            Al_eject_mass_list[Z_n].append(Al_mass)
            Si_eject_mass_list[Z_n].append(Si_mass)
            S_eject_mass_list[Z_n].append(S_mass)
            Ca_eject_mass_list[Z_n].append(Ca_mass)
            Ti_eject_mass_list[Z_n].append(Ti_mass)
            Cr_eject_mass_list[Z_n].append(Cr_mass)
            Mn_eject_mass_list[Z_n].append(Mn_mass)
            Fe_eject_mass_list[Z_n].append(Fe_mass)
            Ni_eject_mass_list[Z_n].append(Ni_mass)
            Metal_eject_mass_list[Z_n].append(Metal_mass)
            O_over_Mg_list[Z_n].append(O_over_Mg)
            Mg_over_Fe_list[Z_n].append(Mg_over_Fe)
            Mg_over_H_list[Z_n].append(Mg_over_H)
            O_over_H_list[Z_n].append(O_over_H)
            Z_over_H_list[Z_n].append(Z_over_H)
            Fe_over_H_list[Z_n].append(Fe_over_H)
            O_over_Fe_list[Z_n].append(O_over_Fe)
        (i) = (i + 1)
    # print(Z_list)
    # print(M_list)
    # print(eject_mass_list)
    # print(H_eject_mass_list)
    # print(Fe_eject_mass_list)

    ###########################
    ### write data to files ###
    ###########################

    write_data()

    return


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

def function_get_mass_grid(yield_table_name): # read in a grid from 0.08 to 150 Msun
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
    elif yield_table_name == "Limongi_R000" or yield_table_name == "Limongi":
        file_yield = open('yield_tables/agb_and_massive_stars_K10_LC18_R000.txt', 'r')
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

def write_data():
    global M_list, Z_list, eject_mass_list, H_eject_mass_list, He_eject_mass_list, C_eject_mass_list, \
        N_eject_mass_list, O_eject_mass_list, Ne_eject_mass_list, Na_eject_mass_list, Mg_eject_mass_list, Al_eject_mass_list, Si_eject_mass_list, \
        S_eject_mass_list, Ca_eject_mass_list, Ti_eject_mass_list, Cr_eject_mass_list, Mn_eject_mass_list, Fe_eject_mass_list, Ni_eject_mass_list, Metal_eject_mass_list

    os.makedirs('yield_tables/rearranged___/setllar_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_H_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_He_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_C_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_N_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_O_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Ne_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Na_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Mg_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Al_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Si_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_S_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Ca_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Ti_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Cr_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Mn_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Fe_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Ni_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)
    os.makedirs('yield_tables/rearranged___/setllar_Metal_eject_mass_from_{}'.format(yield_table_name), exist_ok=True)

    mass_grid = function_get_mass_grid(yield_table_name)

    for Z in range(len(Z_list)):
        metallicity = Z_list[Z]
        mass = M_list[Z]
        eject_mass = eject_mass_list[Z]
        H_eject_mass = H_eject_mass_list[Z]
        He_eject_mass = He_eject_mass_list[Z]
        C_eject_mass = C_eject_mass_list[Z]
        N_eject_mass = N_eject_mass_list[Z]
        O_eject_mass = O_eject_mass_list[Z]
        Ne_eject_mass = Ne_eject_mass_list[Z]
        Na_eject_mass = Na_eject_mass_list[Z]
        Mg_eject_mass = Mg_eject_mass_list[Z]
        Al_eject_mass = Al_eject_mass_list[Z]
        Si_eject_mass = Si_eject_mass_list[Z]
        S_eject_mass = S_eject_mass_list[Z]
        Ca_eject_mass = Ca_eject_mass_list[Z]
        Ti_eject_mass = Ti_eject_mass_list[Z]
        Fe_eject_mass = Fe_eject_mass_list[Z]
        Cr_eject_mass = Cr_eject_mass_list[Z]
        Mn_eject_mass = Mn_eject_mass_list[Z]
        Ni_eject_mass = Ni_eject_mass_list[Z]
        Metal_eject_mass = Metal_eject_mass_list[Z]

        ### Interpolate the metal yield ###
        # The interpretation of the tables are done in logarithmic scale, below first change to log, then change back.
        eject_mass__ = []
        H_eject_mass__ = []
        He_eject_mass__ = []
        C_eject_mass__ = []
        N_eject_mass__ = []
        O_eject_mass__ = []
        Ne_eject_mass__ = []
        Na_eject_mass__ = []
        Mg_eject_mass__ = []
        Al_eject_mass__ = []
        Si_eject_mass__ = []
        S_eject_mass__ = []
        Ca_eject_mass__ = []
        Ti_eject_mass__ = []
        Cr_eject_mass__ = []
        Mn_eject_mass__ = []
        Ni_eject_mass__ = []
        Metal_eject_mass__ = []
        Fe_eject_mass__ = []
        for i in range(len(eject_mass)):
            eject_mass__.append(math.log(eject_mass[i], 10))
            H_eject_mass__.append(math.log(H_eject_mass[i], 10))
            He_eject_mass__.append(math.log(He_eject_mass[i], 10))
            C_eject_mass__.append(math.log(C_eject_mass[i], 10))
            N_eject_mass__.append(math.log(N_eject_mass[i], 10))
            O_eject_mass__.append(math.log(O_eject_mass[i], 10))
            Ne_eject_mass__.append(math.log(Ne_eject_mass[i], 10))
            Na_eject_mass__.append(math.log(Na_eject_mass[i], 10))
            Mg_eject_mass__.append(math.log(Mg_eject_mass[i], 10))
            Al_eject_mass__.append(math.log(Al_eject_mass[i], 10))
            Si_eject_mass__.append(math.log(Si_eject_mass[i], 10))
            S_eject_mass__.append(math.log(S_eject_mass[i], 10))
            Ca_eject_mass__.append(math.log(Ca_eject_mass[i], 10))
            Ti_eject_mass__.append(math.log(Ti_eject_mass[i], 10))
            Cr_eject_mass__.append(math.log(Cr_eject_mass[i], 10))
            Mn_eject_mass__.append(math.log(Mn_eject_mass[i], 10))
            Ni_eject_mass__.append(math.log(Ni_eject_mass[i], 10))
            Metal_eject_mass__.append(math.log(Metal_eject_mass[i], 10))
            Fe_eject_mass__.append(math.log(Fe_eject_mass[i], 10))

        eject_mass = np.interp(mass_grid, mass, eject_mass__).tolist()
        H_eject_mass = np.interp(mass_grid, mass, H_eject_mass__).tolist()
        He_eject_mass = np.interp(mass_grid, mass, He_eject_mass__).tolist()
        C_eject_mass = np.interp(mass_grid, mass, C_eject_mass__).tolist()
        N_eject_mass = np.interp(mass_grid, mass, N_eject_mass__).tolist()
        O_eject_mass = np.interp(mass_grid, mass, O_eject_mass__).tolist()
        Ne_eject_mass = np.interp(mass_grid, mass, Ne_eject_mass__).tolist()
        Na_eject_mass = np.interp(mass_grid, mass, Na_eject_mass__).tolist()
        Mg_eject_mass = np.interp(mass_grid, mass, Mg_eject_mass__).tolist()
        Al_eject_mass = np.interp(mass_grid, mass, Al_eject_mass__).tolist()
        Si_eject_mass = np.interp(mass_grid, mass, Si_eject_mass__).tolist()
        S_eject_mass = np.interp(mass_grid, mass, S_eject_mass__).tolist()
        Ca_eject_mass = np.interp(mass_grid, mass, Ca_eject_mass__).tolist()
        Ti_eject_mass = np.interp(mass_grid, mass, Ti_eject_mass__).tolist()
        Cr_eject_mass = np.interp(mass_grid, mass, Cr_eject_mass__).tolist()
        Mn_eject_mass = np.interp(mass_grid, mass, Mn_eject_mass__).tolist()
        Ni_eject_mass = np.interp(mass_grid, mass, Ni_eject_mass__).tolist()
        Metal_eject_mass = np.interp(mass_grid, mass, Metal_eject_mass__).tolist()
        Fe_eject_mass = np.interp(mass_grid, mass, Fe_eject_mass__).tolist()

        for i in range(len(eject_mass)):
            eject_mass[i] = 10**eject_mass[i]
            H_eject_mass[i] = 10**H_eject_mass[i]
            He_eject_mass[i] = 10**He_eject_mass[i]
            C_eject_mass[i] = 10**C_eject_mass[i]
            N_eject_mass[i] = 10**N_eject_mass[i]
            O_eject_mass[i] = 10**O_eject_mass[i]
            Ne_eject_mass[i] = 10**Ne_eject_mass[i]
            Na_eject_mass[i] = 10**Na_eject_mass[i]
            Mg_eject_mass[i] = 10**Mg_eject_mass[i]
            Al_eject_mass[i] = 10**Al_eject_mass[i]
            Si_eject_mass[i] = 10**Si_eject_mass[i]
            S_eject_mass[i] = 10**S_eject_mass[i]
            Ca_eject_mass[i] = 10**Ca_eject_mass[i]
            Ti_eject_mass[i] = 10**Ti_eject_mass[i]
            Cr_eject_mass[i] = 10**Cr_eject_mass[i]
            Mn_eject_mass[i] = 10**Mn_eject_mass[i]
            Ni_eject_mass[i] = 10**Ni_eject_mass[i]
            Metal_eject_mass[i] = 10**Metal_eject_mass[i]
            Fe_eject_mass[i] = 10**Fe_eject_mass[i]

        for i in range(len(mass_grid)):
            ####################### Extrapolation for low-mass stars assuming the same yield #######################
            if mass_grid[i] < mass[0]:
                eject_mass[i] = eject_mass[i] * mass_grid[i]/mass[0]
                H_eject_mass[i] = H_eject_mass[i] * mass_grid[i]/mass[0]
                He_eject_mass[i] = He_eject_mass[i] * mass_grid[i]/mass[0]
                C_eject_mass[i] = C_eject_mass[i] * mass_grid[i]/mass[0]
                N_eject_mass[i] = N_eject_mass[i] * mass_grid[i]/mass[0]
                O_eject_mass[i] = O_eject_mass[i] * mass_grid[i]/mass[0]
                Ne_eject_mass[i] = Ne_eject_mass[i] * mass_grid[i]/mass[0]
                Na_eject_mass[i] = Na_eject_mass[i] * mass_grid[i]/mass[0]
                Mg_eject_mass[i] = Mg_eject_mass[i] * mass_grid[i]/mass[0]
                Al_eject_mass[i] = Al_eject_mass[i] * mass_grid[i]/mass[0]
                Si_eject_mass[i] = Si_eject_mass[i] * mass_grid[i]/mass[0]
                S_eject_mass[i] = S_eject_mass[i] * mass_grid[i]/mass[0]
                Ca_eject_mass[i] = Ca_eject_mass[i] * mass_grid[i]/mass[0]
                Ti_eject_mass[i] = Ti_eject_mass[i] * mass_grid[i]/mass[0]
                Cr_eject_mass[i] = Cr_eject_mass[i] * mass_grid[i]/mass[0]
                Mn_eject_mass[i] = Mn_eject_mass[i] * mass_grid[i]/mass[0]
                Ni_eject_mass[i] = Ni_eject_mass[i] * mass_grid[i]/mass[0]
                Metal_eject_mass[i] = Metal_eject_mass[i] * mass_grid[i]/mass[0]
                Fe_eject_mass[i] = Fe_eject_mass[i] * mass_grid[i]/mass[0]
            ####################### Extrapolation for massive stars assuming the same yield #######################
            elif mass_grid[i] > mass[-1]:
                eject_mass[i] = eject_mass[i] * mass_grid[i] / mass[-1]
                H_eject_mass[i] = H_eject_mass[i] * mass_grid[i] / mass[-1]
                He_eject_mass[i] = He_eject_mass[i] * mass_grid[i] / mass[-1]
                C_eject_mass[i] = C_eject_mass[i] * mass_grid[i] / mass[-1]
                N_eject_mass[i] = N_eject_mass[i] * mass_grid[i] / mass[-1]
                O_eject_mass[i] = O_eject_mass[i] * mass_grid[i] / mass[-1]
                Ne_eject_mass[i] = Ne_eject_mass[i] * mass_grid[i] / mass[-1]
                Na_eject_mass[i] = Na_eject_mass[i] * mass_grid[i] / mass[-1]
                Mg_eject_mass[i] = Mg_eject_mass[i] * mass_grid[i] / mass[-1]
                Al_eject_mass[i] = Al_eject_mass[i] * mass_grid[i] / mass[-1]
                Si_eject_mass[i] = Si_eject_mass[i] * mass_grid[i] / mass[-1]
                S_eject_mass[i] = S_eject_mass[i] * mass_grid[i] / mass[-1]
                Ca_eject_mass[i] = Ca_eject_mass[i] * mass_grid[i] / mass[-1]
                Ti_eject_mass[i] = Ti_eject_mass[i] * mass_grid[i] / mass[-1]
                Cr_eject_mass[i] = Cr_eject_mass[i] * mass_grid[i] / mass[-1]
                Mn_eject_mass[i] = Mn_eject_mass[i] * mass_grid[i] / mass[-1]
                Ni_eject_mass[i] = Ni_eject_mass[i] * mass_grid[i] / mass[-1]
                Metal_eject_mass[i] = Metal_eject_mass[i] * mass_grid[i] / mass[-1]
                Fe_eject_mass[i] = Fe_eject_mass[i] * mass_grid[i] / mass[-1]


        # write file eject_mass
        out_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_eject_mass += '{} '.format(mass_grid[n])
        out_eject_mass += '\n# eject_mass\n'
        for n in range(len(eject_mass)):
            out_eject_mass += '{} '.format(eject_mass[n])
        name_eject_mass = 'yield_tables/rearranged___/setllar_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_eject_mass = open(name_eject_mass, 'w')
        file_eject_mass.write(out_eject_mass)
        file_eject_mass.close()

        # write file H_eject_mass
        out_H_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_H_eject_mass += '{} '.format(mass_grid[n])
        out_H_eject_mass += '\n# H_eject_mass\n'
        for n in range(len(H_eject_mass)):
            out_H_eject_mass += '{} '.format(H_eject_mass[n])
        name_H_eject_mass = 'yield_tables/rearranged___/setllar_H_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_H_eject_mass = open(name_H_eject_mass, 'w')
        file_H_eject_mass.write(out_H_eject_mass)
        file_H_eject_mass.close()

        # write file He_eject_mass
        out_He_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_He_eject_mass += '{} '.format(mass_grid[n])
        out_He_eject_mass += '\n# He_eject_mass\n'
        for n in range(len(He_eject_mass)):
            out_He_eject_mass += '{} '.format(He_eject_mass[n])
        name_He_eject_mass = 'yield_tables/rearranged___/setllar_He_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_He_eject_mass = open(name_He_eject_mass, 'w')
        file_He_eject_mass.write(out_He_eject_mass)
        file_He_eject_mass.close()

        # write file C_eject_mass
        out_C_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_C_eject_mass += '{} '.format(mass_grid[n])
        out_C_eject_mass += '\n# C_eject_mass\n'
        for n in range(len(C_eject_mass)):
            out_C_eject_mass += '{} '.format(C_eject_mass[n])
        name_C_eject_mass = 'yield_tables/rearranged___/setllar_C_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_C_eject_mass = open(name_C_eject_mass, 'w')
        file_C_eject_mass.write(out_C_eject_mass)
        file_C_eject_mass.close()

        # write file N_eject_mass
        out_N_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_N_eject_mass += '{} '.format(mass_grid[n])
        out_N_eject_mass += '\n# N_eject_mass\n'
        for n in range(len(N_eject_mass)):
            out_N_eject_mass += '{} '.format(N_eject_mass[n])
        name_N_eject_mass = 'yield_tables/rearranged___/setllar_N_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_N_eject_mass = open(name_N_eject_mass, 'w')
        file_N_eject_mass.write(out_N_eject_mass)
        file_N_eject_mass.close()

        # write file O_eject_mass
        out_O_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_O_eject_mass += '{} '.format(mass_grid[n])
        out_O_eject_mass += '\n# O_eject_mass\n'
        for n in range(len(O_eject_mass)):
            out_O_eject_mass += '{} '.format(O_eject_mass[n])
        name_O_eject_mass = 'yield_tables/rearranged___/setllar_O_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_O_eject_mass = open(name_O_eject_mass, 'w')
        file_O_eject_mass.write(out_O_eject_mass)
        file_O_eject_mass.close()

        # write file Ne_eject_mass
        out_Ne_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Ne_eject_mass += '{} '.format(mass_grid[n])
        out_Ne_eject_mass += '\n# Ne_eject_mass\n'
        for n in range(len(Ne_eject_mass)):
            out_Ne_eject_mass += '{} '.format(Ne_eject_mass[n])
        name_Ne_eject_mass = 'yield_tables/rearranged___/setllar_Ne_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Ne_eject_mass = open(name_Ne_eject_mass, 'w')
        file_Ne_eject_mass.write(out_Ne_eject_mass)
        file_Ne_eject_mass.close()

        # write file Na_eject_mass
        out_Na_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Na_eject_mass += '{} '.format(mass_grid[n])
        out_Na_eject_mass += '\n# Na_eject_mass\n'
        for n in range(len(Na_eject_mass)):
            out_Na_eject_mass += '{} '.format(Na_eject_mass[n])
        name_Na_eject_mass = 'yield_tables/rearranged___/setllar_Na_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Na_eject_mass = open(name_Na_eject_mass, 'w')
        file_Na_eject_mass.write(out_Na_eject_mass)
        file_Na_eject_mass.close()

        # write file Mg_eject_mass
        out_Mg_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Mg_eject_mass += '{} '.format(mass_grid[n])
        out_Mg_eject_mass += '\n# Mg_eject_mass\n'
        for n in range(len(Mg_eject_mass)):
            out_Mg_eject_mass += '{} '.format(Mg_eject_mass[n])
        name_Mg_eject_mass = 'yield_tables/rearranged___/setllar_Mg_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Mg_eject_mass = open(name_Mg_eject_mass, 'w')
        file_Mg_eject_mass.write(out_Mg_eject_mass)
        file_Mg_eject_mass.close()

        # write file Al_eject_mass
        out_Al_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Al_eject_mass += '{} '.format(mass_grid[n])
        out_Al_eject_mass += '\n# Al_eject_mass\n'
        for n in range(len(Al_eject_mass)):
            out_Al_eject_mass += '{} '.format(Al_eject_mass[n])
        name_Al_eject_mass = 'yield_tables/rearranged___/setllar_Al_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Al_eject_mass = open(name_Al_eject_mass, 'w')
        file_Al_eject_mass.write(out_Al_eject_mass)
        file_Al_eject_mass.close()

        # write file Si_eject_mass
        out_Si_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Si_eject_mass += '{} '.format(mass_grid[n])
        out_Si_eject_mass += '\n# Si_eject_mass\n'
        for n in range(len(Si_eject_mass)):
            out_Si_eject_mass += '{} '.format(Si_eject_mass[n])
        name_Si_eject_mass = 'yield_tables/rearranged___/setllar_Si_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Si_eject_mass = open(name_Si_eject_mass, 'w')
        file_Si_eject_mass.write(out_Si_eject_mass)
        file_Si_eject_mass.close()

        # write file S_eject_mass
        out_S_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_S_eject_mass += '{} '.format(mass_grid[n])
        out_S_eject_mass += '\n# S_eject_mass\n'
        for n in range(len(S_eject_mass)):
            out_S_eject_mass += '{} '.format(S_eject_mass[n])
        name_S_eject_mass = 'yield_tables/rearranged___/setllar_S_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_S_eject_mass = open(name_S_eject_mass, 'w')
        file_S_eject_mass.write(out_S_eject_mass)
        file_S_eject_mass.close()

        # write file Ca_eject_mass
        out_Ca_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Ca_eject_mass += '{} '.format(mass_grid[n])
        out_Ca_eject_mass += '\n# Ca_eject_mass\n'
        for n in range(len(Ca_eject_mass)):
            out_Ca_eject_mass += '{} '.format(Ca_eject_mass[n])
        name_Ca_eject_mass = 'yield_tables/rearranged___/setllar_Ca_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Ca_eject_mass = open(name_Ca_eject_mass, 'w')
        file_Ca_eject_mass.write(out_Ca_eject_mass)
        file_Ca_eject_mass.close()


        # write file Ti_eject_mass
        out_Ti_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Ti_eject_mass += '{} '.format(mass_grid[n])
        out_Ti_eject_mass += '\n# Ti_eject_mass\n'
        for n in range(len(Ti_eject_mass)):
            out_Ti_eject_mass += '{} '.format(Ti_eject_mass[n])
        name_Ti_eject_mass = 'yield_tables/rearranged___/setllar_Ti_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Ti_eject_mass = open(name_Ti_eject_mass, 'w')
        file_Ti_eject_mass.write(out_Ti_eject_mass)
        file_Ti_eject_mass.close()

        # write file Cr_eject_mass
        out_Cr_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Cr_eject_mass += '{} '.format(mass_grid[n])
        out_Cr_eject_mass += '\n# Cr_eject_mass\n'
        for n in range(len(Cr_eject_mass)):
            out_Cr_eject_mass += '{} '.format(Cr_eject_mass[n])
        name_Cr_eject_mass = 'yield_tables/rearranged___/setllar_Cr_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Cr_eject_mass = open(name_Cr_eject_mass, 'w')
        file_Cr_eject_mass.write(out_Cr_eject_mass)
        file_Cr_eject_mass.close()

        # write file Mn_eject_mass
        out_Mn_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Mn_eject_mass += '{} '.format(mass_grid[n])
        out_Mn_eject_mass += '\n# Mn_eject_mass\n'
        for n in range(len(Mn_eject_mass)):
            out_Mn_eject_mass += '{} '.format(Mn_eject_mass[n])
        name_Mn_eject_mass = 'yield_tables/rearranged___/setllar_Mn_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Mn_eject_mass = open(name_Mn_eject_mass, 'w')
        file_Mn_eject_mass.write(out_Mn_eject_mass)
        file_Mn_eject_mass.close()

        # write file Ni_eject_mass
        out_Ni_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Ni_eject_mass += '{} '.format(mass_grid[n])
        out_Ni_eject_mass += '\n# Ni_eject_mass\n'
        for n in range(len(Ni_eject_mass)):
            out_Ni_eject_mass += '{} '.format(Ni_eject_mass[n])
        name_Ni_eject_mass = 'yield_tables/rearranged___/setllar_Ni_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Ni_eject_mass = open(name_Ni_eject_mass, 'w')
        file_Ni_eject_mass.write(out_Ni_eject_mass)
        file_Ni_eject_mass.close()


        # write file Fe_eject_mass
        out_Fe_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Fe_eject_mass += '{} '.format(mass_grid[n])
        out_Fe_eject_mass += '\n# Fe_eject_mass\n'
        for n in range(len(Fe_eject_mass)):
            out_Fe_eject_mass += '{} '.format(Fe_eject_mass[n])
        name_Fe_eject_mass = 'yield_tables/rearranged___/setllar_Fe_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Fe_eject_mass = open(name_Fe_eject_mass, 'w')
        file_Fe_eject_mass.write(out_Fe_eject_mass)
        file_Fe_eject_mass.close()

        # write file Metal_eject_mass
        out_Metal_eject_mass = '# metallicity\n{}\n# mass\n'.format(metallicity)
        for n in range(len(mass_grid)):
            out_Metal_eject_mass += '{} '.format(mass_grid[n])
        out_Metal_eject_mass += '\n# Metal_eject_mass\n'
        for n in range(len(Metal_eject_mass)):
            out_Metal_eject_mass += '{} '.format(Metal_eject_mass[n])
        name_Metal_eject_mass = 'yield_tables/rearranged___/setllar_Metal_eject_mass_from_{}/{}_Z={}.txt'.format(yield_table_name, yield_table_name, metallicity)
        file_Metal_eject_mass = open(name_Metal_eject_mass, 'w')
        file_Metal_eject_mass.write(out_Metal_eject_mass)
        file_Metal_eject_mass.close()

    print('yield_tables/rearranged___/setllar_...eject_mass_from_{}/{}_Z=....txt saved'.format(yield_table_name, yield_table_name))
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
    global O_over_Mg_list, Mg_over_Fe_list, Mg_over_H_list, Fe_over_H_list, O_over_H_list, Z_over_H_list, O_over_Fe_list, M_list, Z_list
    j = 0
    while j < len(M_list):
        i = 0
        while i < len(M_list[j]):
            M_list[j][i] = math.log(M_list[j][i], 10)
            (i) = (i+1)
        (j) = (j+1)
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(1, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlim(-0.5, 2.2)
    plt.ylim(0, 2)
    i = 0
    while i < len(M_list):
        plt.plot(M_list[i], O_over_Mg_list[i], label='Z={}'.format(Z_list[i]))
        (i) = (i+1)
    O_mass_eject_SNIa = 0.148  # TNH93 0.148 i99CDD1 0.09, i99CDD2 0.06, i99W7 0.14, ivo12/13 0.09-0.1, t03 0.14, t86 0.13
    Mg_mass_eject_SNIa = 0.009  # TNH93 0.009 i99CDD1 0.0077, i99CDD2 0.0042, i99W7 0.0085, ivo12/13 0.015-0.029, t03 0.013, t86 0.016
    O_num = O_mass_eject_SNIa / 15.9994
    Mg_num = Mg_mass_eject_SNIa / 24.305
    O_over_Mg_SNIa = math.log(O_num / Mg_num, 10) - 8.69 + 7.60
    plt.plot([-0.3, 0.9], [O_over_Mg_SNIa, O_over_Mg_SNIa], ls="--", lw=2, label="SNIa")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'[O/Mg]=log($N_{O}/N_{Mg}$)-[O/Mg]$_\odot$')
    plt.tight_layout()


    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(2, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlim(-0.5, 2.2)
    plt.ylim(-2, 7)
    i = 0
    while i < len(M_list):
        plt.plot(M_list[i], Mg_over_Fe_list[i], label='Z={}'.format(Z_list[i]))
        (i) = (i+1)
    Mg_mass_eject_SNIa = 0.009  # TNH93 0.148 i99CDD1 0.09, i99CDD2 0.06, i99W7 0.14, ivo12/13 0.09-0.1, t03 0.14, t86 0.13
    Fe_mass_eject_SNIa = 0.372  #0.63 # Recchi2009 halfed to 0.372  # TNH93 0.744 i99CDD1 0.56, i99CDD2 0.76, i99W7 0.63, ivo12/13 0.62-0.67, t03 0.74, t86 0.63
    Mg_num = Mg_mass_eject_SNIa / 24.305
    Fe_num = Fe_mass_eject_SNIa / 55.845
    Mg_over_Fe_SNIa = math.log(Mg_num / Fe_num, 10) - 7.60 + 7.50
    plt.plot([-0.3, 0.9], [Mg_over_Fe_SNIa, Mg_over_Fe_SNIa], ls="--", lw=2, label="SNIa")
    plt.plot([-2, 3], [0, 0], lw=0.1)
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'[Mg/Fe]=log($N_{Mg}/N_{Fe}$)-[Mg/Fe]$_\odot$')
    plt.tight_layout()


    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(3, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlim(-0.5, 2.2)
    plt.ylim(-2, 7)
    i = 0
    while i < len(M_list):
        plt.plot(M_list[i], O_over_Fe_list[i], label='Z={}'.format(Z_list[i]))
        (i) = (i+1)
    O_over_Fe_SNIa = math.log(O_num / Fe_num, 10) - 7.60 + 7.50
    plt.plot([-0.3, 0.9], [O_over_Fe_SNIa, O_over_Fe_SNIa], ls="--", lw=2, label="SNIa")
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'[O/Fe]=log($N_{O}/N_{Fe}$)-[O/Fe]$_\odot$')
    plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(4, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlim(-0.5, 2.2)
    plt.ylim(-2, 2)
    i = 0
    while i < len(M_list):
        plt.plot(M_list[i], Mg_over_H_list[i], label='Z={}'.format(Z_list[i]))
        (i) = (i + 1)
    plt.plot([-2, 3], [0, 0], lw=0.1)
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'[Mg/H]=log($N_{Mg}/N_{H}$)-log(same)$_\odot$')
    plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(5, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlim(-0.5, 2.2)
    plt.ylim(-2, 2)
    i = 0
    while i < len(M_list):
        plt.plot(M_list[i], O_over_H_list[i], label='Z={}'.format(Z_list[i]))
        (i) = (i + 1)
    plt.plot([-2, 3], [0, 0], lw=0.1)
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'[O/H]')
    plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(6, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlim(-0.5, 2.2)
    plt.ylim(-2, 2)
    i = 0
    while i < len(M_list):
        plt.plot(M_list[i], Fe_over_H_list[i], label='Z={}'.format(Z_list[i]))
        (i) = (i + 1)
    plt.plot([-2, 3], [0, 0], lw=0.1)
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'[Fe/H]')
    plt.tight_layout()

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(7, figsize=(6, 5.25))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlim(-0.5, 2.2)
    plt.ylim(-2, 2)
    i = 0
    while i < len(M_list):
        plt.plot(M_list[i], Z_over_H_list[i], label='Z={}'.format(Z_list[i]))
        (i) = (i + 1)
    plt.plot([-2, 3], [0, 0], lw=0.1)
    plt.legend(prop={'size': 10}, loc='best')
    plt.xlabel(r'log stellar mass [M$_\odot$]')
    plt.ylabel(r'[Z/H]')
    plt.tight_layout()


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
    Na_eject_mass_list = []
    Mg_eject_mass_list = []
    Al_eject_mass_list = []
    Si_eject_mass_list = []
    S_eject_mass_list = []
    Ca_eject_mass_list = []
    Ti_eject_mass_list = []
    Cr_eject_mass_list = []
    Mn_eject_mass_list = []
    Fe_eject_mass_list = []
    Ni_eject_mass_list = []
    Metal_eject_mass_list = []
    O_over_Mg_list = []
    Mg_over_H_list = []
    Fe_over_H_list = []
    O_over_H_list = []
    Z_over_H_list = []
    Mg_over_Fe_list = []
    O_over_Fe_list = []
    yield_table_name = "Limongi_R000" # being "WW95" or "portinari98" or "marigo01" or "Kobayashi06" or "Limongi_R000" or "Nomoto"
    function_read_file(yield_table_name)
    # funtion_plot_yields()
    print(" - Run time: %s -" % round((time.time() - start_time), 2))
