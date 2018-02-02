######## galIMF ##########


#python3 code, last update Sat 27 May
# This, galIMF, is the main part controling and operating the other two part below, IGIMF and OSGIMF.
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
#importing modules and libraries

import math
from scipy.interpolate import interp1d
from scipy.integrate import simps
import numpy as np
import csv # csv and izip/zip are used to create output files
try:
    from itertools import izip as zip
except ImportError: # will be python 3.x series
    pass



#--------------------------------------------------------------------------------------------------------------------------------

# The star mass resolution is the lower resolution among
# the resolution of histogram (resolution_histogram_relative)
# and the resolution of star generation (resolution_star_... in below)
resolution_histogram_relative = 0.01 # The star mass resolution of histogram is: the star mass * resolution_histogram_relative
#also re-defined in a test file, it scales automatically with the SFR

# function_galIMF takes in I/OS-GMF parameters and create output files
def function_galIMF(IorS, SFR, alpha3_model, delta_t, Fe_over_H, I_ecl, M_ecl_U, M_ecl_L, beta_model,
                         I_str, M_str_L, alpha_1, alpha1_model, M_turn, alpha_2, alpha2_model, M_turn2, M_str_U):
    if IorS == "I":
        global List_xi, List_M_str_for_xi_str
        Function_draw_IGIMF(SFR, alpha3_model, beta_model, delta_t, Fe_over_H,
                                               I_ecl, M_ecl_U, M_ecl_L, I_str, M_str_L, alpha_1, alpha1_model,
                                               M_turn, alpha_2, alpha2_model, M_turn2, M_str_U)
        # write data
        with open('GalIMF_IGIMF.txt', 'w') as f:
            writer = csv.writer(f, delimiter=' ')
            f.write("# IGIMF output file. It gives the IGIMF. The columns are:\n# mass xi\n\n")
            writer.writerows(
                zip(List_M_str_for_xi_str, List_xi))
        print("\n### IGIMF data generated in the file GalIMF_IGIMF.txt ###\n")
        return
    elif IorS =="OS":
        global mass_range_center, mass_range, mass_range_upper_limit, mass_range_lower_limit, star_number
        sample_for_one_epoch(SFR, alpha3_model, delta_t, I_ecl, M_ecl_U, M_ecl_L, beta_model,
                         I_str, M_str_L, alpha_1, alpha1_model, M_turn, alpha_2, alpha2_model, M_turn2, Fe_over_H, M_str_U)
        Function_draw(SFR, M_str_L, M_str_U, M_ecl_L, resolution_histogram_relative)
        function_make_drop_line()
        # write data
        function_draw_histogram()
        with open('GalIMF_OSGIMF.txt', 'w') as f:
            writer = csv.writer(f, delimiter=' ')
            f.write("# OSGIMF output file. It gives the star number in each mass range. The columns are:\n# mass_range_center mass_range mass_range_upper_limit mass_range_lower_limit star_number_in_the_mass_range\n\n")
            writer.writerows(
                zip(mass_range_center, mass_range, mass_range_upper_limit, mass_range_lower_limit, star_number))
        print("\n### OSGIMF data generated in the file GalIMF_OSGIMF.txt ###\n")
        return
    else:
        print("Input parameter 'IorS' wrong!")
    return













######## IGIMF #########

#python3 code, last update Sat 27 May
# IGIMF is part computing IGIMF as described in Yan et al 2017
# all physical quantities which are input in some function are described in test_gimf.py scrip or in readme file
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------------
#initialization of floating length arrays
List_M_ecl_for_xi_ecl = []
List_xi_ecl = []
List_M_str_for_xi_str = []
List_xi_str = []
List_xi = []
#--------------------------------------------------------------------------------------------------------------------------------
# Function_dar_IGIMF computes the IGIMF by combining  Function_ECMF (embedded cluster mass
# function) and Function_IMF (stellar mass function in individual embedded clusters)
# equation (1) from Yan, Jerabkova, Kroupa (2017, A&A)
# function returns values of global lists:
# List_M_ecl_for_xi_ecl - list of masses, M_ecl, of embedded clusters for ECMF
# List_xi IGIMF (xi_IGIMF = dN/dm, dN number of star in a mass bin dm) values
#   by default normalized to total mass in Msun units (= SFR*10Myr)
# List_M_str_for_xi_str list of stellar masses for stellar IMF in Msun units
# List_xi_L logarithmic IGIMF (xi_IGIMF_L = dN/d log_10 m)
# List_Log_M_str - natural logarithm
def Function_draw_IGIMF(SFR, alpha3_model, beta_model, delta_t, Fe_over_H, I_ecl, M_ecl_U, M_ecl_L,
                        I_str, M_str_L, alpha_1, alpha1_model, M_turn, alpha_2, alpha2_model, M_turn2, M_str_U):
    global List_M_ecl_for_xi_ecl, List_xi, List_M_str_for_xi_str, List_xi_L, List_Log_M_str, x_IMF, y_IMF
    Function_ECMF(SFR, beta_model, delta_t, I_ecl, M_ecl_U, M_ecl_L)
    x_IMF = []
    y_IMF = []
    alpha_1_change = Function_alpha_1_change(alpha_1, alpha1_model, Fe_over_H)
    alpha_2_change = Function_alpha_2_change(alpha_2, alpha2_model, Fe_over_H)
    alpha_3_change = Function_alpha_3_change(alpha3_model, List_M_ecl_for_xi_ecl[-1], Fe_over_H)
    function_draw_xi_str(M_str_L, List_M_ecl_for_xi_ecl[-1], I_str, M_str_L, alpha_1_change,
                                                 M_turn, alpha_2_change, M_turn2, alpha_3_change, M_str_U)
    List_xi = [0] * len(x_IMF)
    number_of_ecl = len(List_M_ecl_for_xi_ecl)-1
    Function_IMF(alpha3_model, Fe_over_H, I_str, M_str_L, alpha_1_change, M_turn, alpha_2_change, M_turn2, M_str_U, number_of_ecl, 0)
    x_IMF = []
    y_IMF = []
    function_draw_xi_str(M_str_L, List_M_ecl_for_xi_ecl[-1], I_str, M_str_L, alpha_1_change,
                                                 M_turn, alpha_2_change, M_turn2, alpha_3_change, M_str_U)
    List_M_str_for_xi_str = x_IMF
    lenth = len(List_M_str_for_xi_str)
    List_xi_L = [0] * lenth
    List_Log_M_str = [0] * lenth
    Function_xi_to_xiL(lenth-1, List_xi[0])
    return

# Function_ECMF computes IMF of star clusters (ECMF - embedded cluster mass function)
# The assumed shape of ECMF is single powerlaw with slope beta (function of SFR)
# the empyrical lower limit for star cluster mass if 50 Msun
# the hypotetical upper mass limit is 10^9 Msun, but the M_ecl^max is computed, eq (12) in Yan, Jerabkova, Kroupa (2017, A&A)
def Function_ECMF(SFR, beta_model, delta_t, I_ecl, M_ecl_U, M_ecl_L):
    global List_M_ecl_for_xi_ecl, List_xi_ecl, x_ECMF, y_ECMF
    x_ECMF = []
    y_ECMF = []
    beta_change = Function_beta_change(beta_model, SFR)
    function_draw_xi_ecl(M_ecl_L, SFR, delta_t, I_ecl, M_ecl_U, M_ecl_L, beta_change)
    List_M_ecl_for_xi_ecl = x_ECMF
    del List_M_ecl_for_xi_ecl[0]
    del List_M_ecl_for_xi_ecl[-1]
    List_xi_ecl = y_ECMF
    del List_xi_ecl[0]
    del List_xi_ecl[-1]
    return

# Function_IMF computes stellar IMF in individual embedded star clusters
def Function_IMF(alpha3_model, Fe_over_H, I_str, M_str_L, alpha_1_change, M_turn, alpha_2_change, M_turn2, M_str_U, number_of_ecl, i):
    while i < number_of_ecl:
        global List_M_str_for_xi_str, List_xi_str, List_M_ecl_for_xi_ecl, x_IMF, y_IMF
        x_IMF = []
        y_IMF = []
        M_ecl = List_M_ecl_for_xi_ecl[i]
        alpha_3_change = Function_alpha_3_change(alpha3_model, M_ecl, Fe_over_H)
        # Here only alpha_3_change is recalculated as alpha1(2)_change do not depend on M_ecl thus do not change.
        function_draw_xi_str(M_str_L, M_ecl, I_str, M_str_L, alpha_1_change, M_turn,
                                                     alpha_2_change, M_turn2, alpha_3_change, M_str_U)
        List_M_str_for_xi_str = x_IMF
        List_xi_str = y_IMF
        number_of_str = len(List_M_str_for_xi_str)
        Function_update_List_xi(i, number_of_str, 0)
        (i) = (i+1)
    return


def Function_update_List_xi(i, number_of_str, j):
    while j < number_of_str:
        global List_xi, List_xi_str, List_xi_ecl, List_M_ecl_for_xi_ecl
        List_xi[j] += List_xi_str[j] * List_xi_ecl[i] * (List_M_ecl_for_xi_ecl[i+1] - List_M_ecl_for_xi_ecl[i])
        (j) = (j+1)
    return


def Function_xi_to_xiL(i, unit):
    global List_xi_L, List_xi, List_M_str_for_xi_str, List_Log_M_str
    while i > -1:
        if List_xi[i] == 0:
            List_xi[i] = 10**(-5)
        List_xi_L[i] = math.log((List_xi[i] * math.log(10) * List_M_str_for_xi_str[i] / unit * 1800), 10)
        List_Log_M_str[i] = math.log(List_M_str_for_xi_str[i] , 10)
        (i) = (i-1)
    return













############ OSGIMF #############


#-----------------------------------------------------------------------------------------
#initialization of open-lenght arrays
#-----------------------------------------------------------------------------------------
List_M_str_all_i = []
List_n_str_all_i = []
List_mass_grid_x_axis = []
List_star_number_in_mass_grid_y_axis = []
List_star_number_in_mass_grid_y_axis2 = []
List_star_number_in_mass_grid_y_axis3 = []
List_star_number_in_mass_grid_y_axis4 = []
List_mass_grid = []
List_star_number_in_mass_grid = []
#-----------------------------------------------------------------------------------------

# This function gives the stellar masses in entire galaxy in unsorted manner
# i.e. the stars are grouped in parent clusters
def sample_for_one_epoch(SFR, alpha3_model, delta_t, I_ecl, M_ecl_U, M_ecl_L, beta_model,
                         I_str, M_str_L, alpha_1, alpha1_model, M_turn, alpha_2, alpha2_model, M_turn2, Fe_over_H, M_str_U):
    global List_M_str_all_i, List_n_str_all_i, list_M_ecl_i
    beta_change = Function_beta_change(beta_model, SFR)
    Function_sample_cluster(SFR, delta_t, I_ecl, M_ecl_U, M_ecl_L, beta_change)
    len_of_M_ecl_list = len(list_M_ecl_i)
    List_M_str_all_i = []
    List_n_str_all_i = []
    Function_sample_star_from_clusters(alpha3_model, I_str, M_str_L, alpha_1, alpha1_model, M_turn, alpha_2, alpha2_model,
                                       M_turn2, Fe_over_H, M_str_U, len_of_M_ecl_list, 0)
    return

# Masses of formed clusters
def Function_sample_cluster(SFR, delta_t, I_ecl, M_ecl_U, M_ecl_L, beta_change):
    global list_m_ecl_i, list_n_ecl_i, list_M_ecl_i, M_max_ecl
    list_m_ecl_i = []
    list_n_ecl_i = []
    list_M_ecl_i = []
    M_max_ecl = 0
    function_sample_from_ECMF(SFR, delta_t, I_ecl, M_ecl_U, M_ecl_L, beta_change)
    return
# Stellar masses in a given star cluster
def Function_sample_star_from_clusters(alpha3_model, I_str, M_str_L, alpha_1, alpha1_model, M_turn, alpha_2, alpha2_model,
                                       M_turn2, Fe_over_H, M_str_U, len_of_M_ecl_list, i):
    while i < len_of_M_ecl_list:
        global List_M_str_all_i, List_n_str_all_i, list_m_str_i, list_n_str_i, list_M_str_i
        list_m_str_i = []
        list_n_str_i = []
        list_M_str_i = []
        alpha_1_change = Function_alpha_1_change(alpha_1, alpha1_model, Fe_over_H)
        alpha_2_change = Function_alpha_2_change(alpha_2, alpha2_model, Fe_over_H)
        alpha_3_change = Function_alpha_3_change(alpha3_model, list_M_ecl_i[i], Fe_over_H)
        function_sample_from_IMF(list_M_ecl_i[i],
                                            I_str, M_str_L, alpha_1_change, M_turn, alpha_2_change, M_turn2, alpha_3_change, M_str_U)
        List_M_str_all_i += [list_M_str_i]
        List_n_str_all_i += [list_n_str_i]
        (i) = (i+1)
    return


##################################################################################
## The sampling is finished here. Below are just sorting, binning, and plotting.##
##################################################################################

# Now star mass are recorded in individual star clusters in the "List_M_str_all_i" and "List_n_str_all_i"
# we have for the whole galaxy: cluster mass, number of cluster with certain mass
#     and for     each cluster: star    mass, number of stars   with certain mass
# Sort out all star mass in a epoch into a mass grid

# Main purporpose here is the sorting of the stellar masses and preparation for
# plotting output
def Function_draw(SFR, M_str_low, M_str_up, M_ecl_low, resolution_histogram_relative):
    M_low = min(M_str_low, M_ecl_low)
    global List_mass_grid, List_star_number_in_mass_grid, List_mass_grid_x_axis, List_star_number_in_mass_grid_y_axis
    # for all stars
    List_mass_grid = []
    Function_mass_grid(SFR, M_str_up, M_low, resolution_histogram_relative)
    List_mass_grid += [M_low]
    List_star_number_in_mass_grid = [0] * (len(List_mass_grid) - 1)
    Function_sort_out_star_mass(0)
    ##########
    List_mass_grid_x_axis = [M_str_up]
    make_mass_grid_x_axis(1)
    List_mass_grid_x_axis += [M_low]
    List_star_number_in_mass_grid_y_axis = []
    make_star_number_in_mass_grid_y_axis(0)
    List_mass_grid_x_axis = [List_mass_grid_x_axis[0]] + List_mass_grid_x_axis
    List_mass_grid_x_axis += [List_mass_grid_x_axis[-1]]
    List_star_number_in_mass_grid_y_axis = [0.0000001] + List_star_number_in_mass_grid_y_axis
    List_star_number_in_mass_grid_y_axis += [0.0000001]
    # for most massive star
    global List_mass_grid2, List_star_number_in_mass_grid2, List_mass_grid_x_axis2, List_star_number_in_mass_grid_y_axis2
    List_mass_grid2 = List_mass_grid
    List_star_number_in_mass_grid2 = [0] * (len(List_mass_grid2) - 1)
    Function_sort_out_star_mass2(0)
    ##########
    List_star_number_in_mass_grid_y_axis2 = []
    make_star_number_in_mass_grid_y_axis2(0)
    List_star_number_in_mass_grid_y_axis2 = [0.0000001] + List_star_number_in_mass_grid_y_axis2
    List_star_number_in_mass_grid_y_axis2 += [0.0000001]
    ###################################
    global List_mass_grid3, List_star_number_in_mass_grid3, List_mass_grid_x_axis3, List_star_number_in_mass_grid_y_axis3
    List_mass_grid3 = List_mass_grid
    List_star_number_in_mass_grid3 = [0] * (len(List_mass_grid3) - 1)
    Function_sort_out_star_mass3(0)
    ##########
    List_star_number_in_mass_grid_y_axis3 = []
    make_star_number_in_mass_grid_y_axis3(0)
    List_star_number_in_mass_grid_y_axis3 = [0.0000001] + List_star_number_in_mass_grid_y_axis3
    List_star_number_in_mass_grid_y_axis3 += [0.0000001]
    ###################################
    global List_mass_grid4, List_star_number_in_mass_grid4, List_mass_grid_x_axis4, List_star_number_in_mass_grid_y_axis4
    List_mass_grid4 = List_mass_grid
    List_star_number_in_mass_grid4 = [0] * (len(List_mass_grid4) - 1)
    Function_sort_out_star_mass4(0)
    ##########
    List_star_number_in_mass_grid_y_axis4 = []
    make_star_number_in_mass_grid_y_axis4(0)
    List_star_number_in_mass_grid_y_axis4 = [0.0000001] + List_star_number_in_mass_grid_y_axis4
    List_star_number_in_mass_grid_y_axis4 += [0.0000001]
    return


### make a mass grid ###

def Function_mass_grid(SFR, mass, M_str_low, resolution_histogram_relative):
    while mass > M_str_low:
        global List_mass_grid
        List_mass_grid += [mass]
        (mass) = (mass * (1-resolution_histogram_relative))
        # we find it is useful to use the following form of mass grid sometimes.
        # One can apply this alternative form by quote the above line (add a # in front of the line) and unquote the below.
        #(mass) = (mass * (0.967 + math.log(SFR, 10) / 400) / (math.log(mass + 1) ** 2 / (2 ** (math.log(SFR, 10) + 6.85) - 1) + 1))
    return

# Count the number of star in each grid
def Function_sort_out_star_mass(i):
    while i < len(List_M_str_all_i):
        global l
        l = 0
        SubFunction_sort_out(i, 0)
        (i)=(i+1)
    return
def Function_sort_out_star_mass2(i):
    while i < len(List_M_str_all_i):
        global l
        l = 0
        SubFunction_sort_out2(i, 0)
        (i)=(i+1)
    return
def Function_sort_out_star_mass3(i):
    while i < len(List_M_str_all_i):
        global l
        l = 0
        SubFunction_sort_out3(i, 1)
        (i)=(i+1)
    return
def Function_sort_out_star_mass4(i):
    while i < len(List_M_str_all_i):
        global l
        l = 0
        SubFunction_sort_out4(i, 2)
        (i)=(i+1)
    return

def SubFunction_sort_out(i, j):
    while j < len(List_M_str_all_i[i]):
        global l, List_n_str_all_i
        Function_find_k(i, j, l)
        List_star_number_in_mass_grid[l] += List_n_str_all_i[i][j] * list_n_ecl_i[i]
        (j)=(j+1)
    return
def SubFunction_sort_out2(i, j):
    if j < len(List_M_str_all_i[i]):
        global l
        Function_find_k(i, j, l)
        List_star_number_in_mass_grid2[l] += list_n_ecl_i[i]
    return
def SubFunction_sort_out3(i, j):
    if j < len(List_M_str_all_i[i]):
        global l
        Function_find_k(i, j, l)
        List_star_number_in_mass_grid3[l] += list_n_ecl_i[i]
    return
def SubFunction_sort_out4(i, j):
    if j < len(List_M_str_all_i[i]):
        global l
        Function_find_k(i, j, l)
        List_star_number_in_mass_grid4[l] += list_n_ecl_i[i]
    return

def Function_find_k(i, j, k):
    while List_mass_grid[k+1] > List_M_str_all_i[i][j]:
        global l
        l = k+1
        (k) = (k+1)
    return


# Prepare for the breaking line plot
def make_mass_grid_x_axis(i):
    global List_mass_grid_x_axis, List_mass_grid
    while i < len(List_mass_grid)-1:
        List_mass_grid_x_axis += [List_mass_grid[i]]*2
        (i)=(i+1)
    return


def make_star_number_in_mass_grid_y_axis(i):
    global List_star_number_in_mass_grid_y_axis, List_star_number_in_mass_grid, List_mass_grid
    while i < len(List_star_number_in_mass_grid):
        List_star_number_in_mass_grid_y_axis += [
        List_star_number_in_mass_grid[i]/(List_mass_grid[i] - List_mass_grid[i+1])]*2
        (i)=(i+1)
    return
def make_star_number_in_mass_grid_y_axis2(i):
    global List_star_number_in_mass_grid_y_axis2, List_star_number_in_mass_grid2, List_mass_grid2
    while i < len(List_star_number_in_mass_grid2):
        List_star_number_in_mass_grid_y_axis2 += [
        List_star_number_in_mass_grid2[i]/(List_mass_grid2[i] - List_mass_grid2[i+1])]*2
        (i)=(i+1)
    return
def make_star_number_in_mass_grid_y_axis3(i):
    global List_star_number_in_mass_grid_y_axis3, List_star_number_in_mass_grid3, List_mass_grid3
    while i < len(List_star_number_in_mass_grid3):
        List_star_number_in_mass_grid_y_axis3 += [
        List_star_number_in_mass_grid3[i]/(List_mass_grid3[i] - List_mass_grid3[i+1])]*2
        (i)=(i+1)
    return
def make_star_number_in_mass_grid_y_axis4(i):
    global List_star_number_in_mass_grid_y_axis4, List_star_number_in_mass_grid4, List_mass_grid4
    while i < len(List_star_number_in_mass_grid4):
        List_star_number_in_mass_grid_y_axis4 += [
        List_star_number_in_mass_grid4[i]/(List_mass_grid4[i] - List_mass_grid4[i+1])]*2
        (i)=(i+1)
    return


def function_make_drop_line1(i):
    while i < len(List_star_number_in_mass_grid_y_axis)-1:
        if List_star_number_in_mass_grid_y_axis[i] == 0:
            List_star_number_in_mass_grid_y_axis[i] = 0.0000001
        (i) = (i+1)

def function_make_drop_line2(i):
    while i < len(List_star_number_in_mass_grid_y_axis2)-1:
        if List_star_number_in_mass_grid_y_axis2[i] == 0:
            List_star_number_in_mass_grid_y_axis2[i] = 0.0000001
        (i) = (i+1)

def function_make_drop_line3(i):
    while i < len(List_star_number_in_mass_grid_y_axis3)-1:
        if List_star_number_in_mass_grid_y_axis3[i] == 0:
            List_star_number_in_mass_grid_y_axis3[i] = 0.0000001
        (i) = (i+1)

def function_make_drop_line4(i):
    while i < len(List_star_number_in_mass_grid_y_axis4)-1:
        if List_star_number_in_mass_grid_y_axis4[i] == 0:
            List_star_number_in_mass_grid_y_axis4[i] = 0.0000001
        (i) = (i+1)

def function_make_drop_line():
    function_make_drop_line1(0)
    function_make_drop_line2(0)
    function_make_drop_line3(0)
    function_make_drop_line4(0)
    return


######################## histogram ########################

mass_range_center = []
mass_range = []
mass_range_upper_limit = []
mass_range_lower_limit = []
star_number = []

def function_draw_histogram():
    global mass_range_center, mass_range, mass_range_upper_limit, mass_range_lower_limit, star_number
    mass_range_center = []
    i = 0
    while i < len(List_mass_grid) - 1:
        mass_range_center += [
            0.5 * (List_mass_grid[i] + List_mass_grid[i + 1])]
        i = i + 1
    mass_range = []
    i = 0
    while i < len(List_mass_grid) - 1:
        mass_range += [List_mass_grid[i] - List_mass_grid[i + 1]]
        i = i + 1
    mass_range_upper_limit = []
    i = 0
    while i < len(List_mass_grid):
        mass_range_upper_limit += [List_mass_grid[i]]
        i = i + 1
    mass_range_lower_limit = []
    i = 0
    while i < len(List_mass_grid) - 1:
        mass_range_lower_limit += [List_mass_grid[i + 1]]
        i = i + 1
    star_number = List_star_number_in_mass_grid + []
    return










############## IMF #################

# This part use equations in "supplementary-document-galimf.pdf"

# The star mass resolution is the lower resolution among "relative resolution" and "absolute resolution" where
# The relative resolution = star mass * resolution_star_relative
# The absolute resolution = resolution_star_absolute
resolution_star_relative = 0.001
resolution_star_absolute = 0.001

list_m_str_i = []
list_n_str_i = []
list_M_str_i = []

def function_sample_from_IMF(M_ecl, I_str, M_L, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U):
    global list_m_str_i, list_n_str_i, list_M_str_i, M_max, M_max_function, k3, k2, k1, resolution_star_relative, resolution_star_absolute
    M_max = 0
    M_max_function = 0
    function_M_max(M_ecl, I_str, M_L, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U)
    k3 = 0
    k2 = 0
    k1 = 0
    function_k321(I_str, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U)
    list_m_str_i = []
    list_n_str_i = []
    function_m_i_str(k1, k2, k3, M_L, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_max, resolution_star_relative, resolution_star_absolute)  # equation 16
    list_M_str_i = []
    length_n = len(list_n_str_i)
    function_M_i(k1, k2, k3, M_L, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U, length_n)  # equation 18
    del list_n_str_i[0]
    return

# M_max is computed by solving simultaneously equations (3) and (4) from Yan, Jerabkova, Kroupa (2017, A&A)
def function_M_max(M_ecl, I_str, M_L, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U):
    global M_max_function, M_max, M_max_function
    M_constant = M_ecl * M_U ** (1 - alpha_3) / I_str / (1 - alpha_3) - M_turn2 ** (alpha_2 - alpha_3) * M_turn ** (
    alpha_1 - alpha_2) * (M_turn ** (2 - alpha_1) - M_L ** (2 - alpha_1)) / (2 - alpha_1) - M_turn2 ** (
    alpha_2 - alpha_3) * (M_turn2 ** (2 - alpha_2) - M_turn ** (
        2 - alpha_2)) / (2 - alpha_2) + M_turn2 ** (2 - alpha_3) / (2 - alpha_3)  # equation 14
    function_M_max_1(M_constant, M_ecl, I_str, alpha_3, M_U, M_L, 100, 10, -1)  # equation 14
    M_max_function = 1
    if M_max < M_turn2:
        M_constant2 = M_ecl * M_turn2 ** (1 - alpha_2) / I_str / (1 - alpha_2) + M_ecl * M_turn2 ** (
        alpha_3 - alpha_2) * (M_U ** (
            1 - alpha_3) - M_turn2 ** (1 - alpha_3)) / I_str / (1 - alpha_3) - M_turn ** (alpha_1 - alpha_2) * (
        M_turn ** (2 - alpha_1) - M_L ** (
            2 - alpha_1)) / (2 - alpha_1) + M_turn ** (2 - alpha_2) / (2 - alpha_2)  # equation 23
        function_M_max_2(M_constant2, M_ecl, I_str, alpha_2, M_U, M_L, 0.75, 0.1, -1)  # equation 23
        M_max_function = 2
    if M_max < M_turn:
        M_constant3 = M_ecl * M_turn ** (1 - alpha_1) / I_str / (1 - alpha_1) + M_ecl * M_turn ** (
        alpha_2 - alpha_1) * (M_turn2 ** (
            1 - alpha_2) - M_turn ** (1 - alpha_2)) / I_str / (1 - alpha_2) + M_ecl * M_turn2 ** (
        alpha_3 - alpha_2) * M_turn ** (
            alpha_2 - alpha_1) * (M_U ** (1 - alpha_3) - M_turn2 ** (1 - alpha_3)) / I_str / (1 - alpha_3) + M_L ** (
        2 - alpha_1) / (2 - alpha_1)
        # equation 27
        function_M_max_3(M_constant3, M_ecl, I_str, alpha_1, M_U, M_L, 100, 10, -1)  # equation 27
        M_max_function = 3
    if M_max < M_L:
        M_max_function = 0
        print("M_max < M_L")
    return

def function_k321(I_str, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U):
    global M_max_function, k3, k2, k1, M_max
    if M_max_function == 1:
        k3 = I_str*(1-alpha_3)/(M_U**(1-alpha_3)-M_max**(1-alpha_3))
        # equation 12
    elif M_max_function == 2:
        k3 = I_str/(M_turn2**(alpha_2-alpha_3)*(M_turn2**(1-alpha_2)-M_max**(1-alpha_2))/(1-alpha_2) + (
            M_U**(1-alpha_3)-M_turn2**(1-alpha_3))/(1-alpha_3))
        # equation 21
    elif M_max_function == 3:
        k3 = I_str/(M_turn2**(alpha_2-alpha_3) * M_turn**(alpha_1-alpha_2) * (M_turn**(1-alpha_1)-M_max**(1-alpha_1)) / (
            1-alpha_1) + M_turn2**(alpha_2-alpha_3)*(M_turn2**(1-alpha_2)-M_turn**(1-alpha_2))/(1-alpha_2) + (M_U**(
            1-alpha_3)-M_turn2**(1-alpha_3))/(1-alpha_3))
        # equation 25
    else:
        print("function_M_max went wrong")
        return
    k2 = k3*M_turn2**(alpha_2-alpha_3)  # equation 2
    k1 = k2*M_turn**(alpha_1-alpha_2)  # equation 2
    return

def function_M_max_1(M_constant, M_ecl, I_str, alpha_3, M_U, M_L, m_1, step, pm):  # equation 14
    m_1 = round(m_1, 10)  # round
    M_x = m_1**(2-alpha_3)/(2-alpha_3) + M_ecl*m_1**(1-alpha_3)/I_str/(1-alpha_3)
    if abs(M_x-M_constant) < abs(M_constant) * 10 ** (-7):
        global M_max
        M_max = m_1
    elif m_1 - step <= M_L or m_1 + step >= M_U:
        function_M_max_1(M_constant, M_ecl, I_str, alpha_3, M_U, M_L, m_1, step / 2, pm)
    elif M_x > M_constant and pm == -1:
        function_M_max_1(M_constant, M_ecl, I_str, alpha_3, M_U, M_L, m_1 - step, step, -1)
    elif M_x > M_constant and pm == 1:
        function_M_max_1(M_constant, M_ecl, I_str, alpha_3, M_U, M_L, m_1 - step / 2, step / 2, -1)
    elif M_x < M_constant and pm == 1:
        function_M_max_1(M_constant, M_ecl, I_str, alpha_3, M_U, M_L, m_1 + step, step, 1)
    elif M_x < M_constant and pm == -1:
        function_M_max_1(M_constant, M_ecl, I_str, alpha_3, M_U, M_L, m_1 + step / 2, step / 2, 1)
    return

def function_M_max_2(M_constant2, M_ecl, I_str, alpha_2, M_U, M_L, m_1, step, pm):  # equation 23
    m_1 = round(m_1, 10)  # round
    M_x = m_1 ** (2 - alpha_2) / (2 - alpha_2) + M_ecl * m_1 ** (1 - alpha_2) / I_str / (1 - alpha_2)
    if abs(M_x - M_constant2) < abs(M_constant2) * 10 ** (-7):
        global M_max
        M_max = m_1
    elif m_1 - step <= M_L or m_1 + step >= M_U:
        function_M_max_1(M_constant2, M_ecl, I_str, alpha_2, M_U, M_L, m_1, step / 2, pm)
    elif M_x > M_constant2 and pm == -1:
        function_M_max_1(M_constant2, M_ecl, I_str, alpha_2, M_U, M_L, m_1 - step, step, -1)
    elif M_x > M_constant2 and pm == 1:
        function_M_max_1(M_constant2, M_ecl, I_str, alpha_2, M_U, M_L, m_1 - step / 2, step / 2, -1)
    elif M_x < M_constant2 and pm == 1:
        function_M_max_1(M_constant2, M_ecl, I_str, alpha_2, M_U, M_L, m_1 + step, step, 1)
    elif M_x < M_constant2 and pm == -1:
        function_M_max_1(M_constant2, M_ecl, I_str, alpha_2, M_U, M_L, m_1 + step / 2, step / 2, 1)
    return

def function_M_max_3(M_constant3, M_ecl, I_str, alpha_1, M_U, M_L, m_1, step, pm):  # equation 27
    m_1 = round(m_1, 10)  # round
    M_x = m_1 ** (2 - alpha_1) / (2 - alpha_1) + M_ecl * m_1 ** (1 - alpha_1) / I_str / (1 - alpha_1)
    if abs(M_x-M_constant3) < abs(M_constant3) * 10 ** (-7):
        global M_max
        M_max = m_1
    elif m_1 - step <= M_L or m_1 + step >= M_U:
        function_M_max_1(M_constant3, M_ecl, I_str, alpha_1, M_U, M_L, m_1, step / 2, pm)
    elif M_x > M_constant3 and pm == -1:
        function_M_max_1(M_constant3, M_ecl, I_str, alpha_1, M_U, M_L, m_1 - step, step, -1)
    elif M_x > M_constant3 and pm == 1:
        function_M_max_1(M_constant3, M_ecl, I_str, alpha_1, M_U, M_L, m_1 - step / 2, step / 2, -1)
    elif M_x < M_constant3 and pm == 1:
        function_M_max_1(M_constant3, M_ecl, I_str, alpha_1, M_U, M_L, m_1 + step, step, 1)
    elif M_x < M_constant3 and pm == -1:
        function_M_max_1(M_constant3, M_ecl, I_str, alpha_1, M_U, M_L, m_1 + step / 2, step / 2, 1)
    return

def function_m_i_str(k1, k2, k3, M_L, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_max, resolution_star_relative, resolution_star_absolute):  # equation 16
    global list_m_str_i
    if M_max > 100:
        loop_m_i_first_three(k3, M_turn2, alpha_3, M_max, 0, resolution_star_relative, resolution_star_absolute, 0)
        (m_str_i, n_str_i) = cross_M_turn(k3, k2, M_turn2, alpha_3, alpha_2, list_m_str_i[-1], resolution_star_relative, resolution_star_absolute)
        loop_m_i(k2, M_turn, alpha_2, m_str_i, n_str_i, resolution_star_relative, resolution_star_absolute)
        (m_str_i, n_str_i) = cross_M_turn(k2, k1, M_turn, alpha_2, alpha_1, list_m_str_i[-1], resolution_star_relative, resolution_star_absolute)
        loop_m_i(k1, M_L, alpha_1, m_str_i, n_str_i, resolution_star_relative, resolution_star_absolute)
        cross_M_L(k1, M_L, alpha_1, list_m_str_i[-1])
        return
    elif M_max > M_turn2:
        loop_m_i(k3, M_turn2, alpha_3, M_max, 0, resolution_star_relative, resolution_star_absolute)
        (m_str_i, n_str_i) = cross_M_turn(k3, k2, M_turn2, alpha_3, alpha_2, list_m_str_i[-1], resolution_star_relative, resolution_star_absolute)
        loop_m_i(k2, M_turn, alpha_2, m_str_i, n_str_i, resolution_star_relative, resolution_star_absolute)
        (m_str_i, n_str_i) = cross_M_turn(k2, k1, M_turn, alpha_2, alpha_1, list_m_str_i[-1], resolution_star_relative, resolution_star_absolute)
        loop_m_i(k1, M_L, alpha_1, m_str_i, n_str_i, resolution_star_relative, resolution_star_absolute)
        cross_M_L(k1, M_L, alpha_1, list_m_str_i[-1])
        return
    elif M_max > M_turn:
        loop_m_i(k2, M_turn, alpha_2, M_max, 0, resolution_star_relative, resolution_star_absolute)
        (m_str_i, n_str_i) = cross_M_turn(k2, k1, M_turn, alpha_2, alpha_1, list_m_str_i[-1], resolution_star_relative, resolution_star_absolute)
        loop_m_i(k1, M_L, alpha_1, m_str_i, n_str_i, resolution_star_relative, resolution_star_absolute)
        cross_M_L(k1, M_L, alpha_1, list_m_str_i[-1])
        return
    else:
        loop_m_i(k1, M_L, alpha_1, M_max, 0, resolution_star_relative, resolution_star_absolute)
        cross_M_L(k1, M_L, alpha_1, list_m_str_i[-1])
        return

def function_get_n_new_str(m_i, k, alpha, m_i_plus_n, n_i, resolution_star_relative, resolution_star_absolute):
    while m_i - m_i_plus_n < max(resolution_star_relative * m_i, resolution_star_absolute):
        n_new = round(n_i * 1.05 + 1)
        m_i_plus_n_new = (m_i ** (1 - alpha) - n_new * (1 - alpha) / k) ** (1 / (1 - alpha))
        (m_i_plus_n, n_i) = (m_i_plus_n_new, n_new)
    return m_i_plus_n, n_i

def loop_m_i_first_three(k, M_low, alpha, m_i, n_i, resolution_star_relative, resolution_star_absolute, count):
    while m_i > M_low:
        global list_m_str_i, list_n_str_i, n_turn
        list_m_str_i += [m_i]
        list_n_str_i += [n_i]
        m_i_plus_n = (m_i ** (1 - alpha) - n_i * (1 - alpha) / k) ** (1 / (1 - alpha))
        if count < 3:
            m_i_plus_n = (m_i ** (1 - alpha) - (1 - alpha) / k) ** (1 / (1 - alpha))
            n_turn = n_i
            (m_i, n_i, count) = (m_i_plus_n, 1, (count+1))
        elif m_i - m_i_plus_n > max(resolution_star_relative * m_i, resolution_star_absolute):
            n_turn = n_i
            (m_i, n_i) = (m_i_plus_n, n_i)
        else:
            (m_i_plus_n_new, n_turn) = function_get_n_new_str(m_i, k, alpha, m_i_plus_n, n_i, resolution_star_relative, resolution_star_absolute)
            (m_i, n_i) = (m_i_plus_n_new, n_turn)

def loop_m_i(k, M_low, alpha, m_i, n_i, resolution_star_relative, resolution_star_absolute):
    while m_i > M_low:
        global list_m_str_i, list_n_str_i, n_turn
        list_m_str_i += [m_i]
        list_n_str_i += [n_i]
        a = m_i ** (1 - alpha) - n_i * (1 - alpha) / k
        if a > 0:
            b = 1 / (1 - alpha)
            m_i_plus_n = a ** b
            if m_i - m_i_plus_n > max(resolution_star_relative * m_i, resolution_star_absolute):
                (m_i, n_i) = (m_i_plus_n, n_i)
            else:
                (m_i_plus_n_new, n_turn) = function_get_n_new_str(m_i, k, alpha, m_i_plus_n, n_i, resolution_star_relative, resolution_star_absolute)
                (m_i, n_i) = (m_i_plus_n_new, n_turn)
        else:
            return


def cross_M_turn(k_before, k_after, M_cross, alpha_before, alpha_after, m_i, resolution_star_relative, resolution_star_absolute):
    global n_turn
    n_before = int(k_before/(1-alpha_before)*(m_i**(1-alpha_before)-M_cross**(1-alpha_before)))
    m_before_cross = (m_i ** (1 - alpha_before) - n_before * (1 - alpha_before) / k_before) ** (1 / (1 - alpha_before))
    a = (M_cross**(1-alpha_after)+k_before/k_after*(1-alpha_after)/(1-alpha_before)*(m_before_cross**(
        1-alpha_before)-M_cross**(1-alpha_before))-(1-alpha_after)/k_after)
    if a > 0:
        m_after_cross = a ** (1/(1-alpha_after))
        n_after = int(0.9*(n_turn - n_before - 1))
        m_after_cross_plus_n_after = (m_after_cross ** (1 - alpha_after) - n_after * (1 - alpha_after) / k_after) ** (1 / (1 - alpha_after))
        if m_i - m_after_cross_plus_n_after > max(resolution_star_relative * m_i, resolution_star_absolute):
            return (m_after_cross_plus_n_after, n_before + 1 + n_after)
        else:
            (m_after_cross_plus_n_new, n_after_new) = function_get_n_new_str_cross(
                m_i, m_after_cross, k_after, alpha_after, m_after_cross_plus_n_after, n_after, resolution_star_relative, resolution_star_absolute)
            return (m_after_cross_plus_n_new, n_before + 1 + n_after_new)
    else:
        return (0, 0)


def function_get_n_new_str_cross(m_i, m_after_cross, k, alpha, m_after_cross_plus_n, n_i, resolution_star_relative, resolution_star_absolute):
    while m_i - m_after_cross_plus_n < max(resolution_star_relative * m_i, resolution_star_absolute):
        n_after_new = round(n_i * 1.05 + 1)
        m_after_cross_plus_n_new = (m_after_cross ** (1 - alpha) - n_after_new * (1 - alpha) / k) ** (1 / (1 - alpha))
        (m_after_cross_plus_n, n_i) = (m_after_cross_plus_n_new, n_after_new)
    return m_after_cross_plus_n, n_i


def cross_M_L(k_1, M_L, alpha_1, m_i):  # equation 19
    global list_m_str_i, list_n_str_i
    n_i = int(k_1 / (1 - alpha_1) * (m_i ** (1 - alpha_1) - M_L ** (1 - alpha_1)))
    list_m_str_i += [M_L]
    list_n_str_i += [n_i]
    return


def function_M_i(k1, k2, k3, M_L, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U, length_n):  # equation 18
    global list_m_str_i, new_i, list_M_str_i, M_max, list_n_str_i
    new_i = 0
    if M_max > M_turn2:
        loop_M_i(k3, M_turn2, alpha_3, new_i)
        cross_M_turn2(k3, k2, M_turn2, alpha_3, alpha_2, new_i)
        if new_i + 1 < len(list_m_str_i):
            loop_M_i(k2, M_turn, alpha_2, new_i)
            if list_n_str_i[new_i + 1] > 0:
                cross_M_turn2(k2, k1, M_turn, alpha_2, alpha_1, new_i)
                if new_i + 1 < len(list_m_str_i):
                    loop_M_i(k1, M_L, alpha_1, new_i)
                    if list_n_str_i[new_i+1] == 0:
                        return
                    else:
                        M_i = k1 / (2 - alpha_1) * (list_m_str_i[new_i] ** (2 - alpha_1) - list_m_str_i[new_i + 1] ** (2 - alpha_1)) / \
                              list_n_str_i[new_i + 1]
                        list_M_str_i += [M_i]
                        return
    elif M_max > M_turn:
        loop_M_i(k2, M_turn, alpha_2, new_i)
        cross_M_turn2(k2, k1, M_turn, alpha_2, alpha_1, new_i)
        loop_M_i(k1, M_L, alpha_1, new_i)
        if list_n_str_i[new_i+1] == 0:
            return
        else:
            M_i = k1 / (2 - alpha_1) * (list_m_str_i[new_i] ** (2 - alpha_1) - list_m_str_i[new_i + 1] ** (
                2 - alpha_1)) / list_n_str_i[new_i + 1]
            list_M_str_i += [M_i]
            return
    else:
        loop_M_i(k1, M_L, alpha_1, new_i)
        if list_n_str_i[new_i+1] == 0:
            return
        else:
            M_i = k1 / (2 - alpha_1) * (list_m_str_i[new_i] ** (2 - alpha_1) - list_m_str_i[new_i + 1] ** (
                2 - alpha_1)) / list_n_str_i[new_i + 1]
            list_M_str_i += [M_i]
            return


def loop_M_i(k, M_low, alpha, i):
    global list_m_str_i, list_n_str_i, list_M_str_i, new_i
    while list_m_str_i[i+1] > M_low:
        M_i = k/(2-alpha)*(list_m_str_i[i]**(2-alpha)-list_m_str_i[i+1]**(2-alpha))/list_n_str_i[i+1]
        list_M_str_i += [M_i]
        new_i = i + 1
        (i)=(new_i)

def cross_M_turn2(k_before, k_after, M_cross, alpha_before, alpha_after, i):
    global list_m_str_i, list_n_str_i, list_M_str_i, new_i
    M_i = k_before / (2 - alpha_before) * (list_m_str_i[i] ** (2 - alpha_before) - M_cross ** (2 - alpha_before)
            ) / list_n_str_i[i + 1] + k_after / (2 - alpha_after) * (M_cross ** (2 - alpha_after
            ) - list_m_str_i[i + 1] ** (2 - alpha_after)) / list_n_str_i[i + 1]
    list_M_str_i += [M_i]
    new_i = i + 1
    return


################# draw IMF without sampling #################

def k_str(M_ecl, I_str, M_L, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U):
    global M_max, M_max_function, k3, k2, k1
    M_max = 0
    M_max_function = 0
    function_M_max(M_ecl, I_str, M_L, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U)
    k3 = 0
    k2 = 0
    k1 = 0
    function_k321(I_str, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U)
    return

x_IMF = []
y_IMF = []

def function_draw_xi_str(M_str, M_ecl, I_str, M_L, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U):
    global x_IMF, y_IMF, k1, k2, k3, M_max
    k_str(M_ecl, I_str, M_L, alpha_1, M_turn, alpha_2, M_turn2, alpha_3, M_U)
    function_draw_xi_str_loop(M_str, alpha_1, M_turn, alpha_2, M_turn2, alpha_3)
    # normalization corection
    mass_int = np.logspace(np.log10(x_IMF[0] + x_IMF[0] / 1000),
                           np.log10(x_IMF[-1] - x_IMF[-1] / 1000), 100)
    i_int = interp1d(x_IMF, y_IMF)
    y_IMF_int = i_int(mass_int)
    m_tot = simps(y_IMF_int * mass_int, mass_int)
    for i in range(len(y_IMF)):
        y_IMF[i]=y_IMF[i]/m_tot*M_ecl
    return

def function_draw_xi_str_loop(M_str, alpha_1, M_turn, alpha_2, M_turn2, alpha_3):
    global x_IMF, y_IMF, k1, k2, k3, M_max
    while M_str < M_max:
        x_IMF += [M_str]
        if M_str > M_turn2:
            xi = k3 * M_str ** (-alpha_3)
        elif M_str > M_turn:
            xi = k2 * M_str ** (-alpha_2)
        else:
            xi = k1 * M_str ** (-alpha_1)
        y_IMF += [xi]
        (M_str) = (1.02 * M_str)
    return














########### alpha ###########

# This part specify the IMF slope dependence on given property of the embedded star cluster (e.g, [Fe/H]).
# The IMF is devided into three different slopes with power-index alpha_1, alpha_2, and alpha_3.
# By defult, alpha_1 ranging from 0.08 to 0.5 solar mass,
#            alpha_2 ranging from 0.5 to 1 solar mass,
#            alpha_3 ranging from 1 to 150 solar mass.

def Function_alpha_1_change(alpha_1, alpha1_model, Fe_over_H):
    if (alpha1_model == 0):
        return alpha_1
    elif (alpha1_model == 1):  # This model is based on Marks et al.(2012), MNRAS, 422, 2246
        alpha_1_change = alpha_1 + 0.5 * Fe_over_H  # see 
        return alpha_1_change
    else:
        print("alpha1_model: %s, do not exist.\nCheck file 'galIMF.py'" % (alpha1_model))
        return


def Function_alpha_2_change(alpha_2, alpha2_model, Fe_over_H):
    if (alpha2_model == 0):
        return alpha_2
    elif (alpha2_model == 1):  # This model is based on Marks et al.(2012), MNRAS, 422, 2246
        alpha_2_change = alpha_2 + 0.5 * Fe_over_H
        return alpha_2_change
    else:
        print("alpha2_model: %s, do not exist.\nCheck file 'galIMF.py'" % (alpha2_model))
        return


def Function_alpha_3_change(alpha3_model, M_ecl, Fe_over_H):
    if (alpha3_model == 0):
        default_alpha3 = 2.3
        # print("alpha_3 is set to be a constant: %s, as this is the default alpha_3 value for alpha3_model 0.\nFor more options regarding alpha_3 variation, please check file 'galIMF.py'" % (default_alpha3))
        return default_alpha3
    elif (alpha3_model == 1):  # This model is based on Marks et al.(2012), MNRAS, 422, 2246
        rho = 10 ** (0.61 * math.log(M_ecl, 10) + 2.85)
        if rho < 9.5 * 10 ** 4:
            alpha_3_change = 2.3
        else:
            alpha_3_change = 1.86 - 0.43 * math.log(rho / 10 ** 6, 10)
        return alpha_3_change
    elif (alpha3_model == 2):  # This model is based on Marks et al.(2012), MNRAS, 422, 2246
        rho = 10 ** (0.61 * math.log(M_ecl, 10) + 2.85)
        x = -0.1405 * Fe_over_H + 0.99 * math.log(rho / 10 ** 6, 10)
        if x < -0.87:
            alpha_3_change = 2.3
        else:
            alpha_3_change = -0.41 * x + 1.94
        return alpha_3_change
    else:
        # print("alpha_3 is set to be a constant: %s, as this is the input value of parameter 'alpha3_model'.\nFor more options regarding alpha_3 variation, please check file 'galIMF.py'" % (alpha3_model))
        return alpha3_model
















########## ECMF #########

# This part gives the cluster masses according to file "supplementary-document-galimf.pdf".

# The code is only valid when SFR > 3 * 10^(-10) solar / year.

# Inputs:
# SFR, delta_t, I, M_U, M_L, \beta

# step 1
# use equation 13 or 17
# give first integration limit m_1 i.e. M_max_ecl

# step 2
# use equation 10 or 14
# give k

# step 3
# use equation 21
# give every integration limit m_i and the number of stars in this region n_i

# step 4
# use equation 22 or 23
# give every cluster mass M_i

# Outputs:
# list of star mass "list_M_ecl_i"
# and the number of star with each mass "list_n_ecl_i"


################### sample cluster from ECMF #####################

resolution_cluster_relative = 0.001 # The mass resolution of a embedded cluster with mass M is: M * resolution_cluster_relative.
list_m_ecl_i = []
list_n_ecl_i = []
list_M_ecl_i = []
M_max_ecl = 0

def function_sample_from_ECMF(SFR, delta_t, I_ecl, M_U, M_L, beta):
    global list_m_ecl_i, list_n_ecl_i, list_M_ecl_i, M_max_ecl, resolution_cluster_relative
    M_tot = SFR * delta_t * 10**6  # units in Myr
    if beta == 2:
        M_max_ecl = 0
        function_M_max_ecl_2(M_tot, I_ecl, M_U, M_L, 10**8, 10**7, -1)  # equation 44
        k = I_ecl/(1/M_max_ecl-1/M_U)  # equation 41
        list_m_ecl_i = [M_max_ecl]
        list_n_ecl_i = []
        function_m_i_ecl(M_max_ecl, M_L, k, beta, 1)  # equation 48
        list_M_ecl_i = []
        length_n = len(list_n_ecl_i)
        function_M_i_2(k, 0, length_n)  # equation 50
    else:
        M_max_ecl = 0
        function_M_max_ecl_not_2(M_tot, I_ecl, M_U, M_L, beta, 10**8, 10**7, -1)  # equation 40
        k = I_ecl * (1 - beta) / (M_U ** (1 - beta) - M_max_ecl ** (1 - beta))  # equation 37
        list_m_ecl_i = [M_max_ecl]
        list_n_ecl_i = []
        function_m_i_ecl(M_max_ecl, M_L, k, beta, 1)  # equation 48
        list_M_ecl_i = []
        length_n = len(list_n_ecl_i)
        function_M_i_not_2(k, beta, 0, length_n)  # equation 49
    return


def function_M_max_ecl_2(M_tot, I_ecl, M_U, M_L, m_1, step, pm):  # equation 44
    m_1 = round(m_1, 10)  # round makes the code only valid when SFR > 3 * 10^(-10) solar / year
    M_x = I_ecl * (math.log(m_1) - math.log(M_L)) / (1 / m_1 - 1 / M_U)
    if M_tot * (1. + 10 ** (-5)) > M_x > M_tot * (1- 10 ** (-5)):
        global M_max_ecl
        M_max_ecl = m_1
    elif m_1 - step < M_L or m_1 + step > M_U:
        function_M_max_ecl_2(M_tot, I_ecl, M_U, M_L, m_1, step/10, pm)
    elif M_x > M_tot and pm == -1:
        function_M_max_ecl_2(M_tot, I_ecl, M_U, M_L, m_1 - step, step, -1)
    elif M_x > M_tot and pm == 1:
        function_M_max_ecl_2(M_tot, I_ecl, M_U, M_L, m_1 - step/10, step/10, -1)
    elif M_x < M_tot and pm == 1:
        function_M_max_ecl_2(M_tot, I_ecl, M_U, M_L, m_1 + step, step, 1)
    elif M_x < M_tot and pm == -1:
        function_M_max_ecl_2(M_tot, I_ecl, M_U, M_L, m_1 + step/10, step/10, 1)


def function_M_max_ecl_not_2(M_tot, I_ecl, M_U, M_L, beta, m_1, step, pm):  # equation 40
    m_1 = round(m_1, 10)  # round makes the code only valid when SFR > 3 * 10^(-10) solar / year
    M_x = I_ecl * (1 - beta) / (2 - beta) * (m_1 ** (2 - beta) - M_L ** (2 - beta)) / (
    M_U ** (1 - beta) - m_1 ** (1 - beta))
    if M_tot * (1.+10**(-5)) > M_x > M_tot * (1-10**(-5)):
        global M_max_ecl
        M_max_ecl = m_1
    elif m_1 - step <= M_L or m_1 + step >= M_U:
        function_M_max_ecl_not_2(M_tot, I_ecl, M_U, M_L, beta, m_1, step/2, pm)
    elif M_x > M_tot and pm == -1:
        function_M_max_ecl_not_2(M_tot, I_ecl, M_U, M_L, beta, m_1 - step, step, -1)
    elif M_x > M_tot and pm == 1:
        function_M_max_ecl_not_2(M_tot, I_ecl, M_U, M_L, beta, m_1 - step/2, step/2, -1)
    elif M_x < M_tot and pm == 1:
        function_M_max_ecl_not_2(M_tot, I_ecl, M_U, M_L, beta, m_1 + step, step, 1)
    elif M_x < M_tot and pm == -1:
        function_M_max_ecl_not_2(M_tot, I_ecl, M_U, M_L, beta, m_1 + step/2, step/2, 1)

def function_m_i_ecl(m_i, M_L, k, beta, n_i):  # equation 48
    while m_i > M_L:
        global list_m_ecl_i, list_n_ecl_i, resolution_cluster_relative
        m_i_plus_n = (m_i**(1-beta) - n_i * (1-beta) / k)**(1/(1-beta))
        if m_i_plus_n < M_L:
            list_m_ecl_i += [M_L]
            n_L = int((m_i**(1-beta) - M_L**(1-beta)) * k / (1-beta))
            if n_L == 0:
                return
            else:
                list_n_ecl_i += [n_L]
                return
        elif m_i - m_i_plus_n > resolution_cluster_relative * m_i:
            list_m_ecl_i += [m_i_plus_n]
            list_n_ecl_i += [n_i]
            (m_i, n_i) = (m_i_plus_n, n_i)
        else:
            (m_i_plus_n_new, n_new) = function_get_n_new_ecl(m_i, k, beta, m_i_plus_n, n_i)
            list_m_ecl_i += [m_i_plus_n_new]
            list_n_ecl_i += [n_new]
            (m_i, n_i) = (m_i_plus_n_new, n_new)
    return


def function_get_n_new_ecl(m_i, k, beta, m_i_plus_n, n_i):
    while m_i - m_i_plus_n < resolution_cluster_relative * m_i:
        n_new = round(n_i * 1.05 + 1)
        m_i_plus_n_new = (m_i ** (1 - beta) - n_new * (1 - beta) / k) ** (1 / (1 - beta))
        (m_i_plus_n, n_i) = (m_i_plus_n_new, n_new)
    return m_i_plus_n, n_i


def function_M_i_2(k, i, length_n):  # equation 50
    while i < length_n:
        global list_m_ecl_i, list_n_ecl_i, list_M_ecl_i
        M_i = k * (math.log(list_m_ecl_i[i]) - math.log(list_m_ecl_i[i+1])) / list_n_ecl_i[i]
        list_M_ecl_i += [M_i]
        (i) = (i+1)
    return


def function_M_i_not_2(k, beta, i, length_n):  # equation 49
    while i < length_n:
        global list_m_ecl_i, list_n_ecl_i, list_M_ecl_i
        M_i = k / (2-beta) * (list_m_ecl_i[i]**(2-beta)-list_m_ecl_i[i+1]**(2-beta)) / list_n_ecl_i[i]
        list_M_ecl_i += [M_i]
        (i) = (i+1)
    return

################### draw ECMF without sampling #####################

def k_ecl(SFR, delta_t, I_ecl, M_U, M_L, beta):
    global M_max_ecl
    M_tot = SFR * delta_t * 10 ** 6  # units in Myr
    if beta == 2:
        M_max_ecl = 0
        function_M_max_ecl_2(M_tot, I_ecl, M_U, M_L, 10**8, 10**7, -1)  # equation 44
        k = I_ecl/(1/M_max_ecl-1/M_U)  # equation 41
    else:
        M_max_ecl = 0
        function_M_max_ecl_not_2(M_tot, I_ecl, M_U, M_L, beta, 10**8, 10**7, -1)  # equation 40
        k = I_ecl * (1 - beta) / (M_U ** (1 - beta) - M_max_ecl ** (1 - beta))  # equation 37
    return k

x_ECMF = []
y_ECMF = []

def function_draw_xi_ecl(M_ecl, SFR, delta_t, I_ecl, M_U, M_L, beta):
    global x_ECMF, y_ECMF
    k = k_ecl(SFR, delta_t, I_ecl, M_U, M_L, beta)
    function_draw_xi_ecl_loop(M_ecl, k, M_U, beta)
    # normalization corection
    mass_int = np.logspace(np.log10(x_ECMF[0] + x_ECMF[0] / 1000),
                           np.log10(x_ECMF[-1] - x_ECMF[-1] / 1000), 100)
    i_int = interp1d(x_ECMF, y_ECMF)
    y_ECMF_int = i_int(mass_int)
    m_tot = simps(y_ECMF_int * mass_int, mass_int)
    for i in range(len(y_ECMF)):
        y_ECMF[i] = y_ECMF[i] / m_tot * SFR * 10**7
    # add boundary
    x_ECMF = [x_ECMF[0]] + x_ECMF
    x_ECMF += [x_ECMF[-1]]
    y_ECMF = [0.000000001] + y_ECMF
    y_ECMF += [0.000000001]
    return

def function_draw_xi_ecl_loop(M_ecl, k, M_U, beta):
        while M_ecl < M_max_ecl:
            global x_ECMF, y_ECMF
            x_ECMF += [M_ecl]
            xi = k * M_ecl ** (-beta)
            y_ECMF += [xi]
            (M_ecl) = (1.01 * M_ecl)
        return









########## beta ###########



def Function_beta_change(beta_model, SFR):
    if (beta_model == 0):
        default_beta = 2.00000001
        return default_beta
    elif (beta_model == 1):
        beta_change = -0.106 * math.log(SFR, 10) + 2.00000001
        return beta_change
    elif (beta_model == 2):
        if SFR > 1:
            beta_change = -0.106 * math.log(SFR, 10) + 2.00000001
        else:
            beta_change = 2.00000001
        return beta_change
    else:
        return beta_model
