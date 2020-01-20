# Python3 code

# An example file that demonstrates how to construct galaxy-wide IMF
# as well as getting each stellar mass in the galaxy applying the IGIMF theory with the galIMF model.

# Made by: Yan Zhiqiang & Tereza Jerabkova

# The outputs of this example are:

#  - the comparison plot of galaxy-wide IMF, canonical IMF, and the histogram of stellar masses;
#  - the txt file containing the galaxy-wide IMF.
#  - the txt file containing the number of stars in each mass bin;

# --------------------------------------------------------------------------------------------------------------------------------
# Import modules and libraries
# --------------------------------------------------------------------------------------------------------------------------------

import galimf  # galIMF containing IGIMF function and OSGIMF function and additional computational modules
import matplotlib.pyplot as plt  # matplotlib for plotting
import numpy as np
import math
import time
from scipy.integrate import quad

# -----------------------------------------------------------------------

print("\n    ===============================\n"
      "    === example_galaxy_wide_IMF ===\n"
      "    ===============================\n")
print("    This test code serves as an example, "
      "demonstrating how to construct and visualize IGIMF and OSGIMF with GalIMF open source code "
      "for a given galaxy-wide SFR and metallicity.\n")

# --------------------------------------------------------------------------------------------------------------------------------
# Set the final out put plot size:
# --------------------------------------------------------------------------------------------------------------------------------

fig0 = plt.figure(figsize=(3.4, 2.5))


# --------------------------------------------------------------------------------------------------------------------------------
# Define code parameters necessary for the computations:
# --------------------------------------------------------------------------------------------------------------------------------

# the most crucial ones, which you most likely might want to change

SFR = float(input("    Please input the galaxy-wide SFR in solar mass per year and ended the input with the return "
                  "key.\n    A typical input SFR is from 0.0001 to 10000 (You can input 1e-4 as 0.0001).\n    "
                  "We recommend a value smaller than 1 for the first run as high SFR calculations take more time.\n\n"
                  "    SFR [Msolar/yr] = "))
# Star Formation Rate [solar mass / yr]
M_over_H = float(input("\n    Please input the metallicity, [M/H]"
                       "\n    A typical input should be smaller than 0, i.e., metal poor."
                       "\n\n    [M/H] = "))
bindw = galimf.resolution_histogram_relative = 10 ** (max((0 - math.log(SFR, 10)), 0) ** 0.2 - 1.9)
# will change the resolution of histogram for optimal sampling automatically adjusted with SFR value.

alpha3_model = 1  # IMF high-mass-end power-index model, see file 'galimf.py'
alpha_2 = 2.3  # IMF middle-mass power-index
alpha_1 = 1.3  # IMF low-mass-end power-index
alpha2_model = 1  # see file 'galimf.py'
alpha1_model = 1  # see file 'galimf.py'
beta_model = 1
M_str_L = 0.08  # star mass lower limit [solar mass]
M_str_U = 150  # star mass upper limit [solar mass]
M_turn = 0.5  # IMF power-index breaking mass [solar mass]
M_turn2 = 1.  # IMF power-index breaking mass [solar mass]
M_ecl_U = 10**9  # embedded cluster mass upper limit [solar mass]
M_ecl_L = 5.  # embedded cluster mass lower limit [solar mass]


# ----------------------------------------------------------------

# Parameters below are internal parameters of the theory.
# Read Yan, Jerabkova, Kroupa (2017, A&A...607A.126Y) carefully before change them!

delta_t = 10.  # star formation epoch [Myr]
I_ecl = 1.  # normalization factor in the Optimal Sampling condition equation
I_str = 1.  # normalization factor in the Optimal Sampling condition equation


# --------------------------------------------------------------------------------------------------------------------------------
# Construct IGIMF:
# --------------------------------------------------------------------------------------------------------------------------------

print("\n    Calculating galaxy-wide IMF......")
start_time = time.time()
galimf.function_galimf(
    "I",  # IorS ### "I" for IGIMF; "OS" for OSGIMF
    'IGIMF',  # 'R14'
    SFR,  # Star Formation Rate [solar mass / yr]
    alpha3_model,  # IMF high-mass-end power-index model, see file 'galimf.py'
    delta_t,  # star formation epoch [Myr]
    M_over_H,
    I_ecl,  # normalization factor in the Optimal Sampling condition equation
    M_ecl_U,  # embedded cluster mass upper limit [solar mass]
    M_ecl_L,  # embedded cluster mass lower limit [solar mass]
    beta_model,  # ECMF power-index model, see file 'galimf.py'
    I_str,  # normalization factor in the Optimal Sampling condition equation
    M_str_L,  # star mass lower limit [solar mass]
    alpha_1,  # IMF low-mass-end power-index
    alpha1_model,  # see file 'galimf.py'
    M_turn,  # IMF power-index change point [solar mass]
    alpha_2,  # IMF middle-mass power-index
    alpha2_model,  # see file 'galimf.py'
    M_turn2,  # IMF power-index change point [solar mass]
    M_str_U,  # star mass upper limit [solar mass]
    printout=True  # save the generated IMF
)
print("    - Galaxy-wide IMF calculation complete -")

masses = np.array(galimf.List_M_str_for_xi_str)
igimf = np.array(galimf.List_xi)

# igimf is normalized by default to a total mass formed in 10 Myr given the SFR
# to change the normalization follow the commented part of a code
# Norm = simps(igimf*masses,masses) #- normalization to a total mass
# Norm = simps(igimf,masses) #- normalization to number of stars
# Mtot1Myr = SFR*10*1.e6 #total mass formed in 10 Myr
# igimf = np.array(igimf)*Mtot1Myr/Norm

plt.plot(np.log10(masses+1.e-50), np.log10(igimf+1.e-50), color='blue', lw=2.5, label='Galaxy-wide IMF')
ylim_min = np.min(igimf+1.e-50)
ylim_max = np.max(igimf+1.e-50)
plt.ylim(np.log10(ylim_min), np.log10(ylim_max))


# --------------------------------------------------------------------------------------------------------------------------------
# Construct OSGIMF if required by interactive input:
# --------------------------------------------------------------------------------------------------------------------------------

OSrequest = input("\n    Do you wants to calculate OSGIMF? "
                  "\n    OSGIMF gives the stellar masses generated in a 10 Myr epoch "
                  "with constant inputted SFR. This may take time for high SFR input"
                  "\n    Input 1 as yes: ")

if OSrequest == "y" or OSrequest == "Y" or OSrequest == "yes" or OSrequest == "Yes" or OSrequest == "1":
    galimf.resolution_histogram_relative = bindw / float(
        input("\n    Please input the result resolution (Input 1 for the first run): \n\n"
              "Resolution = "))
    print("\n    Calculating stellar masses of the galaxy......")
    start_time = time.time()
    galimf.function_galimf(
        "OS",  # IorS ### "I" for IGIMF; "OS" for OSGIMF
        'IGIMF',  # 'R14'
        SFR,  # Star Formation Rate [solar mass / yr]
        alpha3_model,  # IMF high-mass-end power-index model, see file 'galimf.py'
        delta_t,  # star formation epoch [Myr]
        M_over_H,
        I_ecl,  # normalization factor in the Optimal Sampling condition equation
        M_ecl_U,  # embedded cluster mass upper limit [solar mass]
        M_ecl_L,  # embedded cluster mass lower limit [solar mass]
        beta_model,  # ECMF power-index model, see file 'galimf.py'
        I_str,  # normalization factor in the Optimal Sampling condition equation
        M_str_L,  # star mass lower limit [solar mass]
        alpha_1,  # IMF low-mass-end power-index
        alpha1_model,  # see file 'galimf.py'
        M_turn,  # IMF power-index change point [solar mass]
        alpha_2,  # IMF middle-mass power-index
        alpha2_model,  # see file 'galimf.py'
        M_turn2,  # IMF power-index change point [solar mass]
        M_str_U,  # star mass upper limit [solar mass]
        printout=True  # save the generated OSGIMF
    )
    print("    - Stellar masses calculation complete - Run time: %ss -" % round((time.time() - start_time), 2))
    # One can easily import data considering number of stars in each mass bin assuming optimal sampling
    mass_range_center = galimf.mass_range_center
    mass_range = galimf.mass_range
    mass_range_upper_limit = galimf.mass_range_upper_limit
    mass_range_lower_limit = galimf.mass_range_lower_limit
    star_number = galimf.star_number
    mass_range_center, mass_range, mass_range_upper_limit, mass_range_lower_limit, star_number = zip(
        *sorted(zip(mass_range_center, mass_range, mass_range_upper_limit, mass_range_lower_limit, star_number)))
    masses = np.array(galimf.List_mass_grid_x_axis) + 1.e-50
    osgimf = np.array(galimf.List_star_number_in_mass_grid_y_axis) + 1.e-50
    plt.plot(np.log10(masses), np.log10(osgimf), color='green', lw=2.5, label='Stellar mass histogram')


# --------------------------------------------------------------------------------------------------------------------------------
# Make a grid with power-law index -2.3 to compare with the resulting IMFs.
# --------------------------------------------------------------------------------------------------------------------------------

for k in range(20):
    sal_IMF = masses ** (-2.3)
    plt.plot(np.log10(masses), np.log10((1.e5*np.max(igimf)/np.max(sal_IMF))*sal_IMF)-k, c='grey', lw=0.5)

N = 100
can_imf = np.zeros(N)
masses = np.logspace(np.log10(0.08), np.log10(120), N, base=10)

for i, m in enumerate(masses):
    if m <= 0.5:
        can_imf[i] = m ** (-1.3)
    else:
        can_imf[i] = 0.5*m ** (-2.3)


def imf(mass, k_star, alpha):
    return k_star*mass*mass**(-alpha)


Norm = quad(imf, 0.08, 0.5, args=(1, 1.3))[0] + quad(imf, 0.5, 120, args=(0.5, 2.3))[0]
Mtot1Myr = SFR*10*1.e6  # total mass formed in 10 Myr
can_imf = np.array(can_imf)*Mtot1Myr/Norm
plt.plot(np.log10(masses), np.log10(can_imf), color='r', lw=2, label='canonical IMF')

if ylim_max < np.max(can_imf):
    ylim_max = np.max(can_imf)

# --------------------------------------------------------------------------------------------------------------------------------
# Plot settings
# --------------------------------------------------------------------------------------------------------------------------------

plt.xlabel('$\log{(m\,[M_{\odot}])}$')
plt.ylabel('$\log{(\\xi_{\mathrm{gal}}\,[M_{\odot}^{-1}])}$')

plt.ylim(np.log10(ylim_min), np.log10(ylim_max))
plt.xlim(math.log(0.06, 10), math.log(160, 10))

plt.legend(loc='best', ncol=1, fancybox=True, prop={'size': 7})
plt.tight_layout()
fig0.savefig('galaxy_wide_IMF_plot.pdf', dpi=250)

print("\n    Please check the prompted window for the result."
      "\n    IMFs in the plot are normalized by the same galaxy stellar mass."
      "\n    The generated plot is saved in the file galaxy_wide_IMF_plot.pdf."
      "\n    The program will finish when you close it.")

plt.show()

print("\n    Example complete.\n"
      "    ======================\n")
