import galevo

print("\n    ================================\n"
      "    === example_galaxy_evolution ===\n"
      "    ================================\n")
print("    This test code serves as an example, "
      "explaining (see comments in the code) the input parameters of the galaxy chemical evolution model.\n")


Log_SFR = float(input(
    "    Please input the logarithmic star formation rate in the unit of solar mass per yr "
    "and ended the input with the return key.\n"
    "    A typical input SFR is from -4 to 4.\n\n"
    "    log_{10}(SFR [M_sun/yr]) = "))

SFH_shape = input(
    "\n\n    Please input the shape of the SFH "
    "and ended the input with the return key.\n"
    "    The input can only be: 'flat' or 'skewnorm', where the latter cost more calculation time.\n\n"
    "    ")
# Other SFH shape parameters
location = 0
skewness = 10
sfr_tail = 0

SFEN = round(float(input(
    "\n\n    Please input the characteristic star formation timescale in the unit of 10 Myr (integer only) "
    "and ended the input with the return key.\n"
    "    We recommend a value smaller than 10 for 'flat' SFH and smaller than 3 for 'skewnorm' SFH for the first run, "
    "as longer timescale calculations take more time.\n\n"
    "    SFT [10Myr] = ")))
if SFEN < 1:
    print("\n\n### Warning: Wrong input 'SFEN' smaller than 1! Correct SFEN to 1. ###\n\n")
    SFEN = 1

galevo.generate_SFH(SFH_shape, Log_SFR, SFEN, sfr_tail, skewness, location)

galevo.galaxy_evol(
    imf='igimf',
    STR=1,
    SFEN=SFEN,
    Z_0=0.00000001886,
    Z_solar=0.01886,
    str_evo_table='portinari98',
    IMF_name='Kroupa',
    steller_mass_upper_bound=150,
    time_resolution_in_Myr=1,
    mass_boundary_observe_low=1.5,
    mass_boundary_observe_up=8,
    SFH_model='provided',
    SFE=0.013,
    SNIa_ON=True,
    high_time_resolution=None,
    plot_show=True,
    plot_save=None,
    outflow=None,
    check_igimf=True)