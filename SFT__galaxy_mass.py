import galevo
import math
import element_abundances_solar


def simulate(imf, Log_SFR, SFEN, STF):
    # generate SFH:
    SFH_shape = 'flat'
    location = 0
    skewness = 10
    sfr_tail = 0

    galevo.generate_SFH(SFH_shape, Log_SFR, SFEN, sfr_tail, skewness, location)

    Z_0 = 0.00000001886
    Z_solar_table = "Anders1989_mass"
    Z_solar = element_abundances_solar.function_solar_element_abundances(Z_solar_table, 'Metal')

    galevo.galaxy_evol(
        imf=imf,
        STF=STF,  # unrealistic results if more star are forming at a time step than the instantaneous gas mass
        SFEN=SFEN,
        Z_0=Z_0,
        Z_solar_table=Z_solar_table,
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
        plot_show=None,
        plot_save=None,
        outflow=None,
        check_igimf=True)

    log_Z_0 = round(math.log(Z_0 / Z_solar, 10), 2)
    file = open(
        'simulation_results_from_galaxy_evol/imf:{}-STF:{}-log_SFR:{}-SFEN:{}-Z_0:{}.txt'.format(imf, STF, Log_SFR,
                                                                                                 SFEN, log_Z_0), 'r')
    data = file.readlines()
    file.close()

    Alive_stellar_mass = [float(x) for x in data[7].split()]
    dynamical_mass = [float(x) for x in data[11].split()]
    gas_Mg_over_Fe = [float(x) for x in data[23].split()]
    Mass_weighted_stellar_Mg_over_Fe = [float(x) for x in data[25].split()]
    luminosity_weighted_stellar_Mg_over_Fe = [float(x) for x in data[63].split()]
    gas_Z_over_X = [float(x) for x in data[39].split()]
    Mass_weighted_stellar_Z_over_X = [float(x) for x in data[41].split()]
    luminosit_weighted_stellar_Z_over_X = [float(x) for x in data[61].split()]

    file_name = 'Metal_mass_relation'
    if imf == 'igimf':
        file_name = 'Metal_mass_relation_igimf'

    file = open('simulation_results_from_galaxy_evol/{}.txt'.format(file_name), 'r')
    old_lines = file.read()
    file.close()

    file = open('simulation_results_from_galaxy_evol/{}.txt'.format(file_name), 'w')
    if imf == 'Kroupa':
        imf__ = 0
    elif imf == 'igimf':
        imf__ = 1
    else:
        imf__ = imf
    new_line = old_lines + "{} {} {} {} {} {} {} {} {} {} {} {}\n".format(imf__, Log_SFR, SFEN, STF,
                            Alive_stellar_mass[0], dynamical_mass[0],
                            Mass_weighted_stellar_Mg_over_Fe[-1], Mass_weighted_stellar_Z_over_X[-1],
                            gas_Mg_over_Fe[-1], gas_Z_over_X[-1],
                            luminosity_weighted_stellar_Mg_over_Fe[-1], luminosit_weighted_stellar_Z_over_X[-1])
    file.write(new_line)
    file.close()
    return


if __name__ == '__main__':
    imf = 'igimf'
    STF = 0.3
    Log_SFR = 0.0
    SFEN = 50

    for SFEN in [20, 400]:
        for Log_SFR in [1.0, 3.0]:
            for STR in [0.2, 0.3]:
                print("\n", imf, Log_SFR, SFEN, STF)
                simulate(imf, Log_SFR, SFEN, STF)
