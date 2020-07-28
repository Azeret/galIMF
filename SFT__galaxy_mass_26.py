import galevo
import math
import element_abundances_solar
import multiprocessing as mp
from time import time


def simulate(imf, Log_SFR, SFEN, STF):
    Z_0 = 0.0000000142
    solar_mass_component = "Asplund2009_mass"
    Z_solar = element_abundances_solar.function_solar_element_abundances(solar_mass_component, 'Metal')

    galevo.galaxy_evol(
        imf=imf,
        STF=STF,  # unrealistic results if more star are forming at a time step than the instantaneous gas mass
        SFEN=SFEN,
        Z_0=Z_0,
        solar_mass_component=solar_mass_component,
        str_yield_table='Kobayashi06',
        IMF_name='Kroupa',
        steller_mass_upper_bound=150,
        time_resolution_in_Myr=1,
        mass_boundary_observe_low=1.5,
        mass_boundary_observe_up=8,
        SFH_model='provided',
        SFE=0.013,  # This parameter is not applied when SFH_model='provided'.
        SNIa_ON=True,
        SNIa_yield_table='Iwamoto1999',
        solar_abu_table='Asplund2009',
        high_time_resolution=None,
        plot_show=None,
        plot_save=None,
        outflow=None,
        check_igimf=None)
    end_time = time()

    log_Z_0 = round(math.log(Z_0 / Z_solar, 10), 2)
    file = open(
        'simulation_results_from_galaxy_evol/imf{}STF{}log_SFR{}SFEN{}Z_0{}/chemical_and_SN_evolution.txt'.format(imf, STF, Log_SFR,
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
    gas_Fe_over_H = [float(x) for x in data[19].split()]
    Mass_weighted_stellar_Fe_over_H = [float(x) for x in data[21].split()]
    # luminosit_weighted_stellar_Fe_over_H = [float(x) for x in data[??].split()]

    if imf == 'igimf':
        file_name = 'Metal_mass_relation_IGIMFZ'
    elif  imf == 'Kroupa':
        file_name = 'Metal_mass_relation_KroupaIMF'

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
    new_line = old_lines + "{} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(imf__, Log_SFR, SFEN, STF,
                            Alive_stellar_mass[0], dynamical_mass[0],
                            Mass_weighted_stellar_Mg_over_Fe[-1], Mass_weighted_stellar_Z_over_X[-1],
                            gas_Mg_over_Fe[-1], gas_Z_over_X[-1],
                            luminosity_weighted_stellar_Mg_over_Fe[-1], luminosit_weighted_stellar_Z_over_X[-1],
                            gas_Fe_over_H[-1], Mass_weighted_stellar_Fe_over_H[-1])
    file.write(new_line)
    file.close()
    return


# Parallelizing using Pool.map()
def a_pipeline(parameter):
    STF = parameter
    print("\n Start simulation for: SFEN={} STF={} Log_SFR={} imf={}".format(SFEN, STF, Log_SFR, imf))
    simulate(imf, Log_SFR, SFEN, STF)
    return

def a_pipeline_pair(parameters):
    imf = parameters[0]
    STF = parameters[1]
    print("\n Start simulation for: SFEN={} STF={} Log_SFR={} imf={}".format(SFEN, STF, Log_SFR, imf))
    simulate(imf, Log_SFR, SFEN, STF)
    return


if __name__ == '__main__':
    start = time()

    # Parallelizing only work for the same SFEN since SFH.txt file is the same!
    SFH_shape = 'flat'
    location = 0
    skewness = 10
    sfr_tail = 0
    imf = 'Kroupa'
    # SFEN_list = [100]
    # for SFEN in SFEN_list:
    #     Log_SFR_list = [5.0]
    #     for Log_SFR in Log_SFR_list:
    #         galevo.generate_SFH(SFH_shape, Log_SFR, SFEN, sfr_tail, skewness, location)
    #         STF_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    #         pool = mp.Pool(mp.cpu_count())
    #         pool.map(a_pipeline, [STF for STF in STF_list])
    #         pool.close()

    SFEN_list = [2, 5, 10, 20, 50, 100, 150, 200, 250, 300, 350, 400]
    for SFEN in SFEN_list:
        Log_SFR_list = [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0, 3.5, 4.0]
        for Log_SFR in Log_SFR_list:
            galevo.generate_SFH(SFH_shape, Log_SFR, SFEN, sfr_tail, skewness, location)
            STF_list = [0.1, 0.35, 0.6, 0.85, 1.1]
            pool = mp.Pool(mp.cpu_count())
            pool.map(a_pipeline, [STF for STF in STF_list])
            pool.close()

    # SFEN_list = [400]
    # for SFEN in SFEN_list:
    #     Log_SFR_list = [-2.0, 4.0]
    #     for Log_SFR in Log_SFR_list:
    #         galevo.generate_SFH(SFH_shape, Log_SFR, SFEN, sfr_tail, skewness, location)
    #         STF_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    #         pool = mp.Pool(mp.cpu_count())
    #         pool.map(a_pipeline, [STF for STF in STF_list])
    #         pool.close()
    end = time()
    print("Run time:", end - start)
