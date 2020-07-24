from scipy.integrate import quad


def custom_imf_unnormalized(mass):  # there is no time dependence for Salpeter IMF
    if mass < 0.1:
        return 0
    elif mass < 100:
        return mass ** (-2.35)
    else:
        return 0


def mass_function(mass):
    return custom_imf_unnormalized(mass) * mass


integrated_mass = quad(mass_function, 0.08, 150, limit=50)[0]


def custom_imf(mass, time=0):  # normalized to a population with mass = 1 Msun
    if mass < 0.1:
        return 0
    elif mass < 100:
        return mass ** (-2.35)/integrated_mass
    else:
        return 0
