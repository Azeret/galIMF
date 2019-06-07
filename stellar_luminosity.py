# The function here gives the stellar bolometric luminosity relative to the sun [L_sun],
# assuming a simplified form using only the main-sequence luminosity as a function of mass [M_sun].
# See Yan et al. 2019 for details.
# The stellar luminosity should also be a function of Y_for_helium and Z_for_metal, which shall be added later.


def stellar_luminosity_function(mass):
    if mass < 0.23 ** (1 / 1.7):
        lum = 0.23 * mass ** 2.3
    elif mass < 1.96:
        lum = mass ** 4
    elif mass < (32000 / 1.4) ** (1 / 2.5):
        lum = 1.4 * mass ** 3.5
    else:
        lum = 32000 * mass
    return lum
