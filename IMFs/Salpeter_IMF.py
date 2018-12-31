def custom_imf(mass, time):  # there is no time dependence for Salpeter IMF
    if mass < 0.08:
        return 0
    elif mass < 150:
        return mass ** (-2.35)
    else:
        return 0