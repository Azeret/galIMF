def custom_imf(mass, time=0):  # there is no time dependence for Salpeter IMF
    if mass < 0.1:
        return 0
    elif mass < 100:
        return mass ** (-2.35)
    else:
        return 0
