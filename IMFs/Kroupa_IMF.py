def custom_imf(mass, time=0):  # there is no time dependence for Kroupa IMF
    if mass < 0.08:
        return 0
    elif mass < 0.5:
        return 2*mass**(-1.3)
    elif mass < 150:
        return mass**(-2.3)
    else:
        return 0