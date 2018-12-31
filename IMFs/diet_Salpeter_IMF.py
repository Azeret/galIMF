def custom_imf(mass, time):  # there is no time dependence for Salpeter IMF
    # Bell & de Jong (2001). Salpeter IMF x = 1.35 with a flat x = 0 slope below 0.35
    # integrate this function's output xi result in the number of stars in mass limits.
    if mass < 0.35:
        xi = mass ** (-1)
    elif mass < 150:
        xi = mass ** (-2.35) * 0.35**(-1) / 0.35**(-2.35)
    else:
        xi = 0
    return xi