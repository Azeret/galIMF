def custom_imf(mass, time):
    change_time = 10*10**7
    change_limit = 1
    alpha_change = (change_time - time)/change_time
    if alpha_change < 0 - change_limit:
        alpha_change = 0 - change_limit
    if alpha_change > change_limit:
        alpha_change = change_limit

    if mass < 0.08:
        return 0
    elif mass < 150:
        xi = mass ** (-2.35 + alpha_change)
        return xi
    else:
        return 0
