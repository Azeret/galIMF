def custom_imf(mass, time):  # there is no time dependence for Kroupa IMF
    if mass < 0.08:
        return 0
    elif mass < 0.5:
        return 2*mass**(-1.3)
    elif mass < 150:
        return mass**(-2.3)
    else:
        return 0


# def mass_function(mass):
#     m = custom_imf(mass)*mass
#     return m
#
# from scipy.integrate import quad
#
# m1 = quad(mass_function, 0.08, 10, limit=50)[0]
# m2 = quad(mass_function, 10, 150, limit=50)[0]
# print(m1)
# print(m2)
# print(m1/m2)
