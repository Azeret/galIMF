# This function returns the element mass ejected for a type Ia supernova event


def function_mass_ejected(yield_reference_name, element_name):
    mass_ejected = 0
    if yield_reference_name == ' Thielemann1993':
        # Reference: Thielemann et al. (1993)
        # Values adopted from
        # Gibson, B. K., Loewenstein, M., & Mushotzky, R. F. 1997, MNRAS, 290, 623, their TNH93 dataset
        if element_name == "O":
            mass_ejected = 0.148  #
        elif element_name == "Ne":
            mass_ejected = 0.005  #
        elif element_name == "Mg":
            mass_ejected = 0.009  #
        elif element_name == "Si":
            mass_ejected = 0.158  #
        elif element_name == "S":
            mass_ejected = 0.086  #
        elif element_name == "Fe":
            mass_ejected = 0.744  #
        else:
            mass_ejected = 0
    elif yield_reference_name == 'Seitenzahl2013':
        # Reference: Seitenzahl et al. 2013, MNRAS, 429, 1156
        # Below adopt the mean value of all the model results in their table 2
        if element_name == "C":
            mass_ejected = 0.0073  # +-0.0047
        elif element_name == "O":
            mass_ejected = 0.11  # +-0.06
        elif element_name == "Ne":
            mass_ejected = 0.0057  # +-0.004
        elif element_name == "Mg":
            mass_ejected = 0.019  # +-0.01
        elif element_name == "Si":
            mass_ejected = 0.248  # +-0.092
        elif element_name == "S":
            mass_ejected = 0.0935  # +-0.032
        elif element_name == "Ar":
            mass_ejected = 0.0148  # +-0.005
        elif element_name == "Ca":
            mass_ejected = 0.012  # +-0.004
        elif element_name == "Cr":
            mass_ejected = 0.0072  # +-0.0024
        elif element_name == "Mn":
            mass_ejected = 0.0106  # +-0.0025
        elif element_name == "Fe":
            mass_ejected = 0.69  # +-0.21
        elif element_name == "Ni":
            mass_ejected = 0.065  # +-0.010
        else:
            mass_ejected = 0
    return mass_ejected


# # Other yield tables:
# # t86: Thielemann et al. 1986; ivo13: Seitenzahl et al. 201
# Fe_mass_eject = 0.744  # Nomoto 1984 0.613,  TNH93 0.744, i99CDD1/CDD2/W7 0.56   /0.76   /0.63,   ivo12/13 0.62-0.67,   t03 0.74,  t86 0.63
# Si_mass_eject = 0.158
# O_mass_eject = 0.148  # Nomoto 1984 0.140,  TNH93 0.148, i99CDD1/CDD2/W7 0.09   /0.06,  /0.14,   ivo12/13 0.09-0.1,    t03 0.14,  t86 0.13
# S_mass_eject = 0.086
# Mg_mass_eject = 0.009  # Nomoto 1984 0.023,  TNH93 0.009, i99CDD1/CDD2/W7 0.0077 /0.0042 /0.0085, ivo12/13 0.015-0.029, t03 0.013, t86 0.016
# Ne_mass_eject = 0.005
# #     O/Mg_mass =     # Nomoto 1984 6.0869, TNH93 16.44, i99CDD1/CDD2/W7 11.688 /14.28  /16.47,  ivo12/13 6-3.448,     t03 10.77, t86 8.125
