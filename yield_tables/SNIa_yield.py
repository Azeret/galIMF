# This function returns the element mass ejected for a type Ia supernova event

def function_mass_ejected(yield_reference_name, element_name):
    mass_ejected = 0
    if yield_reference_name == 'Thielemann1993':
        # Reference: Thielemann et al. (1993)
        # Values adopted from
        # Gibson, B. K., Loewenstein, M., & Mushotzky, R. F. 1997, MNRAS, 290, 623, their TNH93 dataset
        if element_name == "O":
            mass_ejected = 0.148
        elif element_name == "Ne":
            mass_ejected = 0.005
        elif element_name == "Mg":
            mass_ejected = 0.009
        elif element_name == "Si":
            mass_ejected = 0.158
        elif element_name == "S":
            mass_ejected = 0.086
        elif element_name == "Fe":
            mass_ejected = 0.744
        else:
            print("element {} not included in SNIa yield table {}.".format(element_name, yield_reference_name))
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
        elif element_name == "Na":
            mass_ejected = 6.8288e-5  # +-0.004
        elif element_name == "Mg":
            mass_ejected = 0.01928  # +-0.01
        elif element_name == "Al":
            mass_ejected = 0.000785
        elif element_name == "Si":
            mass_ejected = 0.248  # +-0.092
        elif element_name == "S":
            mass_ejected = 0.0935  # +-0.032
        elif element_name == "Ar":
            mass_ejected = 0.0148  # +-0.005
        elif element_name == "Ca":
            mass_ejected = 0.012  # +-0.004
        elif element_name == "Ti":
            mass_ejected = 2.5535e-4
        elif element_name == "Cr":
            mass_ejected = 0.0072  # +-0.0024
        elif element_name == "Mn":
            mass_ejected = 0.0106  # +-0.0025
        elif element_name == "Fe":
            mass_ejected = 0.68935  # +-0.21
        elif element_name == "Ni":
            mass_ejected = 0.065  # +-0.010
        else:
            print("element {} not included in SNIa yield table {}.".format(element_name, yield_reference_name))
            mass_ejected = 0
    elif yield_reference_name == 'Iwamoto1999':
        # Reference: https://ui.adsabs.harvard.edu/abs/1999ApJS..125..439I/abstract
        # Below adopt the mean value of all models (W, WDD, CDD) in their table 3
        if element_name == "C":
            mass_ejected = 0.0508  #
        elif element_name == "O":
            mass_ejected = 0.091  # (14.3+13.3+8.82+6.58+5.58+9.34+5.83)/7  # 0.133
        elif element_name == "Ne":  # Ne20
            mass_ejected = 0.00229  #
        elif element_name == "Mg":
            mass_ejected = 0.00727  # (8.5+15.8+7.55+4.47+2.62+7.72+4.2)/7  # 0.0158
        elif element_name == "Al":
            mass_ejected = 3.7214e-4  # (9.86+1.13+4.38+2.47+1.41+4.45+2.35)/7 * 1e-4
        elif element_name == "Si":
            mass_ejected = 0.201  # (1.54+1.42+2.72+2.06+1.58+2.77+1.98)/7 # 0.142
        elif element_name == "S":
            mass_ejected = 0.0914  #
        elif element_name == "Ar":
            mass_ejected = 0.0191  #
        elif element_name == "Ca":
            mass_ejected = 0.0228  # (1.19+1.81+3.1+2.43+1.88+3.18+2.38)/7 # 0.0181
        elif element_name == "Ti":
            mass_ejected = 5.3057e-4  # (2.05+3.13+7.10+6.11+5.23+7.32+6.20)e-4/7
        elif element_name == "Cr": # Cr52
            mass_ejected = 0.00773
        elif element_name == "Mn":
            mass_ejected = 0.00666  #
        elif element_name == "Fe":
            mass_ejected = 0.6747  # (6.26+6.8+5.87+7.13+7.95+5.65+7.57)/7  # 0.68
        elif element_name == "Ni":  # Ni 58
            mass_ejected = 0.0834  #
        else:
            print("element {} not included in SNIa yield table {}.".format(element_name, yield_reference_name))
            mass_ejected = 0
    elif yield_reference_name == 'Iwamoto1999_W70':
        # Reference: https://ui.adsabs.harvard.edu/abs/1999ApJS..125..439I/abstract
        # Below adopt the main isotope of W70 model
        if element_name == "C":
            mass_ejected = 0.0508
        elif element_name == "O":
            mass_ejected = 0.133
        elif element_name == "Ne":  # Ne20
            mass_ejected = 0.00229
        elif element_name == "Mg":
            mass_ejected = 0.0158
        elif element_name == "Al":
            mass_ejected = 1.31e-4
        elif element_name == "Si":
            mass_ejected = 0.142
        elif element_name == "S":
            mass_ejected = 0.0914
        elif element_name == "Ar":
            mass_ejected = 0.0191
        elif element_name == "Ca":
            mass_ejected = 0.0181
        elif element_name == "Ti": # Ti48
            mass_ejected = 3.13e-4
        elif element_name == "Cr": # Cr52
            mass_ejected = 0.00773
        elif element_name == "Mn":
            mass_ejected = 0.00666
        elif element_name == "Fe":
            mass_ejected = 0.68
        elif element_name == "Ni":  # Ni 58
            mass_ejected = 0.0834  #
        else:
            print("element {} not included in SNIa yield table {}.".format(element_name, yield_reference_name))
            mass_ejected = 0
    elif yield_reference_name == 'Iwamoto1999_W7':
        # Reference: https://ui.adsabs.harvard.edu/abs/1999ApJS..125..439I/abstract
        # Below adopt the main isotope of W70 model
        if element_name == "C":
            mass_ejected = 0.0483
        elif element_name == "O":
            mass_ejected = 0.143
        elif element_name == "Ne":  # Ne20
            mass_ejected = 0.00202  #
        elif element_name == "Mg":
            mass_ejected = 0.0085 * 5  ### Francois 2004 suggest W7 model * 5
        elif element_name == "Al":
            mass_ejected = 9.86e-4
        elif element_name == "Si":
            mass_ejected = 0.154
        elif element_name == "S":
            mass_ejected = 0.0846
        elif element_name == "Ar":
            mass_ejected = 0.0147
        elif element_name == "Ca":
            mass_ejected = 0.0119
        elif element_name == "Ti": # Ti48
            mass_ejected = 2.05e-4
        elif element_name == "Cr": # Cr52
            mass_ejected = 0.00636
        elif element_name == "Mn":
            mass_ejected = 0.00887
        elif element_name == "Fe":
            mass_ejected = 0.626
        elif element_name == "Ni":  # Ni 58
            mass_ejected = 0.11
        else:
            print("element {} not included in SNIa yield table {}.".format(element_name, yield_reference_name))
            mass_ejected = 0
    elif yield_reference_name == 'Iwamoto1999_WDD3':
        # Reference: https://ui.adsabs.harvard.edu/abs/1999ApJS..125..439I/abstract
        if element_name == "C":
            mass_ejected = 1.66e-2
        elif element_name == "O":
            mass_ejected = 5.58e-2
        elif element_name == "Ne":  # Ne20
            mass_ejected = 4.55e-4
        elif element_name == "Na":  # Na23
            mass_ejected = 3.01e-05
        elif element_name == "Mg":
            mass_ejected = 2.62e-3
        elif element_name == "Al":
            mass_ejected = 1.41e-4
        elif element_name == "Si":
            mass_ejected = 1.58e-1
        elif element_name == "S":
            mass_ejected = 9.37e-2
        elif element_name == "Ar":
            mass_ejected = 1.87e-2
        elif element_name == "Ca":
            mass_ejected = 1.88e-2
        elif element_name == "Ti": # Ti48
            mass_ejected = 5.23e-4
        elif element_name == "Cr":  # Cr52
            mass_ejected = 1.13e-2
        elif element_name == "Mn":
            mass_ejected = 6.16e-3
        elif element_name == "Fe":
            mass_ejected = 7.95e-1
        elif element_name == "Ni":  # Ni 58
            mass_ejected = 4.97e-2
        else:
            print("element {} not included in SNIa yield table {}.".format(element_name, yield_reference_name))
            mass_ejected = 0
    else:
        print('Input yield reference name for SNIa, "{}", not found.'.format(yield_reference_name))
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
