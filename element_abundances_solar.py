# This function returns the customary astronomical scale for logarithmic abundances of the sun,
# that is, log(N_X/N_H)+12

# reference:
# Asplund, Martin; Grevesse, Nicolas; Sauval, A. Jacques; Scott, Pat (2009). ARAA 47 (1): 481–522.
# Anders, E., & Grevesse, N. 1989 is applied in WW95, Geochim. Cosmochim. Acta, 53, 197


def function_solar_element_abundances(reference_name, element_name):
    if reference_name == 'Anders1989':
        if element_name == "H":
            solar_element_abundances = 12
        elif element_name == "He":
            solar_element_abundances = 10.99
        elif element_name == "C":
            solar_element_abundances = 8.56
        elif element_name == "N":
            solar_element_abundances = 8.05
        elif element_name == "O":
            solar_element_abundances = 8.93
        elif element_name == "Ne":
            solar_element_abundances = 8.09
        elif element_name == "Mg":
            solar_element_abundances = 7.58
        elif element_name == "Si":
            solar_element_abundances = 7.55
        elif element_name == "S":
            solar_element_abundances = 7.21
        elif element_name == "Ca":
            solar_element_abundances = 6.36
        elif element_name == "Fe":
            solar_element_abundances = 7.67
        else:
            print("Wrong/unknown element name for function_solar_element_abundances; Anders1989")
            solar_element_abundances = None
    elif reference_name == 'Asplund2009':
        if element_name == "H":
            solar_element_abundances = 12
        elif element_name == "He":
            solar_element_abundances = 10.93
        elif element_name == "C":
            solar_element_abundances = 8.43
        elif element_name == "N":
            solar_element_abundances = 7.83
        elif element_name == "O":
            solar_element_abundances = 8.69
        elif element_name == "Ne":
            solar_element_abundances = 7.93
        elif element_name == "Mg":
            solar_element_abundances = 7.60
        elif element_name == "Si":
            solar_element_abundances = 7.51
        elif element_name == "S":
            solar_element_abundances = 7.12
        elif element_name == "Ca":
            solar_element_abundances = 6.34
        elif element_name == "Fe":
            solar_element_abundances = 7.50
        else:
            print("Wrong/unknown element name for function_solar_element_abundances; Asplund2009")
            solar_element_abundances = None
    elif reference_name == 'Anders1989_mass':
        if element_name == "H":
            solar_element_abundances = 0.70683
        elif element_name == "He":
            solar_element_abundances = 0.27431
        elif element_name == "Metal":
            solar_element_abundances = 0.01886
        else:
            print("Wrong/unknown element name for function_solar_element_abundances; Anders1989_mass")
            solar_element_abundances = None
    elif reference_name == 'Anders1989_mass_according_to_Asplund2009':
        if element_name == "H":
            solar_element_abundances = 0.7096
        elif element_name == "He":
            solar_element_abundances = 0.2691
        elif element_name == "Metal":
            solar_element_abundances = 0.0213
        else:
            print("Wrong/unknown element name for function_solar_element_abundances; Anders1989_mass_according_to_Asplund2009")
            solar_element_abundances = None
    elif reference_name == 'Asplund2009_mass':
        if element_name == "H":
            solar_element_abundances = 0.7154
        elif element_name == "He":
            solar_element_abundances = 0.2703
        elif element_name == "Metal":
            solar_element_abundances = 0.0142
        else:
            print("Wrong/unknown element name for function_solar_element_abundances; Asplund2009_mass")
            solar_element_abundances = None
    else:
        print('Wrong input reference_name for element_abundances_solar.function_solar_element_abundances.')
        solar_element_abundances = None
    return solar_element_abundances
