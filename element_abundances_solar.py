# This function returns the customary astronomical scale for logarithmic abundances of the sun,
# that is, log(N_X/N_H)+12

# reference:
# Asplund, Martin; Grevesse, Nicolas; Sauval, A. Jacques; Scott, Pat (2009). ARAA 47 (1): 481–522.
# Anders, E., & Grevesse, N. 1989 is applied in WW95, Geochim. Cosmochim. Acta, 53, 197


def function_solar_element_abundances(element_name):
    if element_name == "H":
        solar_element_abundances = 12
    elif element_name == "He":
        solar_element_abundances = 10.99  # Anders 1989: 10.99, Asplund 2009: 10.93
    elif element_name == "C":
        solar_element_abundances = 8.56  # Anders 1989: 8.556, Asplund 2009: 8.43
    elif element_name == "N":
        solar_element_abundances = 8.05  # Anders 1989: 8.0536, Asplund 2009: 7.83
    elif element_name == "O":
        solar_element_abundances = 8.93  # Anders 1989: 8.932, Asplund 2009: 8.69
    elif element_name == "Ne":
        solar_element_abundances = 8.09  # Anders 1989: , Asplund 2009: 7.93
    elif element_name == "Mg":
        solar_element_abundances = 7.58  # Anders 1989: 7.4807, Asplund 2009: 7.60
    elif element_name == "Si":
        solar_element_abundances = 7.55  # Anders 1989: , Asplund 2009: 7.51
    elif element_name == "S":
        solar_element_abundances = 7.21  # Anders 1989: , Asplund 2009: 7.12
    elif element_name == "Ca":
        solar_element_abundances = 6.36  # Anders 1989: 6.329, Asplund 2009: 6.34
    elif element_name == "Fe":
        solar_element_abundances = 7.67  # Anders 1989: 7.4758, Asplund 2009: 7.50
    else:
        print("Wrong element name for function_solar_element_abundances")
        solar_element_abundances = None
    return solar_element_abundances
