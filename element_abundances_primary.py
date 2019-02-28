import element_weight_table, element_abundances_solar

H_weight = element_weight_table.function_element_weight("H")

primary_H_mass_fraction = 0.7381
primary_He_mass_fraction = 0.2485
# Reference: 


def function_element_mass_primary_fraction(element_name, Z_0, Z_solar):
    if element_name == "H":
        element_mass_fraction = primary_H_mass_fraction
    elif element_name == "He":
        element_mass_fraction = primary_He_mass_fraction
    elif element_name == "C":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances("C") - 12) \
                                        * element_weight_table.function_element_weight("C") * Z_0 / Z_solar
    elif element_name == "N":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances("N") - 12) \
                                        * element_weight_table.function_element_weight("N") * Z_0 / Z_solar
    elif element_name == "O":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances("O") - 12) \
                                        * element_weight_table.function_element_weight("O") * Z_0 / Z_solar
    elif element_name == "Ne":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances("Ne") - 12) \
                                        * element_weight_table.function_element_weight("Ne") * Z_0 / Z_solar
    elif element_name == "Mg":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances("Mg") - 12) \
                                        * element_weight_table.function_element_weight("Mg") * Z_0 / Z_solar
    elif element_name == "Si":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances("Si") - 12) \
                                        * element_weight_table.function_element_weight("Si") * Z_0 / Z_solar
    elif element_name == "S":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances("S") - 12) \
                                        * element_weight_table.function_element_weight("S") * Z_0 / Z_solar
    elif element_name == "Ca":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances("Ca") - 12) \
                                        * element_weight_table.function_element_weight("Ca") * Z_0 / Z_solar
    elif element_name == "Fe":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances("Fe") - 12) \
                                        * element_weight_table.function_element_weight("Fe") * Z_0 / Z_solar
    else:
        print("Wrong element name for function_element_mass_primary_fraction")
        element_mass_fraction = None
    return element_mass_fraction
