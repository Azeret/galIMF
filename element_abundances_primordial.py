import element_weight_table, element_abundances_solar

H_weight = element_weight_table.function_element_weight("H")

primary_He_mass_fraction = 0.247
primary_H_mass_fraction_roughly = 1 - primary_He_mass_fraction  # Corrected in below
primary_D_mass_fraction = primary_H_mass_fraction_roughly * 2.58 * 10**-5
primary_He3_mass_fraction = primary_H_mass_fraction_roughly * 10**-4
primary_L_mass_fraction = primary_H_mass_fraction_roughly * 5 * 10**-10
# Reference: Cyburt+ 2016, Big bang nucleosynthesis: Present status, DOI: 10.1103/RevModPhys.88.015004

Z_0 = 10**-6
primary_H_mass_fraction = 1 - primary_He_mass_fraction - primary_D_mass_fraction - primary_He3_mass_fraction\
                          - primary_L_mass_fraction - Z_0


def function_element_mass_primary_fraction(solar_abu_reference_name, element_name, Z_0, Z_solar):
    if element_name == "H":
        element_mass_fraction = primary_H_mass_fraction
    elif element_name == "He":
        element_mass_fraction = primary_He_mass_fraction
    elif element_name == "C":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances(solar_abu_reference_name, "C") - 12) \
                                        * element_weight_table.function_element_weight("C") * Z_0 / Z_solar
    elif element_name == "N":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances(solar_abu_reference_name, "N") - 12) \
                                        * element_weight_table.function_element_weight("N") * Z_0 / Z_solar
    elif element_name == "O":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances(solar_abu_reference_name, "O") - 12) \
                                        * element_weight_table.function_element_weight("O") * Z_0 / Z_solar
    elif element_name == "Ne":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances(solar_abu_reference_name, "Ne") - 12) \
                                        * element_weight_table.function_element_weight("Ne") * Z_0 / Z_solar
    elif element_name == "Mg":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances(solar_abu_reference_name, "Mg") - 12) \
                                        * element_weight_table.function_element_weight("Mg") * Z_0 / Z_solar
    elif element_name == "Si":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances(solar_abu_reference_name, "Si") - 12) \
                                        * element_weight_table.function_element_weight("Si") * Z_0 / Z_solar
    elif element_name == "S":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances(solar_abu_reference_name, "S") - 12) \
                                        * element_weight_table.function_element_weight("S") * Z_0 / Z_solar
    elif element_name == "Ca":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances(solar_abu_reference_name, "Ca") - 12) \
                                        * element_weight_table.function_element_weight("Ca") * Z_0 / Z_solar
    elif element_name == "Fe":
        element_mass_fraction = primary_H_mass_fraction / H_weight\
                                        * 10**(element_abundances_solar.function_solar_element_abundances(solar_abu_reference_name, "Fe") - 12) \
                                        * element_weight_table.function_element_weight("Fe") * Z_0 / Z_solar
    else:
        print("Wrong element name for function_element_mass_primary_fraction")
        element_mass_fraction = None
    return element_mass_fraction
