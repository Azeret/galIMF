# this function returns the element weight


def function_element_weight(element_name):
    # element weight: https://www.lenntech.com/periodic/mass/atomic-mass.htm
    if element_name == "H":
        element_weight = 1.0079
    elif element_name == "He":
        element_weight = 4.0026
    elif element_name == "C":
        element_weight = 12.0107
    elif element_name == "N":
        element_weight = 14.0067
    elif element_name == "O":
        element_weight = 15.9994
    elif element_name == "Ne":
        element_weight = 20.1797
    elif element_name == "Mg":
        element_weight = 24.305
    elif element_name == "Si":
        element_weight = 28.0855
    elif element_name == "S":
        element_weight = 32.065
    elif element_name == "Ca":
        element_weight = 40.078
    elif element_name == "Fe":
        element_weight = 55.845
    else:
        print("Wrong element name for function_element_weight")
        element_weight = None
    return element_weight
