# this function returns the element weight


def function_element_weight(element_name):
    # element weight: https://www.lenntech.com/periodic/mass/atomic-mass.htm
    if element_name == "H":
        element_weight = 1.0079
    elif element_name == "He":
        element_weight = 4.0026
    elif element_name == "C":
        element_weight = 12.0107
    elif element_name == "C12":
        element_weight = 12
    elif element_name == "C13":
        element_weight = 13.003355
    elif element_name == "N":
        element_weight = 14.0067
    elif element_name == "N14":
        element_weight = 14.003074
    elif element_name == "N15":
        element_weight = 15.000109
    elif element_name == "O":
        element_weight = 15.9994
    elif element_name == "O16":
        element_weight = 15.99491461956
    elif element_name == "O17":
        element_weight = 16.999131
    elif element_name == "O18":
        element_weight = 17.99916
    elif element_name == "Ne":
        element_weight = 20.1797
    elif element_name == "Ne20":
        element_weight = 19.99244018
    elif element_name == "Ne21":
        element_weight = 20.99384669
    elif element_name == "Ne22":
        element_weight = 21.991385114
    elif element_name == "Na":
        element_weight = 22.989
    elif element_name == "Mg":
        element_weight = 24.305
    elif element_name == "Mg24":
        element_weight = 23.9850419
    elif element_name == "Mg25":
        element_weight = 24.985837
    elif element_name == "Mg26":
        element_weight = 25.9825930
    elif element_name == "Al":
        element_weight = 26.98154
    elif element_name == "Si":
        element_weight = 28.0855
    elif element_name == "Si28":
        element_weight = 27.9769265
    elif element_name == "Si29":
        element_weight = 28.9764947
    elif element_name == "Si30":
        element_weight = 29.9737702
    elif element_name == "S":
        element_weight = 32.065
    elif element_name == "S32":
        element_weight = 31.9720707
    elif element_name == "S33":
        element_weight = 32.9714585
    elif element_name == "S34":
        element_weight = 33.9678668
    elif element_name == "S36":
        element_weight = 35.9670809
    elif element_name == "Ar":
        element_weight = 39.948
    elif element_name == "Ar36":
        element_weight = 35.9675463
    elif element_name == "Ar38":
        element_weight = 37.9627322
    elif element_name == "Ar40":
        element_weight = 39.962383
    elif element_name == "Ca":
        element_weight = 40.078
    elif element_name == "Ca40":
        element_weight = 39.9625909
    elif element_name == "Ca42":
        element_weight = 41.958618
    elif element_name == "Ca43":
        element_weight = 42.958766
    elif element_name == "Ca44":
        element_weight = 43.955482
    elif element_name == "Ca46":
        element_weight = 45.95369
    elif element_name == "Ca48":
        element_weight = 47.9525229
    elif element_name == "Ti":
        element_weight = 47.867
    elif element_name == "Ti46":
        element_weight = 45.952627
    elif element_name == "Ti47":
        element_weight = 46.9517577
    elif element_name == "Ti48":
        element_weight = 47.9479409
    elif element_name == "Ti49":
        element_weight = 48.9478646
    elif element_name == "Ti50":
        element_weight = 48.9478646
    elif element_name == "Cr":
        element_weight = 51.9961
    elif element_name == "Cr50":
        element_weight = 49.946041
    elif element_name == "Cr52":
        element_weight = 51.940505
    elif element_name == "Cr53":
        element_weight = 52.940647
    elif element_name == "Cr54":
        element_weight = 53.938878
    elif element_name == "Mn":
        element_weight = 54.938044
    elif element_name == "Fe":
        element_weight = 55.845
    elif element_name == "Fe54":
        element_weight = 53.939608
    elif element_name == "Fe56":
        element_weight = 55.934936
    elif element_name == "Fe57":
        element_weight = 56.935392
    elif element_name == "Fe58":
        element_weight = 57.933274
    elif element_name == "Ni":
        element_weight = 58.6934
    elif element_name == "Ni58":
        element_weight = 57.935342
    elif element_name == "Ni60":
        element_weight = 59.930785
    elif element_name == "Ni61":
        element_weight = 60.931055
    elif element_name == "Ni62":
        element_weight = 61.928345
    elif element_name == "Ni64":
        element_weight = 63.927966
    elif element_name == "Zn":
        element_weight = 65.37736
    elif element_name == "Zn64":
        element_weight = 63.929142
    elif element_name == "Zn66":
        element_weight = 65.926034
    elif element_name == "Zn67":
        element_weight = 66.927127
    elif element_name == "Zn68":
        element_weight = 67.924844
    elif element_name == "Zn70":
        element_weight = 69.92532
    elif element_name == "Ce":
        element_weight = 140.116
    else:
        print("Wrong element name for function_element_weight")
        element_weight = None
    return element_weight
