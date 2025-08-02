import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import math

# Open the FITS file
filename = "J_A+A_562_A71_tablec3.dat.fits"
with fits.open(filename) as hdul:
    # Access the data in the BinTableHDU
    data = hdul[1].data
    # print(data[0])
    # Temp = []
    # logg = []
    mass = []
    # age = []
    fe_h = []
    o_fe = []
    # mg_fe = []
    si_fe = []
    # ca_fe = []
    e_fe_h = []
    e_o_fe = []
    # e_mg_fe = []
    e_si_fe = []
    # e_ca_fe = []
    # td_d = []
    # td_h = []
    for i in data:
        if i[10] < 1 and i[4]>4.1 and i[34] > 0 and i[30] > 0:
        # if td_h[i] > 1 and 10**td_d[i] < 0.5 and mass[i] < 1 and logg[i]>4.1 and e_si_fe[i]>0 and e_o_fe[i]>0:
            # Temp.append(i[2])
            # logg.append(i[4])
            mass.append(i[10])
            # age.append(i[13])
            fe_h.append(i[16])
            o_fe.append(i[17])
            # mg_fe.append(i[19])
            si_fe.append(i[21])
            # ca_fe.append(i[22])
            e_fe_h.append(i[29])
            e_o_fe.append(i[30])
            # e_mg_fe.append(i[32])
            e_si_fe.append(i[34])
            # e_ca_fe.append(i[35])
            # td_d.append(math.log(i[-10], 10))
            # td_h.append(i[-9])

print(min(mass))
# plt.scatter(Temp, logg, s=3, alpha=1)
# print("len(data)", len(data))  # 714
# print("len(masked)", len(fe_h))  # 268

# import statistics
# mean = statistics.mean(e_fe_h)
# print("Mean:", mean)
# # std_dev = statistics.stdev(fe_h)
# # print("Standard Deviation:", std_dev)
# mean = statistics.mean(e_o_fe)
# print("Mean:", mean)
# mean = statistics.mean(e_si_fe)
# print("Mean:", mean)

# plt.hist(fe_h, bins=30, density=True, alpha=0.5, label="masked2")
# # plt.scatter(fe_h, si_fe, s=3, alpha=0.5, c='r', label="masked2")
# # plt.scatter(fe_h, o_fe, s=3, alpha=0.5, c='r', label="masked")
# plt.legend()
# plt.show()
#
