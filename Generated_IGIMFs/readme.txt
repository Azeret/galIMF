Here saves all the galaxy-wide IMFs for galaxy star formation epochs
as a function of their metallicity, [Z/X], and star formation rate, SFR.

igimf_file_name =
"igimf_SFR_{}_Fe_over_H_{}".format(round(math.log(SFR, 10) * 100000), round(Z_over_X * 100000))
