#!/usr/bin/python3
'''
Abstract:
    This is a program for calculating extinction curves
Usage:
    calculation_extinction_curve.py [source table] [spectral index table] [spectral table] [ratio table]
Editor:
    Jacob975

##################################
#   Python3                      #
#   This code is made in python3 #
##################################

20181213
####################################
update log
20181213 version alpha 1:
    1. The code works.
'''
import matplotlib.pyplot as plt
import time
from sys import argv
import numpy as np
import phys_const
from uncertainties import ufloat, umath

# Calculate the ratio of A_lambda and Av
def extinction_ratio(u_source, u_fitting, u_Av):
    A_lambda = -2.5 * (umath.log(u_source, 10.0) - umath.log(u_fitting, 10.0))
    ratio = A_lambda/u_Av
    return A_lambda 

class process_Av_ratio():
    def __init__(self, index_table, spectral_table, ratio_table, Av_table, chi_table):
        self.index_table = index_table
        self.spectral_table = spectral_table
        self.ratio_table = ratio_table
        self.Av_table = Av_table
        self.chi_table = chi_table
        self.sqrt_chi_table = np.sqrt(chi_table)
        return
    def do(self, source_array, Q_array, SED_template):
        # Wipe out all Upper limit
        non_UL = Q_array != 'U'
        source_array = source_array[non_UL]
        ratio_table = self.ratio_table[non_UL]
        index_tabel = self.index_table[non_UL]
        ratio = np.zeros(len(source_array))
        err_ratio = np.zeros(len(source_array))
        for i, source in enumerate(source_array):
            u_source = ufloat(source[0], source[1])
            u_Av = ufloat(self.Av_table[i,0], self.Av_table[i,1])
            fitting = ratio_table[i] * SED_template[index_table[i]]
            err_fitting = fitting * self.sqrt_chi_table[i]
            u_fitting = ufloat(fitting, err_fitting)
            u_ratio = extinction_ratio(u_source, u_fitting, u_Av)
            ratio[i] = u_ratio.n
            err_ratio[i] = u_ratio.s
        Av_ratio = np.average(ratio, weights = 1/err_ratio)
        var_Av_ratio = np.sum(np.divide(np.power(ratio - Av_ratio, 2), err_ratio))/np.sum(1/err_ratio)
        err_Av_ratio = np.sqrt(var_Av_ratio)
        return Av_ratio, err_Av_ratio

#--------------------------------------------
# Main code
if __name__ == "__main__":
    VERBOSE = 0
    # Measure time
    start_time = time.time()
    #-----------------------------------
    # Load arguments
    if len(argv) != 8:
        print("The number of arguments is wrong.")
        print("Usage: calculation_extinction_curve.py   [ source table] [ spectral index table] [ spectral table] [ ratio table] [ Av table] [ Q table] [chi table]") 
        exit(1)
    source_table_name = argv[1]
    index_table_name = argv[2]
    spectral_table_name = argv[3]
    ratio_table_name = argv[4]
    Av_table_name = argv[5]
    Q_table_name = argv[6]
    chi_table_name = argv[7]
    #-----------------------------------
    # Load data
    source_table = np.loadtxt(source_table_name)
    index_table = np.loadtxt(index_table_name, dtype = int)
    spectral_table = np.loadtxt(spectral_table_name, dtype = str)
    ratio_table = np.loadtxt(ratio_table_name)
    Av_table = np.loadtxt(Av_table_name)
    Q_table = np.loadtxt(Q_table_name, dtype = str)
    chi_table = np.loadtxt(chi_table_name)
    # Load template
    star_template = np.loadtxt("/mazu/users/Jacob975/ISM/20181227final/ISM_final_program/star_koornneef.dat", dtype = str)
    stage = star_template[:,0]
    Spectral_type = star_template[:,1]
    SED_template = np.array(star_template[:,2:], dtype = float)
    band_array = phys_const.band_array
    # Load theoretical extinction curve
    WD55A = np.loadtxt('/mazu/users/Jacob975/ISM/20181227final/ISM_final_program/Synthetic_Extinction_Curves/WD55A.txt')
    WD31A = np.loadtxt('/mazu/users/Jacob975/ISM/20181227final/ISM_final_program/Synthetic_Extinction_Curves/WD31A_test.txt')
    WD55B = np.loadtxt('/mazu/users/Jacob975/ISM/20181227final/ISM_final_program/Synthetic_Extinction_Curves/WD55B.txt')
    WD31B = np.loadtxt('/mazu/users/Jacob975/ISM/20181227final/ISM_final_program/Synthetic_Extinction_Curves/WD31B.txt')
    # Only take the data with Av > 0
    Av_larger_than_0 = (Av_table[:,0] > 12.0)
    source_table = source_table[Av_larger_than_0]
    index_table = index_table[Av_larger_than_0]
    spectral_table = spectral_table[Av_larger_than_0]
    ratio_table = ratio_table[Av_larger_than_0]
    Av_table = Av_table[Av_larger_than_0]
    Q_table = Q_table[Av_larger_than_0]
    chi_table = chi_table[Av_larger_than_0]
    # Take only ratio > 0
    ratio_larger_than_0 = ratio_table > 0.0
    source_table = source_table[ratio_larger_than_0]
    index_table = index_table[ratio_larger_than_0]
    spectral_table = spectral_table[ratio_larger_than_0]
    ratio_table = ratio_table[ratio_larger_than_0]
    Av_table = Av_table[ratio_larger_than_0]
    Q_table = Q_table[ratio_larger_than_0]
    chi_table = chi_table[ratio_larger_than_0]
    print ('The number of source: {0}'.format(len(source_table)))
    #-----------------------------------
    # Calculate the extinction curve
    stu = process_Av_ratio(index_table, spectral_table, ratio_table, Av_table, chi_table)
    Av_ratio = np.zeros(len(band_array))
    err_Av_ratio = np.zeros(len(band_array))
    for i, band in enumerate(band_array):
        print (band)
        source_array = np.transpose(np.array([source_table[:,i], source_table[:,i+8]]))
        SED_template_lite = SED_template[:,i]
        Q_table_lite = Q_table[:,i]
        Av_ratio[i], err_Av_ratio[i] = stu.do(source_array, Q_table_lite, SED_template_lite)
    Av_ratio = Av_ratio/Av_ratio[2]
    err_Av_ratio = err_Av_ratio/Av_ratio[2]
    print (Av_ratio)
    print (err_Av_ratio)
    fig, axs = plt.subplots(1, 1)
    axs.errorbar(band_array[:,1], \
                Av_ratio, \
                yerr = err_Av_ratio, \
                fmt = 'ro',
                markersize = 2)
    # Normalized by K band
    WD31A_curve = WD31A[:,3]/5.566E-23
    WD55A_curve = WD55A[:,3]/7.547E-23
    WD31B_curve = WD31B[:,3]/5.925E-23
    WD55B_curve = WD55B[:,3]/5.331E-23 
    #axs.plot(WD55A[:,0], WD55A_curve, label = 'W&D, $R_{v}$=5.5A')
    #axs.plot(WD31A[:,0], WD31A_curve, linestyle = '--', label = 'W&D, $R_{v}$=3.1A')
    axs.plot(WD55B[:,0], WD55B_curve, label = 'W&D(2001), $R_{v}$=5.5B')
    axs.plot(WD31B[:,0], WD31B_curve, label = 'W&D(2001), $R_{v}$=3.1B')
    axs.set_ylabel('$A_{\lambda}$/$A_{K}$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1, 25)
    plt.ylim(1e-2, 1e1)
    plt.xlabel('wavelength($\mu$m)')
    plt.grid(True)
    plt.legend()
    plt.savefig('predicted_extinction_curve.png')
    #-----------------------------------
    # Measure time
    elapsed_time = time.time() - start_time
    print("Exiting Main Program, spending ", elapsed_time, "seconds.")
