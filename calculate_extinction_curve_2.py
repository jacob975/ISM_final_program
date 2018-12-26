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
def extinction_ratio(u_source, fitting):
    A_lambda = -2.5 * (umath.log(u_source, 10.0) - umath.log(fitting, 10.0))
    return A_lambda 

class process_Av_ratio():
    def __init__(self, index_table, spectral_table, ratio_table):
        self.index_table = index_table
        self.spectral_table = spectral_table
        self.ratio_table = ratio_table
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
            fitting = ratio_table[i] * SED_template[index_table[i]]
            #print (u_source, fitting)
            u_ratio = extinction_ratio(u_source, fitting)
            ratio[i] = u_ratio.n
            err_ratio[i] = u_ratio.s
        Av_ratio = np.average(ratio, weights = 1/err_ratio)
        err_Av_ratio = np.divide(np.sqrt(np.sum(np.power(err_ratio, 2))), len(source))
        return Av_ratio, err_Av_ratio

#--------------------------------------------
# Main code
if __name__ == "__main__":
    VERBOSE = 0
    # Measure time
    start_time = time.time()
    #-----------------------------------
    # Load arguments
    if len(argv) != 7:
        print("The number of arguments is wrong.")
        print("Usage: calculation_extinction_curve.py   [ source table] [ spectral index table] [ spectral table] [ ratio table] [ Av table] [ Q table]") 
        exit(1)
    source_table_name = argv[1]
    index_table_name = argv[2]
    spectral_table_name = argv[3]
    ratio_table_name = argv[4]
    Av_table_name = argv[5]
    Q_table_name = argv[6]
    #-----------------------------------
    # Load data
    source_table = np.loadtxt(source_table_name)
    index_table = np.loadtxt(index_table_name, dtype = int)
    spectral_table = np.loadtxt(spectral_table_name, dtype = str)
    ratio_table = np.loadtxt(ratio_table_name)
    Av_table = np.loadtxt(Av_table_name)
    Q_table = np.loadtxt(Q_table_name, dtype = str)
    # Load template
    star_template = np.loadtxt("/mazu/users/Jacob975/ISM/20181227final/ISM_final_program/star_koornneef.dat", dtype = str)
    stage = star_template[:,0]
    Spectral_type = star_template[:,1]
    SED_template = np.array(star_template[:,2:], dtype = float)
    band_array = phys_const.band_array
    # Load theoretical extinction curve
    WD55A = np.loadtxt('/mazu/users/Jacob975/ISM/20181227final/ISM_final_program/Synthetic_Extinction_Curves/WD55A.txt')
    WD31A = np.loadtxt('/mazu/users/Jacob975/ISM/20181227final/ISM_final_program/Synthetic_Extinction_Curves/WD31A_test.txt')
    WD31B = np.loadtxt('/mazu/users/Jacob975/ISM/20181227final/ISM_final_program/Synthetic_Extinction_Curves/WD31B.txt')
    # Only take the data with Av > 0
    Av_larger_than_0 = (Av_table[:,0] > 8.0)
    source_table = source_table[Av_larger_than_0]
    index_table = index_table[Av_larger_than_0]
    spectral_table = spectral_table[Av_larger_than_0]
    ratio_table = ratio_table[Av_larger_than_0]
    Av_table = Av_table[Av_larger_than_0]
    Q_table = Q_table[Av_larger_than_0]
    # Take only ratio > 0
    ratio_larger_than_0 = ratio_table > 0.0
    source_table = source_table[ratio_larger_than_0]
    index_table = index_table[ratio_larger_than_0]
    spectral_table = spectral_table[ratio_larger_than_0]
    ratio_table = ratio_table[ratio_larger_than_0]
    Av_table = Av_table[ratio_larger_than_0]
    Q_table = Q_table[ratio_larger_than_0]
    print ('The number of source: {0}'.format(len(source_table)))
    #-----------------------------------
    # Calculate the extinction curve
    stu = process_Av_ratio(index_table, spectral_table, ratio_table)
    Av_ratio = np.zeros(len(band_array))
    err_Av_ratio = np.zeros(len(band_array))
    for i, band in enumerate(band_array):
        print (band)
        source_array = np.transpose(np.array([source_table[:,i], source_table[:,i+8]]))
        SED_template_lite = SED_template[:,i]
        Q_table_lite = Q_table[:,i]
        Av_ratio[i], err_Av_ratio[i] = stu.do(source_array, Q_table_lite, SED_template_lite)
    Av_ratio = Av_ratio/Av_ratio[2]
    print (Av_ratio)
    print (err_Av_ratio)
    fig, axs = plt.subplots(1, 1)
    axs.errorbar(band_array[:,1], \
                Av_ratio, \
                yerr = err_Av_ratio, \
                fmt = 'ro')
    # Normalized by K band
    WD31A_curve = WD31A[:,3]/5.566E-23
    WD55A_curve = WD55A[:,3]/7.547E-23
    WD31B_curve = WD31B[:,3]/5.925E-23
    axs.plot(WD55A[:,0], WD55A_curve, label = 'W&D, $R_{v}$=5.5 A')
    axs.plot(WD31A[:,0], WD31A_curve, linestyle = '--', label = 'W&D, $R_{v}$=3.1A')
    #axs.plot(WD31B[:,0], WD31B_curve, label = 'W&D, $R_{v}$=3.1 B')
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
