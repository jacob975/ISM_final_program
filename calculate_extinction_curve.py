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

def extinction_ratio(source, fitting, Av):
    result_value = np.zeros(len(fitting))
    result_error = np.zeros(len(fitting))
    for i in range(8):
        u_source = ufloat(source[i], source[i+8])
        u_Av = ufloat(Av[0], Av[1])
        print (u_source)
        print (fitting[i])
        print (u_Av)
        result = -2.5 * (umath.log(u_source, 10.0) - umath.log(fitting[i], 10.0))
        result_value[i] = result.n
        result_error[i] = result.s
    return result_value, result_error

#--------------------------------------------
# Main code
if __name__ == "__main__":
    VERBOSE = 0
    # Measure time
    start_time = time.time()
    #-----------------------------------
    # Load arguments
    if len(argv) != 6:
        print("The number of arguments is wrong.")
        print("Usage: calculation_extinction_curve.py [source table] [spectral index table] [spectral table] [ratio table] [ Av table]") 
        exit(1)
    source_table_name = argv[1]
    index_table_name = argv[2]
    spectral_table_name = argv[3]
    ratio_table_name = argv[4]
    Av_table_name = argv[5]
    #-----------------------------------
    # Load data
    source_table = np.loadtxt(source_table_name)
    index_table = np.loadtxt(index_table_name, dtype = int)
    spectral_table = np.loadtxt(spectral_table_name, dtype = str)
    ratio_table = np.loadtxt(ratio_table_name)
    Av_table = np.loadtxt(Av_table_name)
    # Load template
    star_template = np.loadtxt("/mazu/users/Jacob975/ISM/20181227final/ISM_final_program/star_koornneef.dat", dtype = str)
    stage = star_template[:,0]
    Spectral_type = star_template[:,1]
    SED_template = np.array(star_template[:,2:], dtype = float)
    band_array = phys_const.band_array
    # Calculate the extinction curve
    for i, source in enumerate(source_table[100:]):
        selected_fitted = ratio_table[i] * SED_template[index_table[i]] 
        # Plot the result
        fig, axs = plt.subplots(1, 1)
        ratio_value, ratio_error = extinction_ratio(source, selected_fitted, Av_table[i])
        axs.errorbar(band_array[:,1], \
                    ratio_value, \
                    yerr = ratio_error, \
                    fmt = 'ro')
        axs.set_ylabel('$A_{\lambda}$/$A_{v}$')
        ax2 = axs.twinx()
        ax2.plot(band_array[:,1], source[:8], label = 'observation data', c = 'b')
        ax2.errorbar(band_array[:,1], source[:8], yerr = source[8:], fmt = 'bo')
        ax2.plot(band_array[:,1], selected_fitted, label = 'fitting result', c = 'g')
        ax2.scatter(band_array[:,1], selected_fitted, c = 'g')
        ax2.set_ylabel('Flux(mJy)')
        plt.xscale('log')
        plt.xlabel('wavelength($\mu$m)')
        plt.grid(True)
        plt.legend()
        plt.show()
    #-----------------------------------
    # Measure time
    elapsed_time = time.time() - start_time
    print("Exiting Main Program, spending ", elapsed_time, "seconds.")
