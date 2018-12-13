#!/usr/bin/python3
'''
Abstract:
    This is a program for showing planck func 
Usage:
    planck_func.py
Editor:
    Jacob975

##################################
#   Python3                      #
#   This code is made in python3 #
##################################

20181206
####################################
update log
20181206 version alpha 1:
    1. The code works.
'''
import matplotlib.pyplot as plt
import time
import numpy as np
from scipy import optimize
from scipy.stats import chisquare
from sys import argv

# Constants in cgs unit
k = 1.381e-16
h = 6.626e-27
c = 2.998e10
m_H = 1.674e-24
mu = 2.8
pixel_size = 0.0035/3600 # degree
pc = 3.08567758e18

#-------------------------------------------------
# Planck function
#
#               2hc^2               1
#  B_l(l, T) = ------- ----------------------------
#                l^5    e^(hc/(l * k_B * T)) - 1
#

def planck(wavelength_mu, T, A):
    wav  = wavelength_mu*1e-4
    a = 2.0*h*c**2
    b = np.array(np.divide(h*c, (wav*k*T)), dtype = float)
    intensity = np.divide(a, wav**5 * (np.exp(b) - 1.0) )
    return np.array(A * intensity, dtype = float)

# Find the best template of this observation data.
def match_process(source, templates):
    # Initialize 
    ref_chi = 999.0
    ref_p = 0.0
    index = 0
    match_template = None
    # Take each row, calculate the chi-square value.
    for i in range(len(templates)):
        selected_template = templates[i,:3]
        chi, p, ratio = minimum_chi(source, selected_template)
        # If better, take it
        if ref_chi > chi:
            ref_chi = chi
            ref_p = p
            index = i
            ref_ratio = ratio
    # Pick the minimum
    return ref_chi, ref_p, index, ref_ratio 

# Find the best fit of this template
def minimum_chi(source, template):
    # Initialize
    x_linspace = np.logspace(-0.3, 0.3, 30)
    ratios = x_linspace * source[2]
    table = np.outer(ratios, template)
    chi, p = chisquare(source, f_exp = table, axis = 1)
    index = np.argmin(chi)
    ref_chi = chi[index]
    ref_p = p[index]
    ref_ratio = ratios[index]
    return ref_chi, ref_p, ref_ratio 
#--------------------------------------------
# Main code
if __name__ == "__main__":
    VERBOSE = 0
    # Measure time
    start_time = time.time()
    #-----------------------------------
    # Load argv
    if len(argv) != 2:
        print ('The number of arguments is wrong.')
        print ("Usage: planck_func.py [sed table]")
        print ("Example: planck_func.py source_sed_MaxLoss0.txt")
        exit()
    sed_table_name = argv[1]
    #-----------------------------------
    # Load data in J, H, K bands.
    sed_data = np.loadtxt(sed_table_name)
    JHK_sed_data = sed_data[:,:3]
    # fitting planck func with near-infrared data and given temprature.
    band = [ [ 'J',   1.248], 
             [ 'H',   1.631],
             [ 'Ks',  2.201], 
             [ 'IR1', 3.6  ], 
             [ 'IR2', 4.5  ], 
             [ 'IR3', 5.8  ], 
             [ 'IR4', 8.0  ], 
             [ 'MP1', 24.0 ]]
    band_array = np.array(band, dtype = object)
    star_template = np.loadtxt("/mazu/users/Jacob975/ISM/20181227final/ISM_final_program/star_koornneef.dat", dtype = str)
    stage = star_template[:,0]
    Spectral_type = star_template[:,1]
    SED_template = np.array(star_template[:,2:], dtype = float)
    match_index = np.zeros(len(sed_data))
    match_ratio = np.zeros(len(sed_data))
    num_data = len(JHK_sed_data)
    for i, source in enumerate(JHK_sed_data):
        # Fit the sed data to find the temperature.
        chi, p, template_index, ratio = match_process(source, SED_template)
        match_index[i] = template_index
        match_ratio[i] = ratio
        '''
        # Show the plots  
        plt.plot(band_array[:3,1], source, label = 'observation data', c = 'b')
        plt.scatter(band_array[:3,1], source, c = 'b')
        plt.plot(band_array[:3,1], match_template, label = 'fitted template', c = 'r')
        plt.scatter(band_array[:3,1], match_template, c = 'r')
        plt.xlabel('wavelength $\lambda(\mu)$')
        plt.ylabel('Flux(mJy)')
        plt.xscale("log")
        plt.legend()
        plt.show()
        '''
        if i%100 == 0:
            print ("Progress: {0}/{1}".format(i+1, num_data))
    # Save the result
    match_spectral = Spectral_type[match_index]
    np.savetxt('match_spectral.txt', match_spectral)
    np.savetxt('match_index.txt', match_index)
    np.savetxt('match_ratio.txt', match_ratio)
    #-----------------------------------
    # Measure time
    elapsed_time = time.time() - start_time
    print("Exiting Main Program, spending ", elapsed_time, "seconds.")
