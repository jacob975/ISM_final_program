#!/usr/bin/python3
'''
Abstract:
    This is a program for saving the constant I use. 
Usage:
    1. Find a file you like.
    2. import phys_const.py
Editor:
    Jacob975

##################################
#   Python3                      #
#   This code is made in python3 #
##################################

20181220
####################################
update log
20181220 version alpha 1:
    1. The code works.
'''
import numpy as np

# The band and the corresponding wavelength.
band = [ [ 'J',   1.248], 
         [ 'H',   1.631],
         [ 'Ks',  2.201], 
         [ 'IR1', 3.6  ], 
         [ 'IR2', 4.5  ], 
         [ 'IR3', 5.8  ], 
         [ 'IR4', 8.0  ], 
         [ 'MP1', 24.0 ]]
band_array = np.array(band, dtype = object)

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

