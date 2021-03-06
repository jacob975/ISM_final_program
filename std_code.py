#!/usr/bin/python3
'''
Abstract:
    This is a program for showing sources on image. 
Usage:
    std_code.py [coord table] [fits image file]

    coord table is a txt file with the following form.
        [[RA, DEC], 
         [RA, DEC],
         ...
        ]
    fits image file is an image file with WCS
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
from sys import argv
import numpy as np
from astropy.io import fits as pyfits
from astropy import wcs
from astropy.coordinates import *

#--------------------------------------------
# Main code
if __name__ == "__main__":
    VERBOSE = 0
    # Measure time
    start_time = time.time()
    #-----------------------------------
    # Load arguments
    if len(argv) != 3:
        print("The number of arguments is wrong.")
        print("Usage: show_sources.py [coord table] [fits image file]")
        exit(1)
    coord_table_name = argv[1]
    image_name = argv[2]
    #-----------------------------------
    # Load data
    world_coord = np.loadtxt(coord_table_name)
    image = pyfits.getdata(image_name) 
    header = pyfits.getheader(image_name)
    #-----------------------------------
    # Convert WCS coord to pixel coord
    w = wcs.WCS(header)
    pixel_coord = w.wcs_world2pix(world_coord, 1)
    # Plot and show
    fig = plt.figure(figsize = (8, 8))
    plt.subplot(111, projection = w)
    plt.title("Source on {0}".format(image_name))
    plt_image = plt.imshow(image)
    plt.colorbar()
    plt.scatter(pixel_coord[:,0], pixel_coord[:,1], s= 2, c= 'r' )
    plt.show()
    #-----------------------------------
    # Measure time
    elapsed_time = time.time() - start_time
    print("Exiting Main Program, spending ", elapsed_time, "seconds.")
