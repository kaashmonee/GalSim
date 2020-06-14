# Copyright (c) 2012-2019 by the GalSim developers team on GitHub
# https://github.com/GalSim-developers
#
# This file is part of GalSim: The modular galaxy image simulation toolkit.
# https://github.com/GalSim-developers/GalSim
#
# GalSim is free software: redistribution and use in source and binary forms,
# with or without modification, are permitted provided that the following
# conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions, and the disclaimer given in the accompanying LICENSE
#    file.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions, and the disclaimer given in the documentation
#    and/or other materials provided with the distribution.
#
"""
Demo #1

This is the first script in our tutorial about using GalSim in python scripts: examples/demo*.py.
(This file is designed to be viewed in a window 100 characters wide.)

Each of these demo*.py files are designed to be equivalent to the corresponding demo*.yaml file
(or demo*.json -- found in the json directory).  If you are new to python, you should probably
look at those files first as they will probably have a quicker learning curve for you.  Then you
can look through these python scripts, which show how to do the same thing.  Of course, experienced
pythonistas may prefer to start with these scripts and then look at the corresponding YAML files.

To run this script, simply write:

    python demo1.py


This first script is about as simple as it gets.  We draw an image of a single galaxy convolved
with a PSF and write it to disk.  We use a circular Gaussian profile for both the PSF and the
galaxy, and add a constant level of Gaussian noise to the image.

In each demo, we list the new features introduced in that demo file.  These will differ somewhat
between the .py and .yaml (or .json) versions, since the two methods implement things in different
ways.  (demo*.py are python scripts, while demo*.yaml and demo*.json are configuration files.)

New features introduced in this demo:

- obj = galsim.Gaussian(flux, sigma)
- obj = galsim.Convolve([list of objects])
- image = obj.drawImage(scale)
- image.added_flux  (Only present after a drawImage command.)
- noise = galsim.GaussianNoise(sigma)
- image.addNoise(noise)
- image.write(file_name)
- image.FindAdaptiveMom()
"""

import sys
import os
import math
import logging
import galsim
import time
import matplotlib.pyplot as plt
import numpy as np


def main(argv):
    """
    Times that we're keeping track of:
    gal_setup_times = [] => keeps track of galaxy setup times
    gal
    """

    # Initialize the logger
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logger = logging.getLogger("demo7")

    random_seed = 1534225

    gal_flux = 1.e5 # counts
    rng = galsim.BaseDeviate(random_seed + 1)
    gal_r0 = 2.7 # arcseconds
    psf_beta = 5
    pixel_scale = 0.2
    psf_re = 1.0 # arcsec

    nx = 64
    ny = 64

    # Obtaining setup start time to obtain fixed costs for setting up 
    # a particular type of profile. 

    num_gals = 4

    # Creates linearly and logarithmically spaced flux values
    fluxs = np.linspace(1.e3, 1.e7, 15)
    log_fluxs = np.logspace(3, 7, 15)

    # Identify the flux scale used
    # Uncomment one or the other to use one or the other throughout the code.

    flux_scale = fluxs
    # flux_scale = log_fluxs

    # An array to store fluxs. This array is 15x4, where 15 is the number
    # flux values that we are using. There are 4 columns because that is the number of 
    # galaxies.
    setup_times_vary_flux = np.zeros((num_gals, len(flux_scale)))
    gals_flux = []

    for i, gal_flux in enumerate(flux_scale):
        # Storing the setup times for each galaxy
        
        # Using different galaxy profiles
        gal_exp, time_exp_gal = timeit(galsim.Exponential, logger) (half_light_radius=1, flux=gal_flux)
        gal_gauss, time_gauss_gal = timeit(galsim.Gaussian, logger) (half_light_radius=1, flux=gal_flux)
        gal_devauc, time_devauc_gal = timeit(galsim.DeVaucouleurs, logger) (half_light_radius=1, flux=gal_flux)
        gal_sers, time_sers_gal = timeit(galsim.Sersic, logger) (half_light_radius=1, flux=gal_flux, n=2.5)

        # Store the generated galaxy for each flux value
        gals_flux.append([gal_exp, gal_gauss, gal_devauc, gal_sers])

        # Adding all the times to the setup time list
        setup_times_vary_flux[0, i] = time_exp_gal
        setup_times_vary_flux[1, i] = time_gauss_gal
        setup_times_vary_flux[2, i] = time_devauc_gal
        setup_times_vary_flux[3, i] = time_sers_gal


    # Plotting setup times
    galaxy_names = ["Exponential", "Gaussian", "DeVaucouleurs", "Sersic"]

    # Plotting setup times
    plt.title("Setup Time vs. Flux")
    plt.xlabel("Flux")
    plt.ylabel("Setup Time (s)")
    plt.plot(flux_scale, setup_times_vary_flux[0], label=galaxy_names[0])
    plt.plot(flux_scale, setup_times_vary_flux[1], label=galaxy_names[1])
    plt.plot(flux_scale, setup_times_vary_flux[2], label=galaxy_names[2])
    plt.plot(flux_scale, setup_times_vary_flux[3], label=galaxy_names[3])
    plt.legend()
    plt.show()
    plt.figure()

    # Not shearing the galaxy for right now
    
    # Define the Moffat PSF
    moffat_psf, moffat_psf_time = timeit(galsim.Moffat, logger) (beta=psf_beta, flux=1., half_light_radius=psf_re)

    # final profile (these are galsim.Convolution objects)
    finals = []

    # Stores the time it takes to do a convolution with various galaxies at 
    # different fluxs. This is a 15x4 array where at each galaxy generated at a 
    # given flux level, we compute the time required to do a convolution with the Moffat PSF.
    convolution_times = np.zeros((num_gals, len(flux_scale)))

    # For each kind of galaxy, we time the convolution
    for flux_ind, flux_val in enumerate(flux_scale):

        for gals in gals_flux:
            finals_at_flux = []

            for gal_ind, gal in enumerate(gals):

                cnvl_img, time = timeit(galsim.Convolve, logger)([gal, moffat_psf])
                convolution_times[gal_ind, flux_ind] = time
                finals_at_flux.append(cnvl_img)

        finals.append(finals_at_flux)


    # Plotting results...
    plt.title("Time to Convolve with Moffat vs. Flux")
    plt.xlabel("Flux")
    plt.ylabel("Time for Convolution")
    plt.plot(flux_scale, convolution_times[0], label=galaxy_names[0])
    plt.plot(flux_scale, convolution_times[1], label=galaxy_names[1])
    plt.plot(flux_scale, convolution_times[2], label=galaxy_names[2])
    plt.plot(flux_scale, convolution_times[3], label=galaxy_names[3])
    plt.legend()
    plt.show()
    plt.figure()

    # Obtaining the fixed setup costs for a given type of profile
    # Currently, we're only doing Exponential galaxy profiles with Moffat psfs

    # create the image
    images = []
    draw_image_times = []
    for final in finals:
        image = galsim.ImageF(2*nx+2, ny, scale=pixel_scale)
        phot_image = image[galsim.BoundsI(nx+3, 2*nx+2, 1, ny)]

        images.append(image) # Not sure why I need this at the moment. Delete if unnecessary

        img, time = timeit(final.drawImage, logger) (phot_image, method="phot", rng=rng)
        draw_image_times.append(time)

    plt.title("Time to Draw w/Photon Shooting vs. Galaxy")
    plt.xlabel("Galaxy")
    plt.ylabel("Time (Duration)")
    plt.bar(galaxy_names, draw_image_times)
    plt.show()




def timeit(func, logger):
    """
    Takes in a function func and runs and times it on the arguments in *args
    Then outputs the output of func(*args) and the time taken.
    """
    def timeit_wrapper(*args, **kwargs):
        tstart = time.time()
        res = func(*args, **kwargs)
        tend = time.time()
        duration = tend-tstart
        logger.info("Time taken: %f", duration)
        return res, duration

    return timeit_wrapper





if __name__ == "__main__":
    main(sys.argv)
