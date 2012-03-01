// This just includes the relevant header files from the galsim directory.

#ifndef GALSIM_H
#define GALSIM_H

// The basic SBProfile stuff:
#include "galsim/SBProfile.h"
#include "galsim/SBDeconvolve.h"
#include "galsim/SBParse.h"
#include "galsim/SBPixel.h"

// An interface for dealing with images (FITSImage is deprecated)
#include "galsim/Image.h"
#include "galsim/FITSImage.h"

// FFT's
#include "galsim/FFT.h"

// Noise stuff
#include "galsim/Random.h"
#include "galsim/Poisson.h"

// An integration package by Mike Jarvis
#include "galsim/integ/Int.h"

// Adaptive moments code by Hirata, Seljak, and Mandelbaum
#include "galsim/hsm/psfcorr.h"

#endif
