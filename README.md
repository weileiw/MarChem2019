# MarChem2019
Code for MarChem2019 paper

Box_model_eigen_2nd.m is the main driver

neglogpost_eigen.m computes the negative of the logarithm of the posterior pdf for the parameters under steady-state assumption.

buildPFD_l.m computes the sinking particle flux divergence operator for large particles.

buildPFD_s.m computes the sinking particle flux divergence operator for small particles.

cost_2nd.m computes the negative of the logrithm of the posterior pdf for the parameters including slowly-decaying eigen modes.

Box_model_eigen_2nd.m | +--> neglogpost_eigen.m | +--> cost_2nd.m | | | +-->[buildPFD_l.m/buildPFD_s.m]

unitlity scripts:

nsnew.m C.T. Kelly's Newton Armijo solver

mfactor.m Timothy A. Davis' LINFACTOR VERSION 1.1.0, Nov 1, 2007 Copyright 2007, Timothy A. Davis, University of Florida

d0.m makes a sparse diagonal matrix given a 3d field
