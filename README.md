# Monte Carlo integrator for dark matter annihilations

## Introduction
This project aims at the evaluation of <\sigma v> of annihilations/scatterings of dark matter models with a given thermal spectrum, where \sigma is the cross-section and v is velocity of dark matter particles.

## Requirements
Python 3 with numpy and matplotlib.

## How to use it
By modifying integrand.py, you are able to insert the function of your interest as the integrand. The integrand provided in this program is a naive 1/s * v with a Maxwell-Boltzmann velocity distribution, which shows a single peak at s = 5.46.

Importance sampling requires a function g(s) to fit the shape of the integrand. For the exampling integrand we provide an exponential fitting approach in expo_fit.py. It is a piecewise function which consists of 2 exponential functions both peaking at s=5.46. There is also a function in peak25.py that fits a Lorentzian peak at s = 25.

integration.py is a single-channel classical Monte Carlo integration implementation.

integration_multipeaks_optimization.py is a two-channel Monte Carlo integrator. Here the g(s) is constructed by linear combination of 2 base functions. It is useful for integration of functions with another sharp peak. The coefficients of base functions can be automatically modified to minimise the error. Please refer to https://arxiv.org/pdf/hep-ph/9405257.pdf for this self-adjusting approach.
