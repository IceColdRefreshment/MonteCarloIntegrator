# Monte Carlo integrator for dark matter annihilations

## Introduction
This project aims at the evaluation of <\sigma v> of annihilations/scatterings of dark matter models with a given thermal spectrum, where \sigma is the cross-section and v is velocity of dark matter particles.

## Requirements
Python with numpy and matplotlib.

## How to use it
By modifying integrand.py, you are able to insert the function of your interest as the integrand. The integrand.py provided in this program is a naive 1/s * v with a Maxwell-Boltzmann velocity distribution, which shows a single peak at s = 5.46. The integrand_multipeaks.py file contains an integrand with an additional peak at s = 25.

Importance sampling requires a function g(s) to fit the shape of the integrand. We provide an exponential fitting in expo_fit.py for the example integrand.py. It is a piecewise function which consists of 2 exponential functions both peaking at s=5.46. Note that an integral of g(s), rho(s), and its inverse s(rho) are also required to complete the importance sampling. 

integration.py is a single-channel classical Monte Carlo integration implementation.

integration_multipeaks_optimization.py is a two-channel Monte Carlo integrator. Here the g(s) is constructed by linear combination of 2 base functions g1(s) and g2(s). It is useful for integration of functions with another sharp peak. The coefficients of base functions can be automatically modified in order to minimise the error. Please refer to https://arxiv.org/pdf/hep-ph/9405257.pdf for this self-adjusting approach. In the example, we use expo_fit.py as g1(s) and peak25.py as g2(s) and successfully integrate integrand_multipeaks.py.
