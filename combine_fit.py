import numpy as np
from matplotlib import pyplot as plt
import integrand as inte
import cauchy_fit
import expo_fit

peak_0 = 2.0+2.0*np.sqrt(3.0)
height = 0.1295
b1 = 0.7
b2 = 0.15

def g1(s):
    '''for zero_0<x<peak_0'''
    return cauchy_fit.g1(s)

def pho1(s):
    return cauchy_fit.pho1(s)

def s1(pho):
    return cauchy_fit.s1(pho)

def g2(s):
    return expo_fit.g2(s)

def pho2(s):
    return expo_fit.pho2(s)

def s2(pho):
    return expo_fit.s2(pho)

if(__name__=='__main__'):
    ss = np.arange(4.0,50.0,0.01)
    ss1 = np.arange(4.0,peak_0,0.01)
    ss2 = np.arange(peak_0,50.0,0.01)
    
    plt.plot(ss,inte.f(ss),color='black')
    plt.plot(ss1,g1(ss1),color='green')
    plt.plot(ss2,g2(ss2),color='red')
    plt.savefig('combine_fit.png')
    plt.show()

    plt.figure()
    plt.plot(ss1, inte.f(ss1)/g1(ss1))
    plt.plot(ss2, inte.f(ss2)/g2(ss2))
    plt.yscale('log')
    plt.savefig('combine_fit_goodness.png')
    plt.show()