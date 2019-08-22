import numpy as np
from matplotlib import pyplot as plt
import integrand_multipeaks as inte

peak = 30.0
height = 0.00002
b1 = 7.0

def g1(s):
    '''for zero_0<x<peak'''
    return height / ( b1 * (s-peak)**2 + 1.0 )

def rho1(s):
    return height / np.sqrt(b1) * np.arctan((s-peak)*np.sqrt(b1))

def s1(rho):
    return peak + 1.0/np.sqrt(b1) * np.tan(np.sqrt(b1) * rho / height)

if(__name__=='__main__'):
    ss = np.arange(4.0,50.0,0.01)
    
    plt.plot(ss,inte.f(ss),color='black')
    plt.plot(ss,g1(ss),color='green')
    plt.savefig('mp_fit.png')
    plt.show()

    plt.figure()
    plt.plot(ss, inte.f(ss)/g1(ss))
    plt.plot(ss,g1(ss),color='green')
    plt.yscale('log')
    plt.savefig('mp_fit_goodness.png')
    plt.show()