import numpy as np
from matplotlib import pyplot as plt
import integrand as inte

peak_0 = 2.0+2.0*np.sqrt(3.0)
height = 0.1295
b1 = 0.7
b2 = 0.15

def g1(s):
    '''for zero_0<x<peak_0'''
    return height / ( b1 * (s-peak_0)**2 + 1.0 )

def pho1(s):
    return height / np.sqrt(b1) * np.arctan((s-peak_0)*np.sqrt(b1))

def s1(pho):
    return peak_0 + 1.0/np.sqrt(b1) * np.tan(np.sqrt(b1) * pho / height)

def g2(s):
    '''for x>peak_0'''
    return height / ( b2 * (s-peak_0)**2 + 1.0 )

def pho2(s):
    return height / np.sqrt(b2) * np.arctan((s-peak_0)*np.sqrt(b2))

def s2(pho):
    return peak_0 + 1.0/np.sqrt(b2) * np.tan(np.sqrt(b2) * pho / height)

if(__name__=='__main__'):
    ss = np.arange(4.0,50.0,0.01)
    ss1 = np.arange(4.0,peak_0,0.01)
    ss2 = np.arange(peak_0,50.0,0.01)
    
    plt.plot(ss,inte.f(ss),color='black')
    plt.plot(ss1,g1(ss1),color='green')
    plt.plot(ss2,g2(ss2),color='red')
    plt.savefig('cauchy0_fit.png')
    plt.show()

    plt.figure()
    plt.plot(ss1, inte.f(ss1)/g1(ss1))
    plt.plot(ss2, inte.f(ss2)/g2(ss2))
    plt.yscale('log')
    plt.savefig('cauchy0_fit_goodness.png')
    plt.show()