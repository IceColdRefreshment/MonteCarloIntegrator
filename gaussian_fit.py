import numpy as np
from matplotlib import pyplot as plt
import integrand as inte

peak_0 = 2.0+2.0*np.sqrt(3.0)
height = 0.129
mean = peak_0

def g1(s):
    '''for zero_0<x<peak_0'''
    sigma = 1.3
    return height * np.exp(-0.5 * ((s-mean)/sigma)**2)
def g2(s):
    '''for x>peak_0'''
    sigma = 2.4
    return height * np.exp(-0.5 * ((s-mean)/sigma)**2)

'''def MonteCarloIntegrate():
    xs = np.arange(start,end,interval)
''' 

if(__name__=='__main__'):
    ss = np.arange(4.0,50.0,0.01)
    ss1 = np.arange(4.0,peak_0,0.01)
    ss2 = np.arange(peak_0,50.0,0.01)
    
    plt.plot(ss,inte.f(ss),color='black')
    plt.plot(ss1,g1(ss1),color='green')
    plt.plot(ss2,g2(ss2),color='red')
    plt.savefig('gaussian_fit.png')
    plt.show()

    plt.figure()
    plt.plot(ss1, inte.f(ss1)/g1(ss1))
    plt.plot(ss2, inte.f(ss2)/g2(ss2))
    plt.yscale('log')
    plt.savefig('gaussian_fit_goodness.png')
    plt.show()