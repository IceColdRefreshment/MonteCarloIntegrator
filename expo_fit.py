import numpy as np
from matplotlib import pyplot as plt
import integrand as inte

peak_0 = 2.0+2.0*np.sqrt(3.0)
height = 0.1295
k2=0.3
k1=1.0
def g1(s):
    '''for zero_0<x<peak_0'''
    return height*np.exp(k1*(s-peak_0))

def pho1(s):
    '''pho1(s)'''
    return height/k1 * np.exp(k1*(s-peak_0))

def s1(pho):
    '''s1(pho)'''
    return np.log(k1*pho/height)/k1 + peak_0

def g2(s):
    '''for x>peak_0'''
    return height*np.exp(k2*(peak_0-s))

def pho2(s):
    '''pho2(s)'''
    return -height/k2 * np.exp(k2*(peak_0-s))

def s2(pho):
    '''s2(pho)'''
    return -np.log(-k2*pho/height)/k2 + peak_0

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
    plt.savefig('expo_fit.png')
    plt.show()

    plt.figure()
    plt.plot(ss1, inte.f(ss1)/g1(ss1))
    plt.plot(ss2, inte.f(ss2)/g2(ss2))
    plt.yscale('log')
    plt.savefig('expo_fit_goodness.png')
    plt.show()