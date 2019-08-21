import numpy as np
from matplotlib import pyplot as plt
import integrand as inte

peak_0 = 2.0+2.0*np.sqrt(3.0)
height = 0.1295
k2=0.3
k1=1.0                                                                               # parameters for exponential fitting

def g1(s,height=height,k1=k1,k2=k2,peak_0=peak_0):
    """
    g(s). Here it is 2 connected exponential functions with the same peak at peak_0
    s: ndArray or float, s = (p1+p1)^2
    """
    if (isinstance(s,float)):                                                        # Datatype test 
        if(s<peak_0):                                                                # Piecewise function
            return height*np.exp(k1*(s-peak_0))                                      # Left side of the peak
        else:
            return height*np.exp(k2*(peak_0-s))                                      # Right side of the peak
    else:
        ret = np.zeros(len(s),dtype=float)
        ret[s<peak_0] = height*np.exp(k1*(s[s<peak_0]-peak_0))
        ret[s>=peak_0] = height*np.exp(k2*(peak_0-s[s>=peak_0]))
        return ret

def rho1(s,height=height,k1=k1,k2=k2,peak_0=peak_0):
    """
    rho(s). It is essentially indefinite integral of g(s)
    s: ndArray or float, s = (p1+p2)^2
    """
    if(isinstance(s,float)):
        if (s<peak_0):
            return height/k1 * np.exp(k1*(s-peak_0))
        else:
            return height/k1 - height/k2 * np.exp(k2*(peak_0-s)) + height/k2
    else:
        ret = np.zeros(len(s),dtype=float)
        ret[s<peak_0] = height/k1 * np.exp(k1*(s[s<peak_0]-peak_0))
        ret[s>=peak_0] = height/k1 - height/k2 * np.exp(k2*(peak_0-s)) + height/k2
        return ret

def s1(rho,height=height,k1=k1,k2=k2,peak_0=peak_0):
    """
    s(rho). The inverse function of rho(s). It is used to generate importance sampling on s.
    rho: ndArray or float, rho = integral of g(s)
    """
    crit = height/k1                                                                        # value of rho at the peak
    if(isinstance(rho,float)):                                                              # Piecewise function
        if(rho < crit):                                                                    
            return np.log(k1*rho/height)/k1 + peak_0
        else:
            return -np.log(-k2*(rho-crit)/height)/k2 + peak_0
    else:
        ret = np.zeros(len(rho),dtype=float)                                                # The array of s values
        ret[rho < crit] = np.log(k1*rho[rho<crit]/height)/k1 + peak_0
        ret[rho > crit] = -np.log(-k2*(rho[rho>crit]-crit-height/k2)/height)/k2 + peak_0
        return ret

if(__name__=='__main__'):
    print(rho1(2.0))
    ss = np.arange(4.0,50.0,0.01)

    plt.plot(ss,inte.f(ss),color='black')
    plt.plot(ss,g1(ss),color='red')
    plt.show()