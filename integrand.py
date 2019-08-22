import numpy as np
from matplotlib import pyplot as plt

def f(s):
    m=1.0
    v0=1.0
    #factor0 = np.sqrt(2.0/np.pi)/(v0**3 * m**4)                            # Scalar coefficient
    factor0 = 1.0
    factor1 = (s - 4.0*(m**2)) * np.exp(-(s-4.0*m**2)/(2.0*m**2*v0**2))     # velocity and its distribution
    factor2 = 1.0 / s                                                       # matrix element
    return factor0*factor1*factor2

if __name__=="__main__":
    ss = np.arange(4.0,100.0,0.1)
    #print(ss)

    plt.xlabel('s')
    plt.ylabel(r'$\sigma v$')
    ys = f(ss)
    plt.plot(ss,ys)
    plt.savefig('integrand.png')
    plt.show()

    plt.figure()
    plt.plot(ss,ys)
    plt.yscale('log')
    plt.savefig('integrand_log.png')
    plt.show()