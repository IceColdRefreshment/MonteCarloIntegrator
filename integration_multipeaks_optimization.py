import numpy as np
from matplotlib import pyplot as plt
import integrand_multipeaks as inte
import expo_fit
import peak25 as peak_fit

start = 4.0
end = 40.0

alpha_1 = 0.5
alpha_2 = 0.5                                                                   # alpha coefficients of two channels

int1_1, int2_1, area_1, norm_factor_1 = 1.0, 1.0, 1.0, 1.0
int1_2, int2_2, area_2, norm_factor_2 = 1.0, 1.0, 1.0, 1.0                      # initialisation

def initialise(fit_1=expo_fit,fit_2=peak_fit,start=start,end=end):
    '''
    initialistion
    fit_1: function, the first function for fitting
    fit_2: function, the second function for fitting
    start, end: floats, range of integration
    '''
    global int1_1, int2_1, area_1, norm_factor_1, int1_2, int2_2, area_2, norm_factor_2
    
    int1_1, int2_1 = fit_1.rho1(start), fit_1.rho1(end)         
    area_1 = int2_1 - int1_1
    norm_factor_1 = 1.0 / area_1                                                  # First fitting function g1
    
    int1_2, int2_2 = fit_2.rho1(start), fit_2.rho1(end)
    area_2 = int2_2 - int1_2
    norm_factor_2 = 1.0 / area_2                                                  # Second fitting function g2

def g(s,fit_1=expo_fit,fit_2=peak_fit):                                           # The two-channel g = alpha1 g1 + alpha2 g2
    ret = alpha_1*norm_factor_1*fit_1.g1(s) + alpha_2*norm_factor_2*fit_2.g1(s)
    return ret

def MonteCarloIntegrate(N,rounds,fit_1=expo_fit,fit_2=peak_fit,start=start,end=end):
    '''
    Monte Carlo integrator by importance sampling
    N: number of points per round
    rounds: number of rounds
    fit_1: the first fitting function
    fit_2: the second fitting function
    inte: function, the integrand
    start, end: floats, range of integration
    '''
    global alpha_1,alpha_2

    initialise(start=start,end=end)

    f = open('multi_peak_opti.csv','w+')
    f.write('No.,Result,Error,W,alpha_1,alpha_2\n')

    ans = 0.0
    sum_y = 0.0
    sum_ysq = 0.0

    #print('norm_factor_1 = ',norm_factor_1)
    #print('norm_factor_2 = ',norm_factor_2)

    for i in range(rounds):
 
        rhos1 = np.random.uniform(int1_1,int2_1,N)                                  # Sampling from g1
        rhos2 = np.random.uniform(int1_2,int2_2,N)                                  # Sampling from g2
        channels = np.random.uniform(0.0,1.0,N)                                     # Sampling of channels
        s = np.zeros(N,dtype=float)
        y = np.zeros(N,dtype=float)
        
        s[channels<alpha_1] = fit_1.s1(rhos1[channels<alpha_1])
        s[channels>alpha_1] = fit_2.s1(rhos2[channels>alpha_1])                     # Compute s values for each channel   
        y[channels<alpha_1] = inte.f(s[channels<alpha_1]) / g(s[channels<alpha_1])
        y[channels>alpha_1] = inte.f(s[channels>alpha_1]) / g(s[channels>alpha_1])  # Compute y = f/g for each channel
        
        sum_y += np.sum(y)
        sum_ysq += np.sum(y*y)
        ans = sum_y / ((i+1)*N)                                                     # Answer = average of y

        err = np.sqrt((sum_ysq / ((i+1)*N) - (sum_y/((i+1)*N)) **2) / ((i+1)*N))    # Error

        #print('W=',np.sum(y*y)/N)

        dw_1 = np.sum(fit_1.g1(s)*y*y*norm_factor_1 / g(s))                         # dW / (d alpha1)
        dw_2 = np.sum(fit_2.g1(s)*y*y*norm_factor_2 / g(s))                         # dW / (d alpha2)
        
        #print('dW/da_1=',dw_1,'dW/da_2',dw_2)
        #print('alpha_1 =',alpha_1,'alpha_2 =',alpha_2)

        alpha1n = alpha_1 * np.sqrt(dw_1)
        alpha2n = alpha_2 * np.sqrt(dw_2)                                           # New alphas
        alpha_1 = alpha1n / (alpha1n+alpha2n)                         
        alpha_2 = alpha2n / (alpha1n+alpha2n)                                       # Normalisation

        print('No. =',(i+1)*N,'Ans =',ans,'Error=',err)
        f.write(str((i+1)*N)+','+str(ans)+','+str(err)+','+str(np.sum(y*y)/N)+','+str(alpha_1)+','+str(alpha_2)+'\n')
        
        #print('#######################\n')

    
    f.close()
    return ans,sum_y

if(__name__=='__main__'):
    mc, sum_y = MonteCarloIntegrate(10000000,100)
    print(mc,sum_y)