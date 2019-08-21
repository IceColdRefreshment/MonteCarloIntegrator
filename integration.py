import numpy as np
from matplotlib import pyplot as plt
import integrand as inte
import expo_fit_new as fit

start = 4.0
end = 25.0

def MonteCarloIntegrate(N,rounds,fit=fit,inte=inte,start=start,end=end):
    """ 
        Monte Carlo Integration by Importance Sampling
        N : number of intergration per round
        rounds: number of rounds. Note: total points of integration = N * rounds
        inte: function of integrand
        fit: function used to fit
        start, end: range of integration
    """
    f = open('combined.csv','w+')
    f.write('No.,Result,Error\n')                                      # File operations

    ans = 0.0                                                          # Result of integration
    sum_y = 0.0                                                        # Sum of y = f/g
    sum_ysq = 0.0                                                      # Sum of y^2 = (f/g)^2

    area = fit.rho1(end) - fit.rho1(start)
    norm_factor = 1.0 / area                                           # Normalization factor

    print('Normalization Factor =',norm_factor)

    for i in range(rounds):
        rhos = np.random.uniform(fit.rho1(start),fit.rho1(end),N)      # Uniform sampling of rho(s)
        s = fit.s1(rhos)                                               # Transfer rho into s by s = rho^-1  
        y = inte.f(s) / (norm_factor*fit.g1(s))                        # y = f/g

        sum_y += np.sum(y)
        sum_ysq += np.sum(y*y)
        ans = sum_y / ((i+1)*N)                                        # Answer = average of y

        err = np.sqrt((sum_ysq / ((i+1)*N) - (sum_y/((i+1)*N)) **2) / ((i+1)*N))    # Compute the error
     
        print("No. =",(i+1)*N,"Ans =",ans,"Error =",err)
        f.write(str((i+1)*N)+','+str(ans)+','+str(err)+'\n')
    f.close()
    return ans,sum_y

if(__name__=='__main__'):
    mc, sum_y = MonteCarloIntegrate(10000000,100)
    print(mc,sum_y)