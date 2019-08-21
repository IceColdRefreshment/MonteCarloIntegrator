import numpy as np
from matplotlib import pyplot as plt
import integrand as inte
import combine_fit as fit
#import cauchy_fit as fit

start = 4.0
end = 25.0
#interval = 5.0
zero_0 = 4.0
peak_0 = 2.0+2.0*np.sqrt(3.0)

def directIntegrate(N):
    xs = np.random.uniform(start,end,N)
    ys = inte.f(xs)
    ans = np.sum(ys)*(end-start)/N
    return ans

def MonteCarloIntegrate(N,rounds):

    f = open('combined.csv','w+')
    f.write('No.,Result,Error\n')

    ans = 0.0
    cnt_1 = 0
    cnt_2 = 0
    sum_y = 0.0
    sum_ysq = 0.0

    int1, int2, int3, int4 = fit.pho1(start), fit.pho1(peak_0), fit.pho2(peak_0), fit.pho2(end)
    area1 = int2 - int1
    area2 = int4 - int3
    area = area1 + area2
    norm_factor = 1.0 / area

    print('norm_factor = ',norm_factor)

    for i in range(rounds):

        rhos = np.random.uniform(0.0,area,N)

        regime = np.zeros(N,dtype=int)
        regime[rhos<area1] = 1
        regime[rhos>area1] = 2

        rhoo1 = rhos + int1
        s1 = fit.s1(rhoo1)
        y1 = inte.f(s1) / (norm_factor*fit.g1(s1))
        rhoo2 = rhos - area1 + int3
        s2 = fit.s2(rhoo2)
        y2 = inte.f(s2) / (norm_factor*fit.g2(s2))

        s = np.zeros(N,dtype=float)
        y = np.zeros(N,dtype=float)
        s[regime==1]=s1[regime==1]
        s[regime==2]=s2[regime==2]
        y[regime==1]=y1[regime==1]
        y[regime==2]=y2[regime==2]

        cnt_1 += np.sum(regime[regime==1])
        cnt_2 += np.sum(regime[regime==2])/2

        sum_y += np.sum(y)
        sum_ysq += np.sum(y*y)
        ans = sum_y / ((i+1)*N)
        #err = ans / np.sqrt((i+1)*N)

        err = np.sqrt((sum_ysq / ((i+1)*N) - (sum_y/((i+1)*N)) **2) / ((i+1)*N))

        print((i+1)*N,sum_y,ans,err)
        f.write(str((i+1)*N)+','+str(ans)+','+str(err)+'\n')
    print('Cnt_1 = ',cnt_1,'Cnt_2=',cnt_2)
    f.close()
    return ans,sum_y

if(__name__=='__main__'):
    mc, sum_y = MonteCarloIntegrate(10000000,100)
    print(mc,sum_y)