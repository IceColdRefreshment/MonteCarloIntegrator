import numpy as np
from matplotlib import pyplot as plt
import integrand_multipeaks as inte
import expo_fit
import peak25 as peak_fit

start = 4.0
end = 40.0
zero_0 = 4.0
peak_0 = 2.0+2.0*np.sqrt(3.0)

alpha_1 = 0.98
alpha_2 = 0.02
#alpha_1 = 0.05
#alpha_2 = 0.95

int1_1, int2_1, int3_1, int4_1 = expo_fit.pho1(start), expo_fit.pho1(peak_0), expo_fit.pho2(peak_0), expo_fit.pho2(end)
area1_1 = int2_1 - int1_1
area2_1 = int4_1 - int3_1
area_1 = area1_1 + area2_1
norm_factor_1 = 1.0 / area_1

int1_2, int2_2 = peak_fit.pho1(start), peak_fit.pho1(end)
area_2 = int2_2 - int1_2
norm_factor_2 = 1.0 / area_2


def g(s):
    ret = np.zeros(len(s),dtype=float)
    ret[s<peak_0] = alpha_1*norm_factor_1*expo_fit.g1(s[s<peak_0]) + alpha_2*norm_factor_2*peak_fit.g1(s[s<peak_0])
    ret[s>peak_0] = alpha_1*norm_factor_1*expo_fit.g2(s[s>peak_0]) + alpha_2*norm_factor_2*peak_fit.g1(s[s>peak_0])
    return ret

def directIntegrate(N):
    xs = np.random.uniform(start,end,N)
    ys = inte.f(xs)
    ans = np.sum(ys)*(end-start)/N
    return ans

def MonteCarloIntegrate(N,rounds):
    global alpha_1,alpha_2

    f = open('multi_peak_opti.csv','w+')
    f.write('No.,Result,Error,W,alpha_1,alpha_2\n')


    ans = 0.0
    cnt_1 = 0
    cnt_2 = 0
    sum_y = 0.0
    sum_ysq = 0.0

    print('norm_factor_1 = ',norm_factor_1)
    print('norm_factor_2 = ',norm_factor_2)

    for i in range(rounds):

        #alpha_1 = 0.9+i * 0.01
        #alpha_2 = 1-alpha_1

        rhos_1 = np.random.uniform(0.0,area_1,N)
        rhos_2 = np.random.uniform(0.0,area_2,N)
        channels = np.random.uniform(0.0,1.0,N)

        regime_1 = np.zeros(N,dtype=int)
        regime_1[rhos_1<area1_1] = 1
        regime_1[rhos_1>area1_1] = 2

        rhoo1_1 = rhos_1 + int1_1
        s1_1 = expo_fit.s1(rhoo1_1)
        y1_1 = inte.f(s1_1) / g(s1_1)
        rhoo2_1 = rhos_1 - area1_1 + int3_1
        s2_1 = expo_fit.s2(rhoo2_1)
        y2_1 = inte.f(s2_1) / g(s2_1)

        rhoo_2 = rhos_2 + int1_2
        s_2 = peak_fit.s1(rhoo_2)
        y_2 = inte.f(s_2) / g(s_2)


        s = np.zeros(N,dtype=float)
        y = np.zeros(N,dtype=float)
        rc = np.zeros(N,dtype=int)

        rc[(regime_1==1)*(channels<alpha_1)==1]=1
        rc[(regime_1==2)*(channels<alpha_1)==1]=2
        rc[channels>alpha_1]=3

        s[rc==1]=s1_1[rc==1]
        s[rc==2]=s2_1[rc==2]
        s[rc==3]=s_2[rc==3]
        
        y[rc==1]=y1_1[rc==1]
        y[rc==2]=y2_1[rc==2]
        y[rc==3]=y_2[rc==3]

        #cnt_1 += np.sum(regime[regime==1])
        #cnt_2 += np.sum(regime[regime==2])/2

        sum_y += np.sum(y)
        sum_ysq += np.sum(y*y)
        ans = sum_y / ((i+1)*N)

        err = np.sqrt((sum_ysq / ((i+1)*N) - (sum_y/((i+1)*N)) **2) / ((i+1)*N))

        print('W=',np.sum(y*y)/N)

        dw_2=0.0
        #dw_2+=np.sum(peak_fit.g1(s[rc==3])*y[rc==3]*y[rc==3]*norm_factor_2)
        dw_2 += np.sum(peak_fit.g1(s)*y*y*norm_factor_2 / g(s))
        dw_1=0.0
        dw_1 += np.sum(expo_fit.g1(s[s<peak_0])*norm_factor_1*y[s<peak_0]*y[s<peak_0] / g(s[s<peak_0]))
        dw_1 += np.sum(expo_fit.g2(s[s>peak_0])*norm_factor_1*y[s>peak_0]*y[s>peak_0] / g(s[s>peak_0]))
        #dw_1+=np.sum(expo_fit.g1(s[rc==1])*y[rc==1]*y[rc==1]*norm_factor_1)
        #dw_1+=np.sum(expo_fit.g2(s[rc==2])*y[rc==2]*y[rc==2]*norm_factor_1)

        #alpha_1 -= 0.2*(dw_1 - dw_2)/(dw_1 + dw_2)
        #alpha_2 += 0.2*(dw_1 - dw_2)/(dw_1 + dw_2)

        
        print('dW/da_1=',dw_1,'dW/da_2',dw_2)
        print('alphas',alpha_1,alpha_2)

        alpha1n = alpha_1 * np.sqrt(dw_1)
        alpha2n = alpha_2 * np.sqrt(dw_2)
        alpha_1 = alpha1n / (alpha1n+alpha2n)
        alpha_2 = alpha2n / (alpha1n+alpha2n)

        print((i+1)*N,sum_y,ans,err)
        f.write(str((i+1)*N)+','+str(ans)+','+str(err)+','+str(np.sum(y*y)/N)+','+str(alpha_1)+','+str(alpha_2)+'\n')
        
        print('#######################\n')
    
    f.close()
    return ans,sum_y

if(__name__=='__main__'):
    mc, sum_y = MonteCarloIntegrate(10000000,100)
    print(mc,sum_y)