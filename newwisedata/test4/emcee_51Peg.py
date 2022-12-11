#Import required libraries
#coding=utf-8
import sys
import numpy as np
import batman #package by Laura Kreidberg:  http://astro.uchicago.edu/~kreidberg/batman/
import emcee  #package by Dan Foreman-Mackey:  http://dan.iel.fm/emcee/current/
import corner #package by Dan Foreman-Mackey:  http://corner.readthedocs.io/en/latest/
import matplotlib
import matplotlib.pyplot as plt
from numpy import sin,cos,tan,arcsin,arccos,arctan,sqrt,pi,radians,degrees,arctan2,median,mean,array,log
from scipy import constants
#%matplotlib inline

#datafile=np.loadtxt(sys.argv[1],unpack=True)
datafile=np.loadtxt("ELODIE_RV_Naef_2004.dat",unpack=True)
rv_time=datafile[0]-2450000
d_rv=datafile[1]
e_rv=datafile[2]


#parafile=np.loadtxt(sys.argv[2],skiprows=1)
parafile=np.loadtxt("para_51Peg_ELODIE_RV.dat",skiprows=1)
K_guess=parafile[0]
sesino_guess=parafile[1]
secoso_guess=parafile[2]
M0_guess=parafile[3]
per_guess=parafile[4]
gamma_guess=parafile[5]

#系统及恒星参数
au=constants.au
sec=86400.
ratio=0.942002e-3
dis=14.7
ms=1.11
MS=1.98855e30
MJ=1.8986e27
yr=365.25


#定义函数，用于求解开普勒方程
def newton(E,e,M):
    eps=1
    while(abs(eps)>1e-8):
        E1=E-(E-e*sin(E)-M)/(1-e*cos(E))
        eps=E1-E
        E=E1
    return E


def func_rv(time,theta):
    K, sesino, secoso, M0, per, gamma= theta
    M0=radians(M0)
    omega=arctan2(sesino,secoso)
    ecc=sesino**2+secoso**2
    E_rv=list(map(lambda x:newton(x,ecc,x),(2*pi)/per*time-M0))
    sinf=(sqrt(1-ecc**2)*sin(E_rv))/(1-ecc*cos(E_rv))
    cosf=(cos(E_rv)-ecc)/(1-ecc*cos(E_rv))
    v0=K*ecc*cos(omega)+gamma
    rv=K*cos(omega)*cosf-K*sin(omega)*sinf+v0
    return rv



#prior
def lnprior(theta):
    K, sesino, secoso, M0, per, gamma= theta
    K_min = max(K_guess[0]-3*max(K_guess[1],K_guess[2]), 0.01)
    K_max= K_guess[0]+3*max(K_guess[1],K_guess[2])
    sesino_min= max(sesino_guess[0]-3*max(sesino_guess[1],sesino_guess[2]),-0.2)
    sesino_max= min(sesino_guess[0]+3*max(sesino_guess[1],sesino_guess[2]),0.2)
    secoso_min= max(secoso_guess[0]-3*max(secoso_guess[1],secoso_guess[2]),-0.2)
    secoso_max= min(secoso_guess[0]+3*max(secoso_guess[1],secoso_guess[2]),0.2)
    M0_min= max(M0_guess[0]-3*max(M0_guess[1],M0_guess[2]),0.)
    M0_max= min(M0_guess[0]+3*max(M0_guess[1],M0_guess[2]), 360.)
    per_min= max(per_guess[0]-3*max(per_guess[1],per_guess[2]),0.2)
    per_max= per_guess[0]+3*max(per_guess[1],per_guess[2])
    gamma_min= gamma_guess[0]-3*max(gamma_guess[1],gamma_guess[2])
    gamma_max= gamma_guess[0]+3*max(gamma_guess[1],gamma_guess[2])
    if K_min<K<K_max and per_min<per<per_max and gamma_min<gamma<gamma_max and\
            sesino_min<sesino<sesino_max and secoso_min<secoso<secoso_max and M0_min<M0<M0_max:
        return -(log(per*log(per_max/per_min))+log(K*log(K_max/K_min)))
        #return 0.0
    return -np.inf


#likelihood function
def lnlike(theta, time_rv, d_rv, e_rv):
    K_b, sesino_b, secoso_b, M0_b, per_b, gamma= theta
    e_b=sesino_b**2+secoso_b**2
    o_b=arctan2(sesino_b,secoso_b)
    rv=func_rv(time_rv,[K_b,sesino_b,secoso_b,M0_b,per_b,gamma])

    residuals_rv = d_rv-rv
    ln_likelihood = -0.5*(np.sum((residuals_rv/e_rv)**2 + np.log(2.0*np.pi*(e_rv)**2)))
    return ln_likelihood


#posterior probability
def lnprob(theta, time_rv, d_rv, e_rv):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, time_rv, d_rv, e_rv)

#Now for some initial parameter guesses
g_K = K_guess[0]
g_sesino= sesino_guess[0]
g_secoso= secoso_guess[0]
g_M0 = M0_guess[0]
g_per = per_guess[0]
g_gamma = gamma_guess[0]


theta = [g_K, g_sesino, g_secoso, g_M0, g_per, g_gamma]


#initialize sampler
ndim, nwalkers = len(theta), 200
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (rv_time, d_rv, e_rv))
pos = [theta + 1e-6*np.random.randn(ndim) for i in range(nwalkers)]


#run mcmc
sampler.run_mcmc(pos,5000);


samples = sampler.chain[:, 500:, :].reshape((-1, ndim))
K_mcmc, sesino_mcmc, secoso_mcmc, M0_mcmc, per_mcmc, gamma_mcmc= list(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0))))


#f = open("chain.dat", "w")
#f.close()
#
#for result in sampler.sample(pos, iterations=500, storechain=False):
#   position = result[0]
#   f = open("chain.dat", "a")
#   for k in range(position.shape[0]):
#       f.write("%15.10f\t%15.10f\t%15.10f\t%15.10f\t%15.10f\t%15.10f\t%15.10f\t%15.10f\t%15.10f\n"%tuple(position[k]))
#   f.close()


#print u3_mcmc
#print u4_mcmc

import corner
fig = corner.corner(samples,bins=40,labels=[r'$\alpha$','$\sqrt{e}\sin\omega$','$\sqrt{e}\cos\omega$','$M0$','$P$','$\gamma$'],\
                    truths=[K_mcmc[0], sesino_mcmc[0], secoso_mcmc[0], M0_mcmc[0], per_mcmc[0],gamma_mcmc[0]],\
                    smooth = 1.0,show_titles = True,title_fmt = '.6f',quantiles = [0.1586,0.5,0.8414],verbose = True)
fig.savefig("triangle_emcee_51Peg_ELODIE_RV.jpg")
#plt.show()

out_para=open('para_51Peg_ELODIE_RV.dat','w')
out_para.write("%15s\t%15s\t%15s\n"%("Median","Upper_sigma","Lower_sigma"))
out_para.write("%15.10f\t%15.10f\t%15.10f\n"%K_mcmc)
out_para.write("%15.10f\t%15.10f\t%15.10f\n"%sesino_mcmc)
out_para.write("%15.10f\t%15.10f\t%15.10f\n"%secoso_mcmc)
out_para.write("%15.10f\t%15.10f\t%15.10f\n"%M0_mcmc)
out_para.write("%15.10f\t%15.10f\t%15.10f\n"%per_mcmc)
out_para.write("%15.8f\t%15.8f\t%15.8f\n"%gamma_mcmc)
out_para.close()


time_con=np.linspace(rv_time[0],rv_time[-1],2000)
rv_cal_dis=func_rv(rv_time,[K_mcmc[0], sesino_mcmc[0], secoso_mcmc[0], M0_mcmc[0], per_mcmc[0],gamma_mcmc[0]])
#rv_cal_dis=func_rv(rv_time,[K_b_t, sesino_b_t, secoso_b_t, M0_b_t, per_b_t])
rv_cal_con=func_rv(time_con,[K_mcmc[0], sesino_mcmc[0], secoso_mcmc[0], M0_mcmc[0], per_mcmc[0],gamma_mcmc[0]])
#rv_cal_con=func_rv(time_con,[K_b_t, sesino_b_t, secoso_b_t, M0_b_t, per_b_t])

rv_res=d_rv-rv_cal_dis
chi_sqaure = np.sqrt(sum(rv_res**2/(e_rv)**2)/(len(rv_time)-len(theta)))
print(chi_sqaure)

plt.figure()
plt.errorbar(rv_time,d_rv,yerr=e_rv,fmt='ro')
plt.plot(time_con,rv_cal_con,lw=0.2)
plt.savefig('51Peg_ELODIE_RV.jpg')

