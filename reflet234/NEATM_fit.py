# -*- coding: utf-8 -*-
from re import DOTALL
import matplotlib as mpl
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import norm
from scipy import integrate
from matplotlib.ticker import FuncFormatter
import NEATM
from NEATM import *
# #loading position files:
wlen = np.loadtxt('ObsWLs.txt')
epoch = np.loadtxt('ObsECs.txt')
#fobs = np.loadtxt('fobs.txt')
w12 = np.loadtxt('w12_2010co1.txt')
w34 = np.loadtxt('w34_2010co1.txt')
w1 = w12[0:32]
w2 = w12[32:64]
w3 = w34[0:32]
w4 = w34[32:64]
fobs = np.zeros((len(w1)*4,3))
for i in range(0,len(w1)):
    fobs[i*4+0,0:3] = [3.4,w1[i],w1[i]*0.1]
    fobs[i*4+1,0:3] = [4.6,w2[i],w2[i]*0.1]
    fobs[i*4+2,0:3] = [11.98,w3[i],w3[i]*0.1]
    fobs[i*4+3,0:3] = [22.0,w4[i],w4[i]*0.1]


def para_fit(yita,Dia,wf,A,epoch,wlen):
    flux = []
    frLamb = []
    frLomm = []
    w_final = []
    for i in range(0,len(wlen)):
        astp = epoch[i,0:3]
        obsp = epoch[i,3:6]
        for j in range(0,len(wlen[0,:])):
            wl = wlen[i,j]
            if wl != 0:
                #print(wl,astp,obsp)
                ff = get_flux_ref(astp,obsp,Dia,wl,yita,A) #flux and reflected sunlight
                flux.append(ff[0])
                frLamb.append(ff[1])
                frLomm.append(ff[2])
                w_final.append(wl)
    flux = np.array(flux)
    frLamb = np.array(frLamb)
    frLomm = np.array(frLomm)
    fref = wf * frLamb + frLomm
    w_final = np.array(w_final)
    return flux,fref,w_final

# for i in range(len(fobs)-1,-1,-1):
#     if fobs[i,0] == 3.4 or fobs[i,0]== 4.6:
#         fobs = np.delete(fobs,i,axis=0)

chiall = []  
Dall = []
yita_all = []
for wf in range(0,11):
    wf = wf * 0.05
    for D in range(280,500,10):
        #variables---------------------------------------------------------------------------
        #slope parameter:
        G = 0.15
        #abs magnitude
        mag = 21.5
        #geometric albedo:
        pv = (1329*10**(-mag/5)/D*1000)**2
        #phase integral
        qph = 0.15 + 0.684 * G
        #bond albedo
        A = pv * qph
        for yita in range(80,450,5):
            yita = yita * 0.01
            fm,fref,wfinal = para_fit(yita,D,wf,A,epoch,wlen)
            chi2 = 0
            for i in range(0,len(fm)):
                chi2 = chi2 +  ((fm[i]+fref[i] - fobs[i,1])/(fobs[i,1]*0.1))**2
            #chi2 = 1/(len(wlen)*2-3) * sum(((fm - fobs[:,1])/(fobs[:,1]*0.1))**2)
            chi2 = chi2/(len(fobs)-3)
            chiall.append(chi2)
            print(chi2,D,yita,wf)
chiall = np.array(chiall)
n = np.where(chiall == chiall.min())

k = 0
chif = []
d1 = []
yita1 = []
for nw in range(0,11):
    for i in range(280,500,10):
        for j in range(80,450,5):
            k = k + 1
            if chiall[k-1] <= chiall.min() + 1:
                d1.append(i+10)
                yita1.append((j+1)*0.05)
            if k == n[0]:
                #print(i,j)
                nnw = nw
                nd = i
                ny = j
                
D_final = nd + 10
yita_final = (ny + 5) * 0.01
nw_final = (nnw+1) * 0.05
pv_final = (1329*10**(-mag/5)/D_final*1000)**2
A_final = pv_final * qph

flux_final,fref,w_final = para_fit(yita_final,D_final,nw_final,A_final,epoch,wlen)
    
    



font1 = {'family': 'Times New Roman', 'weight': 'normal', 'size': 30 }
font2 = {'family': 'Times New Roman', 'weight': 'normal', 'size': 30 }




# wlen = np.arange(5,35,0.2)
# flux_final = np.zeros((len(wlen),1))
# for j in range(0,len(wlen)):
#     flux_final[j] = get_flux_ref(astp,obsp,Dia,wlen[j],yita,A)
# plt.subplots_adjust(left=0.14, right=0.95, top=0.9, bottom=0.16, wspace=0.45)
# plt.plot(wlen,flux_final)
# plt.tick_params(labelsize=30)
# plt.legend(loc='upper right', numpoints=1, fancybox=True, prop = {'family': 'Times New Roman', 'weight': 'normal', 'size': 18 }, ncol = 2)
# plt.xlabel('wavelength (\mu m)',font2)
# plt.ylabel('theoretical flux (mjy)',font2)