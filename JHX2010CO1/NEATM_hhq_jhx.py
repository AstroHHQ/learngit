# -*- coding: utf-8 -*-
import matplotlib as mpl
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import norm
from scipy import integrate
from matplotlib.ticker import FuncFormatter
#constant
#solar constant:
Fsun = 1367.5           # sun constant
sigma = 5.67E-8  #stefan-boltzmann constant
#planck constant:
h = 6.626007015E-34     #plank constant
#q = 0.29+0.684*0.15    #phase integral = 0.29+0.684*G(=0.15)
#emissivity:
epsi = 0.9             #radiance epsilon
#boltzmann constant
kB = 1.38064852E-23     #boltzmann constant  j/k
cl = 299792458.0        #lightspeed m/s
au = 149597870700.0
pi = 3.1415926535
Rsun = 0.00465*au   #Rsun *m
Tsun = 5778         #Tsun  K
#solar radius
sr = 6.9550826e8

def planck(T,wl): #wl in unit of um
    wlen = wl * 1e-6
    jy = 3e-12 / wl**2
    mjy = jy * 1e-3
    if T==0:
        B=0
    else:
        B= 1e-6 * 2 * h * cl**2  / (wlen**5*(math.exp(h*cl/(wlen*kB*T))-1)) 
    pl = B / mjy
    #print(B)
    return pl

def get_angle(a):
    x = a[0]
    y = a[1]
    z = a[2]
    c_theta = x / np.sqrt(x**2 + y**2)
    if y > 0:
        theta = math.acos(c_theta)
    else:
        theta = -math.acos(c_theta) + pi
    return theta

def get_phase_dist(past,pobs):
    angle1 = get_angle(past)
    angle2 = get_angle(pobs)
    r1 = -1 * past
    r2 = pobs
    r3 = pobs - past # vector from ast to obs
    d1 = np.sqrt(r1[0]**2 + r1[1]**2 + r1[2]**2)
    d2 = np.sqrt(r2[0]**2 + r2[1]**2 + r2[2]**2)
    d3 = np.sqrt(r3[0]**2 + r3[1]**2 + r3[2]**2)
    cphase = dot(r1,r3)/(d1 * d3)
    if angle1<angle2:
        phase = -np.abs(math.acos(cphase)) * 180 / pi
    else:
        phase = np.abs(math.acos(cphase)) * 180 / pi
    dast = d1 #distance from ast to the Sun 
    dobs = d2 #distance from obs to the Sun
    dao = d3  #distance from obs to ast
    return phase,dast,dobs,dao


#def get_flux_noref(astp,obsp,Dia,wlenth,yita,A):

def get_flux_ref(astp,obsp,Dia,wlenth,yita,A):
    '''
    astp: position vector of ast
    obsp: position vector of obs
    dast: distance from asteroid to Sun
    dobs: distance from observer to Sun
    dao:  distance from asteroid to observer
    alpha: solar phase
    dia: diameter of asteroid
    wlenth: wavelength
    yita: beaming parameter
    A: bond albedo
    '''
    ph_dis = get_phase_dist(astp,obsp)
    alpha = ph_dis[0] * pi / 180 # phase angle in arc
    dast = ph_dis[1]
    dao = ph_dis[3]
    #print(dast,dao,alpha)
    Nd = int(10) #divide theta and phi into Nd 
    T_ss = ((1 - A) * Fsun / epsi / yita / sigma / dast ** 2) ** 0.25
    #print(T_ss)
    phi = np.zeros((Nd,1))
    theta = np.zeros((Nd,1))
    for i in range(0,Nd):#phi , theta is the angle from subsolar point
            phi[i] = -pi/2.0 + i * pi / Nd
            theta[i] = -pi/2.0 + i * pi / Nd
    nj = floor(((np.abs(alpha) - pi / 2.0) + pi / 2.0) / (pi / Nd)) 
    dphi, dtheta = pi/Nd,pi/Nd
    #print (nj)
	
    wlenth = wlenth * 10 ** (-6)
    temp = np.zeros((len(phi),len(theta)))
    flux = 0.0
    fref = 0.0
    ref_lomm = 0.0
    flux_con = epsi * Dia ** 2 * pi  * h * cl ** 2  / (wlenth ** 5) / (2 * (dao * au) ** 2)
    ref_con = A * Dia**2 * pi * h * cl**2 / (wlenth**5) / (2 * dao * au)**2 * (sr / (dast * au))**2 
    T_solar = 5778 # solar surface temperature
    for j in range(0,len(phi)):
        mu_inci = math.cos(phi[j]) #incident angle of latitude phi
        if alpha < 0:
            mu_emer = np.abs(np.abs(alpha)-phi[j])#emergence angle of latitude phi
        else:
            mu_emer = np.abs(alpha) + np.abs(phi[j])
        for k in range(int(nj),len(theta)):
            temp[j,k] = T_ss * cos(theta[k]) ** 0.25 * cos(phi[j]) ** 0.25
            flux = flux + flux_con * math.cos(phi[j])**2 * np.abs(math.cos(alpha - theta[k]))  / (exp(h * cl / (wlenth * kB * temp[j,k])) - 1) * dphi * dtheta * wlenth ** 2 / cl * 10 ** 29 # obtain flux in unit of mjy
            fref = fref +  ref_con * math.cos(phi[j])**2 * np.abs(math.cos(alpha - theta[k]))  / (exp(h * cl / (wlenth * kB * T_solar)) - 1) * dphi * dtheta * wlenth ** 2 / cl * 10 ** 29 #Lambertian reflection
            ref_lomm = ref_lomm + ref_con * math.cos(phi[j])**2 * np.abs(math.cos(alpha - theta[k]))  / (math.exp(h * cl / (wlenth * kB * T_solar)) - 1) * dphi * dtheta * wlenth ** 2 / cl * 10 ** 29 * 1/(math.cos(mu_inci)+math.cos(mu_emer)) #Lommel-Seeliger reflection
    frLamb = fref
    frLomm = ref_lomm
    #print(shape(frLomm))
    return flux,frLamb,frLomm

