# -*- coding: utf-8 -*-
import matplotlib as mpl
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import norm
from scipy import integrate
from matplotlib.ticker import FuncFormatter



def planck(T,wl): #wl in unit of um
    c = 3.0 * 1e8 
    h = 6.63 * 1e-34
    k = 1.38 * 1e-23
    wlen = wl * 1e-6
    jy = 3e-12 / wl**2
    mjy = jy * 1e-3
    if T==0:
        B=0
    else:
        B= 1e-6 * 2 * h * c**2  / (wlen**5*(math.exp(h*c/(wlen*k*T))-1)) 
    pl = B / mjy
    #print(B)
    return pl

def get_angle(a):
    pi = 3.141592653
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


def get_flux_noref(astp,obsp,Dia,wlenth,yita,A):
    '''
    flux without ref
    '''
    #constants---------------------------------------------------------------------------
    #emissivity:
    epsi = 0.9
    #stefan-boltzman:
    sigmas = 5.67*10**(-8)
    #solar constant:
    f_solar = 1367.5
    #planck constant:
    h = 6.62607015 * 10 ** (-34)
    #speed of light
    c = 3.0 * 10 ** 8 
    #boltzmann constant
    kb = 1.380649 * 10 ** (-23)
    #astronomical unit in unit of meter
    au = 1.496 * 10 ** 11
    #solar radius
    sr = 6.9550826e8
    
    ph_dis = get_phase_dist(astp,obsp)
    alpha = ph_dis[0] * pi / 180 # phase angle in arc
    dast = ph_dis[1]
    dao = ph_dis[3]
    #print(dast,dao,alpha)
    Nd = int(10) #divide theta and phi into Nd 
    T_ss = ((1 - A) * f_solar / epsi / yita / sigmas / dast ** 2) ** 0.25
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
    #fref = 0.0
    #ref_lomm = 0.0
    flux_con = epsi * Dia ** 2 * pi  * h * c ** 2  / (wlenth ** 5) / (2 * (dao * au) ** 2)
    #ref_con = A * Dia**2 * pi * h * c**2 / (wlenth**5) / (2 * dao * au)**2 * (sr / (dast * au))**2 
    T_solar = 5778 # solar surface temperature
    for j in range(0,len(phi)):
        mu_inci = math.cos(phi[j]) #incident angle of latitude phi
        if alpha < 0:
            mu_emer = np.abs(np.abs(alpha)-phi[j])#emergence angle of latitude phi
        else:
            mu_emer = np.abs(alpha) + np.abs(phi[j])
        for k in range(int(nj),len(theta)):
            temp[j,k] = T_ss * cos(theta[k]) ** 0.25 * cos(phi[j]) ** 0.25
            flux = flux + flux_con * math.cos(phi[j])**2 * np.abs(math.cos(alpha - theta[k]))  / (exp(h * c / (wlenth * kb * temp[j,k])) - 1) * dphi * dtheta * wlenth ** 2 / c * 10 ** 29 # obtain flux in unit of mjy
            #fref = fref +  ref_con * math.cos(phi[j])**2 * np.abs(math.cos(alpha - theta[k]))  / (exp(h * c / (wlenth * kb * T_solar)) - 1) * dphi * dtheta * wlenth ** 2 / c * 10 ** 29 #Lambertian reflection
            #ref_lomm = ref_lomm + ref_con * math.cos(phi[j])**2 * np.abs(math.cos(alpha - theta[k]))  / (math.exp(h * c / (wlenth * kb * T_solar)) - 1) * dphi * dtheta * wlenth ** 2 / c * 10 ** 29 * 1/(math.cos(mu_inci)+math.cos(mu_emer)) #Lommel-Seeliger reflection
    #frLamb = fref
    #frLomm = ref_lomm
    #print(shape(frLomm))
    return flux#,frLamb,frLomm
def get_flux_ref(astp,obsp,Dia,wlenth,yita,A,Hv):
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
    
    #constants---------------------------------------------------------------------------
    #emissivity:
    epsi = 0.9
    #stefan-boltzman:
    sigmas = 5.67*10**(-8)
    #solar constant:
    f_solar = 1367.5
    #planck constant:
    h = 6.62607015 * 10 ** (-34)
    #speed of light
    c = 3.0 * 10 ** 8 
    #boltzmann constant
    kb = 1.380649 * 10 ** (-23)
    #astronomical unit in unit of meter
    au = 1.496 * 10 ** 11
    #solar radius
    sr = 6.9550826e8
    
    ph_dis = get_phase_dist(astp,obsp)
    alpha = ph_dis[0] * pi / 180 # phase angle in arc
    dast = ph_dis[1]
    dao = ph_dis[3]
    #print(dast,dao,alpha)
    Nd = int(10) #divide theta and phi into Nd 
    T_ss = ((1 - A) * f_solar / epsi / yita / sigmas / dast ** 2) ** 0.25
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
    flux_con = epsi * Dia ** 2 * pi  * h * c ** 2  / (wlenth ** 5) / (2 * (dao * au) ** 2)
    pv = (1329*pow(10,-Hv/5)/(Dia*0.001))**2
    ref_con = pv * Dia**2 * pi * h * c**2 / (wlenth**5) / (2 * dao * au)**2 * (sr / (dast * au))**2 
    T_solar = 5778 # solar surface temperature
    for j in range(0,len(phi)):
        mu_inci = math.cos(phi[j]) #incident angle of latitude phi
        if alpha < 0:
            mu_emer = np.abs(np.abs(alpha)-phi[j])#emergence angle of latitude phi
        else:
            mu_emer = np.abs(alpha) + np.abs(phi[j])
        for k in range(int(nj),len(theta)):
            temp[j,k] = T_ss * cos(theta[k]) ** 0.25 * cos(phi[j]) ** 0.25
            flux = flux + flux_con * math.cos(phi[j])**2 * np.abs(math.cos(alpha - theta[k]))  / (exp(h * c / (wlenth * kb * temp[j,k])) - 1) * dphi * dtheta * wlenth ** 2 / c * 10 ** 29 # obtain flux in unit of mjy
            fref = fref +  ref_con * math.cos(phi[j])**2 * np.abs(math.cos(alpha - theta[k]))  / (exp(h * c / (wlenth * kb * T_solar)) - 1) * dphi * dtheta * wlenth ** 2 / c * 10 ** 29 #Lambertian reflection
            ref_lomm = ref_lomm + ref_con * math.cos(phi[j])**2 * np.abs(math.cos(alpha - theta[k]))  / (math.exp(h * c / (wlenth * kb * T_solar)) - 1) * dphi * dtheta * wlenth ** 2 / c * 10 ** 29 * 1/(math.cos(mu_inci)+math.cos(mu_emer)) #Lommel-Seeliger reflection
    frLamb = fref
    frLomm = ref_lomm
    #print(shape(frLomm))
    return flux,frLamb,frLomm

