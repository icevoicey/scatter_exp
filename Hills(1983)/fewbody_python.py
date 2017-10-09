#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 21:36:54 2017

@author: icevoicey

python programme to calculate related values used in fewbody, such as units.
"""

import numpy as np
from astropy import units as u
from astropy.constants import G 


G=G.value
Msun=(1*u.Msun).to(u.kg).value
year=(1*u.year).to(u.s).value
pc=(1*u.pc).to(u.m).value
AU=(1*u.AU).to(u.m).value
Rsun=(1*u.Rsun).to(u.m).value

m0=1e2*Msun
m1=1e2*Msun
m2=1*Msun

def sigmaf(m): #velocity dispersion from M-sigma relation in Tremaine(2002)
    return (m/(10**8.13*Msun))**(1/4.02)*2e5

def r_influf(m): #influence radius of single-BH. For BHB, m=(m1+m2)/2 
    return G*m/(sigmaf(m)**2)

def r_hf(ma,mb): #hard radius of BHB, ma is the main BH
    mu=ma*mb/(ma+mb)
    return G*mu/(4*sigmaf(ma)**2)

def r_bf(ma,mb): #bound radius of BHB, ma is the main BH
    return G*(ma+mb)/sigmaf(ma)**2


a01=r_bf(m0,m1)
e01=0

v_c=(G*(m0+m1)/a01)**0.5
vinf3=0.1*v_c

def v_unitf(m0=m0,m1=m1,m2=m2,a01=a01):
    m01=m0+m1
    mu=m01*m2/(m01+m2)
    return (G/mu*(m0*m1/a01))**0.5

def l_unitf(a1=a01):
    return a1

def t_unitf(l_unit=l_unitf(),v_unit=v_unitf()):
    return l_unit/v_unit

def r_tidf(tidal_tol=1e-4,m0=m0,m1=m1,m2=m2,a01=a01,e01=e01):
    m01=m0+m1
    r01=(2*(m01+m2)*m2/(m01*m2*tidal_tol))**(1/3)*a01*(1+e01)
    return r01

#def calt(r,a):
#    M=1e7*Msun
#    r_b=1.9e3*pc
#    A=3.6e6*Msun/(0.16**(-1.8)*np.exp(-(0.16/1900)**2)*pc**1.2)
#    rou=A*r**(-1.8)*np.exp(-(r/r_b)**2)
#    t=0.1*G**(-0.5)*M**(0.5)/(r**0.5*rou*a*np.log((M/Msun)*(a/r)))
#    return t/(1.3e10*year)
#
def a_orbitf(vinf,M): #semi major axis of conic orbit
    return -G*M/vinf**2

def e_orbitf(b,vinf,M): #ecentricity of conic orbit
    return (1+(b**2*vinf**4/(G*M)**2))**0.5
#
def r_pf(b,vinf,M):  #pericenter distance of conic orbit
    a=a_orbitf(vinf,M)
    e=e_orbitf(b,vinf,M)
    return a*(1-e) 

def T_bf(a,m1,m2): #period of binary system
    return 2*np.pi*(a**3/(G*(m1+m2)))**0.5    

v_omega=2*np.pi/T_bf(a01,m0,m1)
v0=v_omega*a01/(1+1)
v1=v_omega*a01/(1+1)  

E=0.5*m0*v0**2+0.5*m1*v1**2
Ep=G*m0*m1/a01 

v_unit=v_unitf()
l_unit=l_unitf()
t_unit=t_unitf()
m_unit=l_unit*v_unit**2/G
e_unit=v_unit**2*m_unit
#
#t_stop=130e8*year/t_unit #hubble time as stoptime in unit of t_unit
#
#def rp(b,alpha=0.1):
#    return ((1+b**2*alpha**4)**0.5-1)/alpha**2

