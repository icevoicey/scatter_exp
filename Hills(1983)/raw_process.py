#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 18:08:38 2017

@author: icevoicey
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.constants import G
import math



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

names=['succ','hier1','hier2','x0','y0','z0','x1','y1','z1','x2','y2','z2','vx0','vy0','vz0','vx1','vy1','vz1','vx2','vy2','vz2','b/pc', 'dL','dE','res','ncount','ttot']

hier1,hier2,x0,y0,z0,x1,y1,z1,x2,y2,z2,vx0,vy0,vz0,vx1,vy1,vz1,vx2,vy2,vz2,b,dL,dE,res,a230,e230=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
for i in ['1']:
    resulti=pd.read_csv('result'+i, sep=',', names=names, nrows=50)
    for j in range(len(resulti)):
        hier1.append(resulti['hier1'][j])
        hier2.append(resulti['hier2'][j])
        x0.append(float(resulti['x0'][j])/100)
        y0.append(float(resulti['y0'][j])/100)
        z0.append(float(resulti['z0'][j])/100)
    
        x1.append(float(resulti['x1'][j])/100)
        y1.append(float(resulti['y1'][j])/100)
        z1.append(float(resulti['z1'][j])/100)
    
        x2.append(float(resulti['x2'][j])/100)
        y2.append(float(resulti['y2'][j])/100)
        z2.append(float(resulti['z2'][j])/100)
       
        vx0.append(float(resulti['vx0'][j])/100)
        vy0.append(float(resulti['vy0'][j])/100)
        vz0.append(float(resulti['vz0'][j])/100)
    
        vx1.append(float(resulti['vx1'][j])/100)
        vy1.append(float(resulti['vy1'][j])/100)
        vz1.append(float(resulti['vz1'][j])/100)
    
        vx2.append(float(resulti['vx2'][j])/100)
        vy2.append(float(resulti['vy2'][j])/100)
        vz2.append(float(resulti['vz2'][j])/100)
    
        dL.append(float(resulti['dL'][j]))
        dE.append(float(resulti['dE'][j]))
        res.append(resulti['res'][j])
        

def my_mod(x,y,z):
    return (x**2+y**2+z**2)**0.5

def r_t(ma,mb):
    return (ma/mb)**(1/3)*R_M(mb)

Elist,dElist,Clist=[],[],[]
E0=-G*m0*m1/(2*a01)
for i in range(len(hier1)):
    d=my_mod(x1[i]-x0[i],y1[i]-y0[i],z1[i]-z0[i])
    vr=my_mod(vx1[i]-vx0[i],vy1[i]-vy0[i],vz1[i]-vz0[i])
    E=0.5*(m1*m0/(m1+m0))*vr**2-G*m1*m0/d
    Elist.append(E)
    dElist.append(E-E0)
Clist=[(dElist[i]/E0)*(m0+m1)/(2*m2) for i in range(len(Elist))]           
final=pd.DataFrame({'E':Elist,'dE':dElist,'C':Clist})
final.to_csv('final1',sep=',',index=False)
