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

m0=1e7*Msun
m1=1e7*Msun
m2=1*Msun
m3=1*Msun


names=['succ','hier1','hier2','x0','y0','z0','x1','y1','z1','x2','y2','z2','x3','y3','z3','vx0','vy0','vz0','vx1','vy1','vz1','vx2','vy2','vz2','vx3','vy3','vz3','a/AU','e','b/pc', 'dL','dE','res']

hier1,hier2,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,vx0,vy0,vz0,vx1,vy1,vz1,vx2,vy2,vz2,vx3,vy3,vz3,a23,e23,b,dL,dE,res,a230,e230=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
#,'2','3','4','5','6','7','8' skipfooter=6+9000000
for i in ['8']:
    resulti=pd.read_csv('result'+i, sep=',', names=names, nrows=10000000)
    for j in range(len(resulti['succ'])):      
        if resulti['succ'][j] == 1 and (resulti['hier2'][j] == 'single-triple' \
                  or resulti['hier2'][j] == 'triple-single'):
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

            x3.append(float(resulti['x3'][j])/100)
            y3.append(float(resulti['y3'][j])/100)
            z3.append(float(resulti['z3'][j])/100)
       
            vx0.append(float(resulti['vx0'][j])/100)
            vy0.append(float(resulti['vy0'][j])/100)
            vz0.append(float(resulti['vz0'][j])/100)

            vx1.append(float(resulti['vx1'][j])/100)
            vy1.append(float(resulti['vy1'][j])/100)
            vz1.append(float(resulti['vz1'][j])/100)

            vx2.append(float(resulti['vx2'][j])/100)
            vy2.append(float(resulti['vy2'][j])/100)
            vz2.append(float(resulti['vz2'][j])/100)

            vx3.append(float(resulti['vx3'][j])/100)
            vy3.append(float(resulti['vy3'][j])/100)
            vz3.append(float(resulti['vz3'][j])/100)            
            
            a23.append(float(resulti['a/AU'][j])*AU)
            e23.append(float(resulti['e'][j]))
            dL.append(float(resulti['dL'][j]))
            dE.append(float(resulti['dE'][j]))
            res.append(resulti['res'][j])

def my_cross(x1,y1,z1,x2,y2,z2):
    x3=y1*z2-y2*z1
    y3=x2*z1-x1*z2
    z3=x1*y2-x2*y1
    return my_mod(x3,y3,z3)

def my_dot(x1,y1,z1,x2,y2,z2):
    x3=x1*x2
    y3=y1*y2
    z3=z1*z2
    return (x3,y3,z3)

def my_mod(x,y,z):
    return (x**2+y**2+z**2)**0.5


def r_t(ma,mb):
    return (ma/mb)**(1/3)*R_M(mb)

def R_M(m):
    return (m/Msun)**0.57*Rsun

alist,elist,a23list,e23list=[],[],[],[]
dLlist,dElist,reslist=[],[],[]
for i in range(len(hier1)):
    if hier1[i]=='nstar=4 nobj=2:  2 [[1 3] 0]':
        d=my_mod(x1[i]-x3[i],y1[i]-y3[i],z1[i]-z3[i])
        if d>r_t(m1,m3):
            vr=my_mod(vx1[i]-vx3[i],vy1[i]-vy3[i],vz1[i]-vz3[i])
            E=0.5*(m1*m3/(m1+m3))*vr**2-G*m1*m3/d
            a=-G*m1*m3/(2*E)
            J=my_cross(x3[i]-x1[i],y3[i]-y1[i],z3[i]-z1[i],vx3[i]-vx1[i],vy3[i]-vy1[i],vz3[i]-vz1[i])*m3
            e=(1-J**2/(a*G*m1*m3**2))**0.5
            if e<1 and e>0 and a*(1-e)>r_t(m1,m3):
                elist.append(e)
                alist.append(float(a)) 
                a23list.append(a23[i])
                e23list.append(e23[i])
                dLlist.append(dL[i])
                dElist.append(dE[i])
                reslist.append(res[i])
                
    if hier1[i]=='nstar=4 nobj=2:  3 [[1 2] 0]':
        d=my_mod(x1[i]-x2[i],y1[i]-y2[i],z1[i]-z2[i])
        if d>r_t(m1,m2):
            vr=my_mod(vx1[i]-vx2[i],vy1[i]-vy2[i],vz1[i]-vz2[i])
            E=0.5*(m1*m2/(m1+m2))*vr**2-G*m1*m2/d
            a=-G*m1*m2/(2*E)
            J=my_cross(x2[i]-x1[i],y2[i]-y1[i],z2[i]-z1[i],vx2[i]-vx1[i],vy2[i]-vy1[i],vz2[i]-vz1[i])*m2
            e=(1-J**2/(a*G*m1*m2**2))**0.5
            if e<1 and e>0 and a*(1-e)>r_t(m1,m2):
                elist.append(e)
                alist.append(float(a)) 
                a23list.append(a23[i])
                e23list.append(e23[i])
                dLlist.append(dL[i])
                dElist.append(dE[i])
                reslist.append(res[i])
                
    if hier1[i]=='nstar=4 nobj=2:  [[0 2] 1] 3':
        d=my_mod(x0[i]-x2[i],y0[i]-y2[i],z0[i]-z2[i])
        if d>r_t(m0,m2):
            vr=my_mod(vx0[i]-vx2[i],vy0[i]-vy2[i],vz0[i]-vz2[i])
            E=0.5*(m0*m2/(m0+m2))*vr**2-G*m0*m2/d
            a=-G*m0*m2/(2*E)
            J=my_cross(x2[i]-x0[i],y2[i]-y0[i],z2[i]-z0[i],vx2[i]-vx0[i],vy2[i]-vy0[i],vz2[i]-vz0[i])*m2
            e=(1-J**2/(a*G*m0*m2**2))**0.5
            if e<1 and e>0 and a*(1-e)>r_t(m0,m2):
                elist.append(e)
                alist.append(float(a)) 
                a23list.append(a23[i])
                e23list.append(e23[i])
                dLlist.append(dL[i])
                dElist.append(dE[i])
                reslist.append(res[i])
                
    if hier1[i]=='nstar=4 nobj=2:  [[0 3] 1] 2':
        d=my_mod(x0[i]-x3[i],y0[i]-y3[i],z0[i]-z3[i])
        if d>r_t(m0,m3):
            vr=my_mod(vx0[i]-vx3[i],vy0[i]-vy3[i],vz0[i]-vz3[i])
            E=0.5*(m0*m3/(m0+m3))*vr**2-G*m0*m3/d
            a=-G*m0*m3/(2*E)
            J=my_cross(x3[i]-x0[i],y3[i]-y0[i],z3[i]-z0[i],vx3[i]-vx0[i],vy3[i]-vy0[i],vz3[i]-vz0[i])*m3
            e=(1-J**2/(a*G*m0*m3**2))**0.5
            if e<1 and e>0 and a*(1-e)>r_t(m0,m3):
                elist.append(e)
                alist.append(float(a))      
                a23list.append(a23[i])
                e23list.append(e23[i])
                dLlist.append(dL[i])
                dElist.append(dE[i])
                reslist.append(res[i])            


final=pd.DataFrame({'a':alist,'e':elist,'a23':a23list,'e23':e23list,'dL':dLlist,'dE':dElist,'res':reslist})
final.to_csv('final25',sep=',',index=False)
