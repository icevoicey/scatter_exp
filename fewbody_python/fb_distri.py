#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 18:30:24 2017

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


data=pd.read_csv('final',sep=',')


e0=np.array(data['e'])
a0=np.array(data['a'])
a230=np.array(data['a23'])
e230=np.array(data['e23'])
dL0=np.array(data['dL'])
dE0=np.array(data['dE'])
res0=np.array(data['res'])

a,e,a23,e23,dL,dE,res=[],[],[],[],[],[],[]

for i in range(len(a0)):
    a.append(a0[i])
    e.append(e0[i])
    a23.append(a230[i])
    e23.append(e230[i])
    dL.append(dL0[i])
    dE.append(dE0[i])
    res.append(res0[i])
    

a_AU=[a[i]/AU for i in range(len(a))]    
a_AU_log=[np.log10(a[i]/AU) for i in range(len(a))]

a23_AU=[a23[i]/AU for i in range(len(a23))]
a23_AU_log=[np.log10(a23[i]/AU) for i in range(len(a23))]
e23_2=[e23[i]**2 for i in range(len(e23))]

a_pc=[a[i]/pc for i in range(len(a))]    
a_pc_log=[np.log10(a[i]/pc) for i in range(len(a))]
e_2=[e[i]**2 for i in range(len(e))]

rmax_AU_log=[np.log10(a[i]*(1+e[i])/AU) for i in range(len(a))]
rmin_AU_log=[np.log10(a[i]*(1-e[i])/AU) for i in range(len(a))]

#plt.figure()        
#ax1=plt.subplot(111)
#ax1.hist(a23_AU_log,bins=len(a23)//30)
#ax1.set_xlabel('log (a23/AU)')
#ax1.set_ylabel('N')

plt.figure()
ax1=plt.subplot(211)
na,bina,pata=ax1.hist(a_AU_log,bins=len(a)//20)
ax1.set_xlabel('log (a /AU)')
ax1.set_ylabel('N')

ax2=plt.subplot(212)
nr,binr,patr=ax2.hist(e,bins=len(e)//20)
ax2.set_xlabel('e')
ax2.set_ylabel('N')
        
#plt.figure()
#ax1=plt.subplot(211)
#na,bina,pata=ax1.hist(rmax_AU_log,bins=len(a)//30)
#ax1.set_xlabel('log (rmax /AU)')
#ax1.set_ylabel('N')
#
#ax2=plt.subplot(212)
#nr,binr,patr=ax2.hist(rmin_AU_log,bins=len(a)//30)
#ax2.set_xlabel('log (rmin /AU)')
#ax2.set_ylabel('N')
#
#a_pre_log=[np.log10(0.56*(5e6)**(2/3)*a23[i]/AU) for i in range(len(a23))]
#
#plt.figure()
#plt.hist(a_pre_log,bins=len(a)//10)

#plt.figure()
#plt.hist(a_AU,bins=len(a)//30)

