#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:06:40 2023

@author: julien
"""

import numpy as np
import matplotlib.pyplot as plt
import random as rd


X=[]
Y=[]
Z=[]
P=[]
T=[]

for i in range(3*10**2):
    phi=2*np.pi*rd.uniform(0,1)
    theta=np.pi*rd.uniform(0,1)
    P.append(phi)
    T.append(theta)
    X.append(np.sin(theta)*np.cos(phi))
    Y.append(np.sin(theta)*np.sin(phi))
    Z.append(np.cos(theta))



fig_ang, ax_ang = plt.subplots(1,1)

ax_ang.scatter(P,T)
ax_ang.grid()
ax_ang.set_xlabel("Phi")
ax_ang.set_ylabel("Theta")

fig_ang.suptitle("Simple random distribution")



fig3d_1=plt.figure()
ax3d_1=plt.axes(projection="3d")

ax3d_1.scatter3D(X,Y,Z)
ax3d_1.set_xlabel("x axis")
ax3d_1.set_ylabel("y axis")
ax3d_1.set_zlabel("z axis")


X2=[]
Y2=[]
Z2=[]
P2=[]
T2=[]
for i in range(3*10**2):
    phi=2*np.pi*rd.uniform(0,1)
    theta=np.arccos(1-2*rd.uniform(0,1))
    P2.append(phi)
    T2.append(theta)
    X2.append(np.sin(theta)*np.cos(phi))
    Y2.append(np.sin(theta)*np.sin(phi))
    Z2.append(np.cos(theta))


fig_ang2, ax_ang2 = plt.subplots(1,1)

ax_ang2.scatter(P2,T2)
ax_ang2.grid()
ax_ang2.set_xlabel("Phi")
ax_ang2.set_ylabel("Theta")

fig_ang2.suptitle("Random distribution 2")



fig3d_2=plt.figure()
ax3d_2=plt.axes(projection="3d")

ax3d_2.scatter3D(X2,Y2,Z2)
ax3d_2.set_xlabel("x axis")
ax3d_2.set_ylabel("y axis")
ax3d_2.set_zlabel("z axis")



A11 = [z1 for z1 in Z if z1>1/2]
A12 = [z1 for z1 in Z if 0<z1<1/2]
A13 = [z1 for z1 in Z if -1/2<z1<0]
A14 = [z1 for z1 in Z if z1<-1/2]
    
        
A21 = [z2 for z2 in Z2 if z2>1/2]
A22 = [z2 for z2 in Z2 if 0<z2<1/2]
A23 = [z2 for z2 in Z2 if -1/2<z2<0]
A24 = [z2 for z2 in Z2 if z2<-1/2]


NBR11 = len(A11) /len(Z)
NBR12 = len(A12) /len(Z)
NBR13 = len(A13) /len(Z)
NBR14 = len(A14) /len(Z)

NBR21 = len(A21) /len(Z2)
NBR22 = len(A22) /len(Z2)
NBR23 = len(A23) /len(Z2)
NBR24 = len(A24) /len(Z2)


print(NBR11)
print(NBR12)
print(NBR13)
print(NBR14)
print("\n")
print(NBR21)
print(NBR22)
print(NBR23)
print(NBR24)