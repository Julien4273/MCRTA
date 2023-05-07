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












