#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 10:51:46 2023

@author: julien
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import random as rd
from prettytable import PrettyTable
import gc



#%% Global variables and necessary functions

markers = ["." , "," , "o" , "v" , "^" , "<", ">", "*"]
colors = ['r','g','b','c','m', 'y', 'k', 'o']

x_lim, y_lim, z_lim = 5, 5, 5
r_0 = 10**-3                    # Radius of scaterers
#n = (10**5*r_0)**-3             # Low density
n = (10**2*r_0)**-3            # Medium density
#n = (2*r_0)**-3                # High density

dx, dy, dz = 0.01, 0.01, 0.01     # Steps

sig = 4 * np.pi * r_0           # Cross section

nc = 100                        # Number of cells
nl = 100                        # Number of layers

nbr_photon = 10

field_names = ["Photon", "Interactions", "Compton", "Rayleigh", "Absorbed", "Mean dist trav", "wl_i", "wl_f"]
photon_table = PrettyTable(field_names)


def choose_cell(x_pos : float, y_pos : float, z_pos : float, nbr_cell : int):
    
    c_x, c_y, c_z = x_lim/nbr_cell, y_lim/nbr_cell, z_lim/nbr_cell
    
    i,j,k = abs(int(x_pos/c_x)),abs(int(y_pos/c_y)),abs(int(z_pos/c_z))
    
    if i > nbr_cell+1 :
        return nbr_cell+1, j, k
    elif j > nbr_cell+1 : 
        return i, nc+1, k
    elif k > nbr_cell+1 :
        return i, j, nbr_cell+1
    elif i > nbr_cell+1 and j > nbr_cell+1 :
        return nbr_cell+1, nbr_cell+1, k
    elif i > nbr_cell+1 and k > nbr_cell+1 :
        return nbr_cell+1, j, nbr_cell+1
    elif j > nbr_cell+1 and k > nbr_cell+1 :
        return i, nbr_cell+1, nbr_cell+1
    
    return i, j, k



def rand_dens(tab : list):
    
    sum_rd=0
    
    for i in range(1,nc+1):
        for j in range(1,nc+1):
            for k in range(1,nc+1):
                rdn=rd.random()
                tab[i][j][k] = rdn * n
                sum_rd += tab[i][j][k] / n
    
    for i in range(1,nc+1):
        for j in range(1,nc+1):
            for k in range(1,nc+1):
                
                tab[i][j][k] /= sum_rd
    print("\n")
    print("Cube initialized")
    print("\n")



def spherical_density(nbr_layer : int):
    
    RDN = [rd.random() for i in range(nbr_layer)]
    S = 0
    for x in RDN:
        S+=x
    RDN_U = [x/S for x in RDN]
    RDN_U.sort()
    Layers = [x * n for x in RDN_U]
    return Layers
        


def choose_layer(x_pos : float, y_pos : float, z_pos : float, nbr_layer : int):
    
    l_x_0, l_y_0, l_z_0 = x_lim/2, y_lim/2, z_lim/2
    l_x, l_y, l_z = x_lim/nbr_layer, y_lim/nbr_layer, z_lim/nbr_layer
    
    r_layer = dist(l_x_0, l_y_0, l_z_0, l_x, l_y, l_z)
    r_photon = dist(l_x_0, l_y_0, l_z_0, x_pos, y_pos, z_pos)
    
    i = abs(int(r_photon/r_layer))
    
    if i > nbr_layer :
        return nbr_layer - 1
    
    return i



def choose_density(geo : str, x_pos : float, y_pos : float, z_pos : float, nbr_cell : int, nbr_layer : int, tbn_cell : list, tbn_layer : list):
    
    if geo == "cell" :
        
        # Collect indices for table_density corresponding to density where photon is
        I = [*choose_cell(x_pos, y_pos, z_pos, nbr_cell)]
        n_cell = tbn_cell[I[0]][I[1]][I[2]]
        return n_cell
    
    elif geo == "layer" :
        
        i = choose_layer(x_pos, y_pos, z_pos, nbr_layer)
        n_layer = tbn_layer[i]
        return n_layer
    
    raise(ValueError("Geometry does not exist. Please choose between 'cell' or 'layer'."))



def dist(x_i : float, y_i : float, z_i : float, x_f : float, y_f : float, z_f : float):
    L=np.sqrt((x_f - x_i)**2 + (y_f - y_i)**2 + (z_f - z_i)**2)
    return L



def find_best_scat(pos_xf : float, pos_yf : float, pos_zf : float, last_tau : float, ndx : float, ndy : float, ndz : float, tau : float, geo : str, nbr_cell : int, nbr_layer : int, tbn_cell : list, tbn_layer : list, epsilon : float) :
    
    delta = 1
    tau_b = last_tau
    x_b, y_b, z_b = pos_xf - ndx, pos_yf - ndy, pos_zf - ndz
    
    while abs(tau_b - tau) > epsilon :
        
        x_b, y_b, z_b = pos_xf - ndx, pos_yf - ndy, pos_zf - ndz
        n_0 = choose_density(geo, x_b, y_b, z_b, nbr_cell, nbr_layer, tbn_cell, tbn_layer)
        tau_b = last_tau
        delta *= 10
        
        while tau_b < tau :
            
            x_b += ndx/delta
            y_b += ndy/delta
            z_b += ndz/delta
            
            n_1 = choose_density(geo, x_b, y_b, z_b, nbr_cell, nbr_layer, tbn_cell, tbn_layer)
            
            delta_l = tau / (n_1 * sig * 100 * delta)
            tau_b += (n_1 + n_0) * sig * delta_l/2
            
            n_0 = n_1
    
    return [x_b, y_b, z_b]


def plot_rt(positions : list, scatt : list ):
    fig = plt.figure(figsize = (10,5))
    ax = plt.axes(projection = "3d")        
    for i in range(nbr_photon):
        ax.plot3D(positions[i][0], positions[i][1], positions[i][2], label = "photon "+str(i+1))
        ax.scatter3D(scatt[i][0], scatt[i][1], scatt[i][2], marker = markers[-1], s = 10, c = "black")

    ax.legend(fontsize = 6, loc = 2)
    ax.set_xlabel("x axis")
    ax.set_ylabel("y axis")
    ax.set_zlabel("z axis")

    fig.suptitle("Radiative transfer")
    fig.savefig("Radiative_transfer")
    rt_plot = Image.open("Radiative_transfer.png")
    rt_plot.show()

    

def plot_tele(tm : list):
    ind_diff_0 = [i for i in range(nbr_photon) if len(tm[i][0]) != 0]
    if len(ind_diff_0) == 0:
        return 404
    
    fig_maps,ax_maps = plt.subplots(1,1)  
    for x in ind_diff_0:
        ax_maps.scatter(tm[x][0], tm[x][1], label = "photon "+str(x+1))
    ax_maps.grid()
    ax_maps.legend(fontsize = 7, loc = 0)
    ax_maps.set_xlabel("x axis")
    ax_maps.set_ylabel("z axis")

    fig_maps.suptitle("Telescope")
    fig_maps.savefig("Telescope")
    tele_plot = Image.open("Telescope.png")
    tele_plot.show()   

def plot_all(positions : list, scatt : list, tm : list):
    plot_rt(positions, scatt)
    plot_tele(tm)
    


#%% Definition of the photon class

class photon:
    
    # Constructor
    def __init__(self, mx : float, my : float, mz : float, wl : float, ref : int):
        
        phi = 2 * np.pi * rd.random()
        theta = np.arccos(1 - 2 * rd.random())
        
        self.x = mx
        self.y = my
        self.z = mz
        
        self.ndx = np.sin(theta) * np.cos(phi) * dx
        self.ndy = np.sin(theta) * np.sin(phi) * dy
        self.ndz = np.cos(theta) * dz
        
        self.wl = wl
        
        self.ref = ref


    # Destructor
    def __del__(self):
        print("Photon "+str(self.ref)+" : destroyed")
        
        
    

    # Define whether the photon is absorbed, scattered by Raygleigh or Compton
    # If not absorbed : True, otherwise : False
    def absornot(self, p_abs : float, p_comp : float):
        
        # Angle of scattering
        phi = 2 * np.pi * rd.random()
        theta = np.arccos(1 - 2 * rd.random())
        
        # Random probability
        rdn = rd.random()
        

        
        if (rdn <= p_abs):
            """print("photon "+str(self.ref)+" : absorbed")"""
            return False, 0
        
        elif (rdn <= p_comp) :
            
            wl_0 = self.wl
            
            # Calculate new wavelength
            self.wl = wl_0 * (1 + 1/wl_0 * (1-rd.random())) 
            
            # New direction
            self.ndx = np.sin(theta) * np.cos(phi) * dx
            self.ndy = np.sin(theta) * np.sin(phi) * dy
            self.ndz = np.cos(theta) * dz
            
            """print("photon "+str(self.ref)+" : Compton")"""
            return True, 1
        
        """print("photon "+str(self.ref)+" : Rayleigh")"""
        
        # New direction
        self.ndx = np.sin(theta) * np.cos(phi) * dx
        self.ndy = np.sin(theta) * np.sin(phi) * dy
        self.ndz = np.cos(theta) * dz
        
        return True, 2


    # Check whether the photon is inside the medium or not
    # If inside : True, otherwise : False
    def in_medium(self):
        
        if self.x > x_lim or self.x < 0 or self.y > y_lim or self.y < 0 or self.z > z_lim or self.z < 0:
            
            """print("Photon "+str(self.ref)+" : escaped")"""
            return False
            
        return True
    
    
    # Stock coordonates of photon inside list and return it
    def coord(self, geo: str, p_abs : float, p_comp : float, nbr_cell : int, nbr_layer : float, tbn_cell : list, tbn_layer : list):
        
        # Coordonates of the photon
        pos_photon = [[self.x],[self.y],[self.z]]
        # Position of interaction (absorbed or scattered)
        pos_scat = [[],[],[]]
        # Position on the telescope
        tele_map = [[],[]]
        
        # Starting position and wavelength of the photon
        starting_x, starting_y, starting_z = self.x, self.y, self.z
        starting_wl = self.wl
        
        # Density corresponding to position of the photon and geometry chosen
        n_0 = choose_density(geo, self.x, self.y, self.z, nbr_layer, nbr_cell, tbn_cell, tbn_layer)
        
        # Optical depth
        tau = -np.log(1 - rd.random())
        
        # Position before and after interaction
        x_i, y_i, z_i = self.x, self.y, self.z
        x_f, y_f, z_f = 0, 0, 0
        
        not_absorbed = True
        # Check whether the photon is inside the medium
        inside = self.in_medium()
        
        nbr_comp, nbr_ray, is_abs = 0, 0, "no"
        av_L = 0
        nbr_interactions = 0
        
        while (inside and not_absorbed):
            
            tau_s = 0
            
            while (tau_s < tau and inside):
                
                # Mooves the photon
                self.x += self.ndx
                self.y += self.ndy
                self.z += self.ndz
                
                # Stock coordonates
                pos_photon[0].append(self.x)
                pos_photon[1].append(self.y)
                pos_photon[2].append(self.z)
                
                
                n_1 = choose_density(geo, self.x, self.y, self.z, nbr_layer, nbr_cell, tbn_cell, tbn_layer)
                
                # Integration step
                delta_l = tau / (n_1 * sig * 100)
                
                # Stock tau for further prescision
                last_tau = tau_s
                
                tau_s += (n_1 + n_0) * sig * delta_l/2
                
                n_0 = n_1

                # Check whether the photon is still inside the medium
                inside = self.in_medium()
                
                if inside == False :
                # Photon escaped from the medium, it hit the telecope   
                    norme_r = dist(starting_x, starting_y, starting_z, self.x, self.y, self.z)
                    x_map = norme_r * np.cos(np.arctan(self.z/norme_r)) * np.sin(np.arctan(self.ndy/self.ndx)-np.arctan(self.y/self.x))
                    tele_map[0].append(x_map)
                    tele_map[1].append(self.z)
            
            if inside :
                
                x_f, y_f, z_f = self.x, self.y, self.z

                # Find where exactly the photon interacted with a precision of espilon
                epsilon = 10**-5
                Xf = find_best_scat(x_f, y_f, z_f, last_tau, self.ndx, self.ndy, self.ndz, tau, geo, nbr_cell, nbr_layer, tbn_cell, tbn_layer, epsilon)
                
                # Stock coordonates of interactions
                pos_scat[0].append(Xf[0])
                pos_scat[1].append(Xf[1])
                pos_scat[2].append(Xf[2])
                # Replace last coordonates of the photon by the new ones, more precise
                pos_photon[0][-1] = Xf[0]
                pos_photon[1][-1] = Xf[1]
                pos_photon[2][-1] = Xf[2]
                
                # Calculate the distance traveled between two interactions
                L = dist(x_i, y_i, z_i, *Xf)
                av_L += L
                
                x_i, y_i, z_i = Xf[0], Xf[1], Xf[2]
                
                # Check whether the photon is absorbed, if not continue with a new optical depth
                nbr_interactions += 1 
                not_absorbed, int_type = self.absornot(p_abs,p_comp)
                if int_type == 0:
                    is_abs = "yes"
                elif int_type == 1:
                    nbr_comp += 1
                else:
                    nbr_ray += 1
                """print("L = "+str(L))"""
                
                if not_absorbed :
                    
                    tau = -np.log(1 - rd.random())
        """
        print("\n")
        # Print ratio of originate and last value of wavelength
        print("Photon "+str(self.ref)+" : lambda_0/last_lambda : "+str(self.wl/starting_wl))
        print("\n")
        """
        row_photon = [self.ref, nbr_interactions, nbr_comp, nbr_ray, is_abs, '{:.2f}'.format(av_L/nbr_interactions), starting_wl, int(self.wl)]
        return [pos_photon, pos_scat, tele_map, row_photon]
    
        
#%% main


# Create table for cartesian geometry
Cells = [[[n/(nc**2*6+nc*16+8) for k in range(nc+2)] for j in range(nc+2)] for i in range(nc+2)]
# Change the density into random density
rand_dens(Cells)

# Create table for spherical geometry
Layers = spherical_density(nc)

     
# Create photon objects, starting from a point source in the middle of the medium

photons=[photon(x_lim/2, y_lim/2, z_lim/2, rd.randint(400,801), i+1) for i in range(nbr_photon)]

# Retrieve coordonates of photons
"""
print("Choose a geometry between cubic or spherical :")
print("For cubic, enter : cell")
print("For spherical, enter : layer")
geo = input()
"""
geo = "cell"
#geo = "layer"
to_plot=[x.coord(geo, 0.01, 0.5, nc, nl, Cells, Layers) for x in photons]

pos = [to_plot[i][0] for i in range(len(to_plot))]
scat = [to_plot[i][1] for i in range(len(to_plot))]
tele_maps = [to_plot[i][2] for i in range(len(to_plot))]
all_rows = [to_plot[i][3] for i in range(len(to_plot))]
photon_table.add_rows(all_rows)
photon_table.title = "Density : "+ str(int(np.ceil(n)))+", geometry : "+ str(geo)

# Plot photons paths and map of the telescope           
plot_all(pos, scat, tele_maps)


im = Image.new("RGB", (550, 255), "white")
draw = ImageDraw.Draw(im)
draw.text((10, 10), str(photon_table), fill="black")
im.show()


# Destroy objects
print("Destroyers called")
print("\n")
for x in photons:
    del x
gc.collect()

print("\n")
print(photon_table)





    

