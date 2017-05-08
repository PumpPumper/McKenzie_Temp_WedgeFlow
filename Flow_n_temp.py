#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 13:09:36 2017

@author: JoeRippke

This module combines slab thermal profile work I did earlier this semester with
inaccurate mantle wedge flow around the subducting slab.

I created a 3d numpy array to record the x and z velocity components of the
mantle flow, and I populated that array with values calculated from equations
published in McKenzie (1969). The plotted streamlines are clearly not correct.

This script was written in Python3 and tested with iPython on the command line.

"""
# module that contains a slab temperature calculating function
import slabtemp as st
import numpy as np
import matplotlib.pyplot as plt

# Set Values

rho = 3 # g/cm^3               set value of rock density
Cp = 0.35 # cal/(g °C)         set value of heat capacity
vx = 10 # cm/yr                set value of slab velocity
l = 50 # km                    set value of slab thickness
kappa = 0.01 # cal/(cm °C s)   set value of thermal conductivity
xmax = 500 # km               set maximum slab length
dip = 35 # °                   set slab dip angle

Tprime = np.ones([150,350])
WedgeFlow = np.zeros([150,350,2])
diprads = (dip)*np.pi/180 # convert slab dip degrees to radians
slope = np.tan(diprads) # calculate slope
sax = 27 # x intercept bottom of slab
sbx = 114 # x intercept top of slab
sx = 114

for x in range(350):
    for z in range(150):

        # loop over each x and z coordinate
        # calculate Tprime for each location
        
        # line equation for the bottom of the slab
        zl1 = x*slope-80
        
        # define the location of the slab
        if z >= zl1 and z <= zl1+50/np.cos(diprads):
            Tprime[z,x] = 0
        else:
            Tprime[z,x] = 1

        # calculate Tprime for each location
        if Tprime[z,x] == 0:
            start = [0,-80,0]
            end = [500,500*slope-80,0]
            dist = st.pnt2line((x,z,0),start,end)
            hx = z/np.sin(diprads)
            hz = dist[0]
            Tprime[z,x] = st.slabtemp(hx,xmax,rho,Cp,vx,l,kappa,hz)

        # don't allow values greater than 1
        if Tprime[z,x] > 1:
            Tprime[z,x] = 1

        # Calculate the x,z velocity components of the wedge flow
        # Region A (below slab)
        if z > zl1+50/np.cos(diprads):
            ra = np.sqrt((x-sax)**2 + z**2)
            thetaa = 0
            if x > sax:
                thetaA = np.arcsin(z/ra) + np.pi/2
            elif x < sax:
                thetaA = np.arcsin(z/ra)
            thetaa = (180-35)*np.pi/180 + np.pi/2
            Vra = -vx*(((thetaA*np.cos(thetaA-thetaa) + np.sin(thetaA))/(thetaa + np.sin(thetaa))) - \
                       ((np.cos(thetaa) + 1)*((thetaa-thetaA)*np.sin(thetaA) - thetaA*np.sin(thetaA-thetaa)))/ \
                       (thetaa + np.sin(thetaa))**2)
            Vthetaa = (vx/(thetaa + np.sin(thetaa)))*((thetaa-thetaA)*np.sin(thetaa) \
                                                      + thetaa*np.sin(thetaa-thetaA))
            WedgeFlow[z,x,0] = Vra*np.cos(thetaA) - Vthetaa*np.sin(thetaA)
            WedgeFlow[z,x,1] = Vra*np.sin(thetaA) + Vthetaa*np.cos(thetaA)



        # Region B (above slab)
        if z < zl1: 
            rb = np.sqrt((x-sbx)**2 + z**2)
            thetaB = np.arcsin(z/rb)
            thetab = 35*np.pi/180
            Vrb = vx*((thetaB*np.sin(thetaB-thetab)+np.sin(thetaB)*np.sin(thetab) - \
                       thetaB*thetab*np.cos(thetaB-thetab) + (thetab-thetaB)*np.sin(thetaB)*np.cos(thetab))/ \
                      (thetab**2 - np.sin(thetab)**2) - \
                      ((thetaB*thetab*np.sin(thetaB-thetab)+np.sin(thetaB)*np.sin(thetab))* \
                       (2*thetab - 2*np.sin(thetab)*np.cos(thetab)))/(thetab**2 - np.sin(thetab)**2)**2)
            Vthetab = (-vx/(thetab**2 - np.sin(thetab)**2)) * \
                      ((thetab-thetaB)*np.sin(thetab)*np.sin(thetaB) - \
                       thetab*thetaB*np.sin(thetab-thetaB))
            WedgeFlow[z,x,0] = Vrb*np.cos(thetaB) - Vthetab*np.sin(thetaB)
            WedgeFlow[z,x,1] = Vrb*np.sin(thetaB) + Vthetab*np.cos(thetaB)


# Create the figure
dist = np.linspace(0,350,350)
depth = np.linspace(0,150,150)
fig = plt.figure()
strm = plt.streamplot(dist, depth, WedgeFlow[:,:,0], WedgeFlow[:,:,1], color = 'k')
im = plt.imshow(Tprime, cmap=plt.get_cmap('hot'),vmin=0,vmax=1)
Tslab = plt.axes().set_aspect('equal')
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
cbar_ax.set_xlabel('T (°C)', labelpad=20)
fig.colorbar(im, cax=cbar_ax)
plt.axes().set_ylabel('Depth (km)')
plt.axes().set_xlabel('Distance (km)')
plt.axes().set_title('Temperature Profile of a Subducting Slab \n and Mantle Flow Around the Slab')
plt.show()
