#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 13:09:36 2017

@author: JoeRippke


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
                thetaa = (np.pi/2) + np.arccos(z/ra)
            elif x < sax:
                thetaa = np.arcsin(z/ra)
            else:
                thetaa = np.pi/2
            thetaA = (180-35)*np.pi/180
##            Vra = -vx*((thetaa-thetaA)*np.cos(thetaA) - thetaA*np.cos(
##                thetaA-thetaa) - np.sin(thetaA) - np.sin(thetaA-thetaa))/(
##                    thetaa + np.sin(thetaa))
            Vra = (-vx/(thetaa+np.sin(thetaa)))*(-np.sin(thetaA-thetaa) + \
                                                 (thetaa - thetaA)*np.cos(thetaA) \
                                                 - thetaA*np.cos(thetaA*thetaa) \
                                                 - np.sin(thetaA))
##            Vthetaa = vx*((thetaa-thetaA)*np.sin(thetaA) + thetaA*np.sin(
##                thetaa - thetaA))/(thetaa + np.sin(thetaa))
            Vthetaa = (vx/(thetaa + np.sin(thetaa)))*((thetaa-thetaA)*np.sin(thetaa) \
                                                      + thetaa*np.sin(thetaa-thetaA))
            WedgeFlow[z,x,0] = Vra*np.cos(thetaa) - Vthetaa*np.sin(thetaa)
            WedgeFlow[z,x,1] = Vra*np.sin(thetaa) + Vthetaa*np.cos(thetaa)

        # Region B (above slab)
        if z < zl1: 
            rb = np.sqrt((x-sbx)**2 + z**2)
            thetab = np.arcsin(z/rb)
            thetaB = diprads
##            Vrb = vx*((thetab*np.sin(thetab-thetaB) - np.sin(thetaB)*np.sin(
##                thetab) + thetaB*thetab*np.cos(thetaB-thetab) + (
##                    thetab-thetaB)*np.cos(thetaB)*np.sin(thetab)
##                       )/thetab**2 - np.sin(thetab)**2)
            Vrb = (vx/(thetab**2 - np.sin(thetab)**2))* \
                  (thetab*np.sin(thetaB-thetab) - np.sin(thetaB)*np.sin(thetab) + \
                   thetaB*thetab*np.cos(thetaB-thetab) + (thetab-thetaB)* \
                   np.cos(thetaB)*np.sin(thetab))
##            Vthetab = -vx*((thetab-thetaB)*np.sin(thetab)*np.sin(
##                thetaB) - thetab*thetaB*np.sin(thetab-thetaB)
##                           )/thetab**2 - np.sin(thetab)**2
            Vthetab = (-vx/(thetab**2 - np.sin(thetab)**2)) * \
                      ((thetab-thetaB)*np.sin(thetab)*np.sin(thetaB) - \
                       thetab*thetaB*np.sin(thetab-thetaB))
            WedgeFlow[z,x,0] = Vrb*np.cos(thetab) - Vthetab*np.sin(thetab)
            WedgeFlow[z,x,1] = Vrb*np.sin(thetab) + Vthetab*np.cos(thetab)

dist = np.linspace(0,350,350)
depth = np.linspace(0,150,150)

strm = plt.streamplot(dist, depth, WedgeFlow[:,:,0], WedgeFlow[:,:,1], color = 'k')
im = plt.imshow(Tprime, cmap=plt.get_cmap('hot'),vmin=0,vmax=1)
Tslab = plt.axes().set_aspect('equal')
#plt.gca().set_ylim([150,0])
plt.colorbar()
plt.show()
