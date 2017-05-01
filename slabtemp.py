#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 16:33:36 2017

@author: JoeRippke

This module exists to contain functions needed to calculate slab temperature.
In the future I may add functions to calculate other slab attributes.
"""
import numpy as np
import vectors as vec

def slabtemp(x,xm,rho,Cp,vx,l,kappa,z):
    v1 = vx / 31536000 #           convert cm/yr to cm/sec
    l1 = l * 100000 #              convert km to cm
    xp = x/xm #      nondimensionalize x
    zp = z/l #                 nondimensionalize z
    R = (rho*Cp*v1*l1)/(2*kappa) #  calculate Reynolds number
    Tprime = 1-(2/np.pi)*np.exp(-1*((np.pi**2)*xp)/(2*R))*np.sin(np.pi*zp)
    return Tprime

def lithospheretemp(x,xm,rho,Cp,vx,l,kappa,z):
    v1 = vx / 31536000
    l1 = l * 100000
    zp = l/(l-z)
    xp = x/xm
    R = (rho*Cp*v1*l1)/(2*kappa)
    Tprime = 1-(2/np.pi)*np.exp(-1*((np.pi**2)*xp)/(2*R))*np.sin(np.pi*zp)
    return Tprime

def pnt2line(pnt,start,end):
    line_vec = vec.vector(start,end)
    pnt_vec = vec.vector(start,pnt)
    line_len = vec.length(line_vec)
    line_unitvec = vec.unit(line_vec)
    pnt_vec_scaled = vec.scale(pnt_vec, 1.0/line_len)
    t = vec.dot(line_unitvec,pnt_vec_scaled)    
    if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0
    nearest = vec.scale(line_vec,t)
    dist = vec.distance(nearest,pnt_vec)
    nearest = vec.add(nearest,start)
    return (dist,nearest)

