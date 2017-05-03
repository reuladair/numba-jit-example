#!/usr/bin/env python

## 
## Created Friday, April 15, 2016 at 22:01:00 reuladair (spinor@yahoo.com)
##
import time
import math
from math import *
import numpy as np
import matplotlib.pyplot as plt

import numba
from numba import jit, int32, float64, float32, void

## atmospheric temperature as a function of altitude...
@jit(float64(float64), nopython=True)
def T(h):
    if h >= 25000:                 # upper stratosphere
       Tofh = -131.21 + 0.00200 * h
    elif h >= 11000 and h < 25000: # lower stratosphere
       Tofh = -56.46
    elif h >= 0 and h < 11000:     # troposphere
       Tofh = 15.04 - 0.00649 * h
    else:
       Tofh = 0.0
    return Tofh

## atmospheric pressure as a function of altitude...
@jit(float64(float64), nopython=True)
def P(h):
    if h >= 25000:                 # upper stratosphere
       pofh = 2.488 * ( ( T(h) + 273.1 )/216.6 )**(-11.388)
    elif h >= 11000 and h < 25000: # lower stratosphere
       pofh = 22.65 * exp(1.73 - 0.000157 * h)
    elif h >= 0 and h < 11000:     # troposphere
       pofh = 101.29 * ( ( T(h) + 273.1 )/288.08 )**5.256
    else:
       pofh = 0.0
    return pofh

## atmospheric density as a function of altitude
@jit(float64(float64), nopython=True)
def rho(h):
    return P(h) / ( 0.2869 * ( T(h) + 273.1) )

## compute drag force (after Wikipedia) needs mass and cross-sectional area...
@jit(float64(float64, float64, float64, float64, float64, float64), nopython=True)
def get_dfactor(x,y,vx,vy,area,mass):
    cd = 0.25 
    v2 = vx * vx + vy * vy
    iv = 1.0 / sqrt(v2)
    mag = 0.5 * rho(y) * v2 * cd * area
    return iv * mag / mass

## acceleration of gravity...
@jit(float64(float64, float64), nopython=True)
def get_fgy(y,mass):
    g0 = 9.8000
    Re = 6371000.0 # m - radius of the earth
    return - g0 * (Re/(Re+y)) ** 2.0

## driver module that actually does the numerical integration
@jit
def driver(x0,y0,v0,thetad, mass, drag=True, init=True, dt=1.0e-5, bore=0.25):

    if (init == True):
        xi  = []
        yi  = []
        vxi = []
        vyi = []
        ti  = []
        begin = False
        theta = thetad * pi / 180.0
        vx0   = v0 * cos(theta)
        vy0   = v0 * sin(theta)
        area  = pi * (bore/2)**2.0
        print("nx0=%f, y0=%f, v0=%f, thd=%f, bore=%f, mass=%f, vx0=%f, vy0=%f, area=%f" % (x0,y0,v0,thetad,bore,mass,vx0,vy0,area))
        x  = x0
        y  = y0
        vx = vx0
        vy = vy0
        t  = 0
       
    while True:
       y += vy * dt
       if y > 0 :
         t += dt
         x += vx * dt
         fgx = 0.0
         fgy = get_fgy(y,mass)
         coeff = get_dfactor(x,y,vx,vy,area,mass)
         fdx = -vx * coeff
         fdy = -vy * coeff
         vx += dt * (fgx + fdx)
         vy += dt * (fgy + fdy)
         xi.append(x)
         yi.append(y)
         vxi.append(vxi)
         vyi.append(vyi)
         ti.append(t)
       else:
         return (xi,yi)

## main program...     
if __name__ == "__main__":
   print('begin...')
   start_time = time.time()
   # So we run the model several times now for cases where the projectile mass is
   # 125, 250, 500, 1000, 2000, 4000, and 8000 kgm respectively...
   curves = []
   masses = [1., 2., 5., 10., 20., 50., 100., 200., 500., 1000., 2000., 5000.]
   start_time = time.time()
   for mass in masses:
       print("computing: mass=%7.2f kgm" % mass)
       (thisx, thisy) = driver(0., 0., 330., 45., mass, drag=True, init=True)
       curves.append((mass, thisx, thisy))
   stop_time = time.time()
   print("WITH jit: %9.3f secs run-time" % (stop_time-start_time))

   # draw and annotate the curves...
   cnt = 0
   plt.figure(figsize=(24,6))
   for (mass, x, y) in curves:
       # draw the curves...
       plt.plot(x,y)
       # compute a shifting location (as a function of curve number) along
       # the curves for the projectile mass to be written...
       ixloc = int((0.95-0.02*cnt)*len(x))
       iyloc = int((0.95-0.02*cnt)*len(y))
       # now print the mass labels
       plt.text(x[ixloc], y[iyloc], ("%d kg" % int(mass)))
       # increment the counter so that we keep shifting the i{x,y}loc values
       cnt += 1
           
   plt.title("Trajectories for Cannon Shells Allowing for Drag")
   plt.xlabel("Down Range Distance (meters)")
   plt.ylabel("Altitude (meters)")
   plt.savefig("with.png")
   #plt.show()
   plt.close()   
