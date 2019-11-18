#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  4 16:20:49 2017

@author: jacquetpierreedouard
"""

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy import random
import scipy.stats as stats

k_B=1.38064852*10**(-23) #m2 kg s-2 K-1     Boltzman Constant
radiuscontainer=10 #m    
masscontainer=10**10  #kg
massH=1    #kh
radiusH=1
def mag(x): #returns magnitude of a 2D Vector
    return np.sqrt(x[0]*x[0]+x[1]*x[1])
def unit(x): #returns unit vector along a given 2D vector
    return x/mag(x)


class Atom:
    container= plt.Circle((0,0), 10, ec = 'b', fill = False, ls = 'solid')
    def __init__(self, mass, radius, r=[0,0], v=[0,0],clr='r'):
         self.P=[]
         self.__m= mass 
         self.__radius = radius
         self.__r=np.array(r, dtype=float)
         self.__v = np.array(v, dtype=float)
         if self.__m>=10:  # only the container has such a mass, returns a big circle 
             self.__patch = plt.Circle(self.__r, radiuscontainer, fill=False)
         else:             # returns the circle associated to our 'Atoms' (atoms)
             self.__patch = plt.Circle(self.__r, self.__radius, fc=clr)
         
     
    
    def __repr__(self):
         return "Atom(%g, %g, %s, %s)" % (self.__m, self.__radius, self.__r, self.__v)
    def __str__(self):
         return "Atom(%g, %g, %s, %s)" % (self.__m, self.__radius, self.__r,self.__v)     
     
    def pos(self):
         return  self.__r  
     
    def mass(self):
        return self.__m
     
    def radius(self):
        return self.__radius
     
    def relpos(self, other):
         dr=self.__r-other.__r
         return mag(dr)
         
    def vel(self):
         return  self.__v   
     
    def KE(self): #return kinetic energy of a atom
         E=0.5*self.__m*(mag(self.__v)**2)
         return E
    
    def speed(self): 
         return mag(self.__v)
     
    def mom(self): #returns vector momentum
         return self.__m*self.__v
     
    def absmom(self): #return scalar momentum 
         return self.__m*mag(self.__v)
    
    def time_to_collision(self, other): # calculate and returns the time to collision of 2 atoms
         dr=self.__r-other.__r
         dv=self.__v-other.__v
         
         c=np.dot(dr,dr)-(self.__radius+other.__radius)**2
         b=2*np.dot(dr, dv)
         a=np.dot(dv, dv)
         
         d=b**2-4*a*c
         
         t0=(-b+np.sqrt(d))/(2*a)
         t1=(-b-np.sqrt(d))/(2*a)
         
         t=float()
         if b<0:
             t=t0
         else:
             t=t1
         return "time to collision is t=%s" % (t)
     
    def collide(self, other): #makes the necessary velocity changes after a collision between 2 atoms
        dr = self.__r - other.__r
        dv=self.__v-other.__v
        if mag(dr)<=self.__radius+other.__radius and np.dot(dr , dv)<0:
            self.__newv = self.__v - (2*other.__m)*(np.dot(self.__v-other.__v,self.__r-other.__r))*(self.__r-other.__r)/((self.__m+other.__m)*(mag(self.__r-other.__r))**2)
            other.__newv = other.__v - (2*self.__m)*(np.dot(other.__v-self.__v,other.__r-self.__r))*(other.__r-self.__r)/((self.__m+other.__m)*(mag(other.__r-self.__r))**2)
            
            self.__v=self.__newv
            other.__v=other.__newv
            
       
           
    def collidecontainer(self): #makes the necessary velocity changes after a collision with the container
        if mag(self.__r)>= radiuscontainer-self.__radius and np.dot(self.__r, self.__v)>0 :
            self.__newv = self.__v - 2*(np.dot(self.__v,self.__r)*(self.__r))/((mag(self.__r)**2))
            self.dp=self.__m*mag(self.__newv-self.__v)
            self.__v=self.__newv
        else: 
            self.dp=0

    def get_patch(self): #returns the patches of an atom to implement it to the animation later
        return self.__patch 
    
    def move(self, dt): 
         
         self.__r += dt*self.__v 
         self.__patch.center = self.__r
   
class Gas:
    
    def __init__(self):
        self.V=[]
        self.T=[]
        self.deltaP=0
        self.sumdp=[]
        self.Sumdp=[]
        self.MeanT=float()
        self.dt=0.05
        self.Maxspeed=5
        self.__Atom1 = Atom(massH,radiusH,[-0,0],[1,0])
        self.__Atom2 = Atom(massH,radiusH,[0,2.5],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom3 = Atom(massH,radiusH,[0,5],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom4 = Atom(massH,radiusH,[0,7.5],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom5 = Atom(massH,radiusH,[0,-2.5],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom6 = Atom(massH,radiusH,[0,-5],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom7 = Atom(massH,radiusH,[0,-7.5],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom8 = Atom(massH,radiusH,[2.5,0],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom9 = Atom(massH,radiusH,[5,0],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom10 = Atom(massH,radiusH,[7.5,-0],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom11 = Atom(massH,radiusH,[-2.5,0],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom12 = Atom(massH,radiusH,[-5,0],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom13 = Atom(massH,radiusH,[-7.5,0],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom14 = Atom(massH,radiusH,[2,2],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom15 = Atom(massH,radiusH,[4,4],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom16 = Atom(massH,radiusH,[6,6],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom17 = Atom(massH,radiusH,[-2,-2],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom18 = Atom(massH,radiusH,[-4,-4],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom19 = Atom(massH,radiusH,[-6,-6],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom20 = Atom(massH,radiusH,[2,-2],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom21 = Atom(massH,radiusH,[4,-4],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom22 = Atom(massH,radiusH,[6,-6],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom23 = Atom(massH,radiusH,[-2,2],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom24 = Atom(massH,radiusH,[-4,4],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__Atom25 = Atom(massH,radiusH,[-6,6],[np.random.uniform(-self.Maxspeed,self.Maxspeed), np.random.uniform(-self.Maxspeed,self.Maxspeed)])
        self.__container= Atom(masscontainer,radiuscontainer, [0,0],[0,0])
        self.text0=None
        self.B=[self.__Atom1, self.__Atom2, self.__Atom3,self.__Atom4,self.__Atom5,self.__Atom6,self.__Atom7,self.__Atom8,self.__Atom9,self.__Atom10, self.__Atom11, self.__Atom12,self.__Atom13,self.__Atom14,self.__Atom15,self.__Atom16,self.__Atom17,self.__Atom18,self.__Atom19,self.__Atom20,self.__Atom21,self.__Atom22,self.__Atom23,self.__Atom24,self.__Atom25]
        
        self.text0=None
        
    def InternalEnergy(self):  #return the total internal energy(since the gas is kinetic) of the system at time t 
        ke=[]
        for e in self.B:
            ke.append(e.KE())
        ET=float(sum(ke))
        return ET
     
     
    def MeanKinE(self): #returns the mean kinetic energy per atom at time t 
         mke=self.InternalEnergy() /len(self.B)
         return mke
     
    def Temp(self): #returns the temperature of the system at time t 
        T=self.MeanKinE()/k_B
        return T
    
    def Vel(self): #for a given framenumber, returns a list of all speeds of individual components 
        for v in self.B:
            self.V.append(v.speed())
        return self.V
         
    def init_figure(self): #defines the initial frame
        ax = plt.axes(xlim=(-10, 10), ylim=(-10, 10))
        ax.add_artist(Atom.container)
        self.__text0 = ax.text(-9.9,9,"f={:4d}".format(0,fontsize=12))
        patches = [self.__text0]
        for c in self.B:
            pch = c.get_patch()
            ax.add_patch(pch)
            patches.append(pch)
        
        return patches
        
        
    def next_frame(self, i): # define the state of the next frame 
            self.__text0.set_text("t={:4d}".format(i))
            patches = [self.__text0]
            for d in self.B:
                d.move(self.dt)
                d.collidecontainer()
                for d2 in self.B:
                    if not d2==d:
                        d.collide(d2)
                patches.append(d.get_patch())

            return patches
        
    def next_frame2(self,i): #same as previous but used to plot graph (without the text0 bit)
         for d in self.B:
             d.move(self.dt)
             d.collidecontainer()
             if d.dp!=0:
                 self.deltaP+=d.dp
             for d2 in self.B:
                 if not d2==d:
                        d.collide(d2)
         self.T.append(self.Temp())  
         self.MeanT=sum(self.T)/len(self.T)
         self.sumdp.append(self.deltaP)
         
         self.MeanKinE()
         self.Vel() 
         
    def Pressure(self, a): #returns pressure for a given frame number 
        P=[]
        for i in range(a):
            gas.next_frame2(i)
            t=i*gas.dt
            if t!=0:
                P.append(self.deltaP/(t*20*np.pi))   
        return P   
     
#booleans act as switch ! If you want to activate a particular switch ,all the other to False et the desired one on True  

Animate=False     #basic animation
if Animate==True:     
    fig = plt.figure()
    ax = plt.axes(xlim=(-10, 10), ylim=(-10, 10))
    ax.axes.set_aspect('equal') 
    movie = Gas()
    a=movie.next_frame
    anim = animation.FuncAnimation(fig,
                                   movie.next_frame, 
                                   init_func = movie.init_figure, 
                                   interval = 20,
                                   blit = True)
    
    plt.show()




#plot internal Energy vs Time
PlotInternalEnergy=False   
if PlotInternalEnergy==True:
    gas=Gas()
    A=[]
    B=[]
    fig=plt.figure()
    for i in range (50):
        gas.next_frame2(i)
        A.append(gas.InternalEnergy())
        print (gas.InternalEnergy())
        B.append(gas.MeanKinE())
    x=np.arange(0, 2.5, 0.05)
    ax.axes.set_aspect('equal') 
    ax = plt.axes(xlim=(0, 2.5), ylim=(0, 300))
    X=np.arange(len(A))
    x=X*gas.dt
    plt.plot(x, A, '.')
    plt.xlabel('time (s)')
    plt.ylabel('Internal Energy (J)')
    plt.title('Internal Energy vs Time')
    plt.grid()
    plt.show()
    




#Plot Pressure vs Time
PlotPressure=False  
if PlotPressure==True: 
    gas=Gas()
    a=1000 #time over which we calculate the variations
    fig=plt.figure()
    A=gas.Pressure(a)
    X=np.arange(len(A))
    x=X*gas.dt
    plt.plot(x, A, '.', lw=1)
    plt.xlabel('time (s)')
    plt.grid()
    plt.ylabel('Pressure (Pa)')
    plt.title('Pressure vs Time')
    plt.show()
 


   
PlotTemperature=False   
if PlotTemperature==True:
    gas=Gas()
    T1=[]
    fig=plt.figure()
    for i in range (50):
        gas.next_frame2(i)
        T1.append(gas.Temp())
        print (T1)
    ax.axes.set_aspect('equal') 
    ax = plt.axes(xlim=(0, 2.5), ylim=(0, 1E24))
    X=np.arange(len(T1))
    x=X*gas.dt
    plt.plot(x, T1, '.', lw=1)
    plt.xlabel('time (s)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature vs Time')
    plt.grid()
    plt.show()
   
PlotBoltzmann=True    
if PlotBoltzmann==True: 
    gas=Gas()
    for i in range (20000):
        gas.next_frame2(i)
    maxwell=stats.maxwell
    params=maxwell.fit(gas.V)
    Mean=np.mean(gas.V)
    Var=np.var(gas.V)
    x=np.linspace(0,15,80)
    plt.hist(gas.V, bins=100 , normed=True)
    plt.plot(x, maxwell.pdf(x, *params), lw=1)
    plt.xlabel('Speed (m/s')
    plt.ylabel('number of occurence')
    plt.title('Normalized Speed Distribution')
    plt.grid()
    plt.show()
    print('Mean velocity', np.mean(gas.V))
    print('Variance', np.var(gas.V))
    
