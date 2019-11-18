# -*- coding: utf-8 -*-
"""
An example module
Pierre-Edouard Jacquet
"""

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy import random






radiuscontainer=10
masscontainer=10**10

def mag(x):
    return np.sqrt(x[0]*x[0]+x[1]*x[1])

def unit(x):
    return x/mag(x)

class Ball:

    container= plt.Circle((0,0), 10, ec = 'b', fill = False, ls = 'solid')
    def __init__(self, mass, radius, r=[0,0], v=[0,0],clr='r'):
         self.__m= mass 
         self.__radius = radius
         self.__r=np.array(r, dtype=float)
         self.__v = np.array(v, dtype=float)
         if self.__m>=10:
             self.__patch = plt.Circle(self.__r, radiuscontainer, fill=False)
         else:
             self.__patch = plt.Circle(self.__r, self.__radius, fc=clr)
         
     
    
    def __repr__(self):
         return "Ball(%g, %g, %s, %s)" % (self.__m, self.__radius, self.__r, self.__v)
    def __str__(self):
         return "Ball(%g, %g, %s, %s)" % (self.__m, self.__radius, self.__r,self.__v)     
     
    def pos(self):
         return  self.__r   
     
    def radius(self):
        return self.__radius
     
    def relpos(self, other):
         dr=self.__r-other.__r
         return mag(dr)
         
    def vel(self):
         return  self.__v      
     
    def time_to_collision(self, other):
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
     
    def collide(self, other):
        dr = self.__r - other.__r
        
        if mag(dr)-(self.__radius+other.__radius)<=1E-15:
            self.__newv = self.__v - (2*other.__m)*(np.dot(self.__v-other.__v,self.__r-other.__r))*(self.__r-other.__r)/((self.__m+other.__m)*(mag(self.__r-other.__r))**2)
            other.__newv = other.__v - (2*self.__m)*(np.dot(other.__v-self.__v,other.__r-self.__r))*(other.__r-self.__r)/((self.__m+other.__m)*(mag(other.__r-self.__r))**2)
            
            self.__v=self.__newv
            other.__v=other.__newv
            
       
           
    def collidecontainer(self):
        if mag(self.__r)>= 9 and np.dot(self.__r, self.__v)>=0 :    
            self.__newv = self.__v - 2*(np.dot(self.__v,self.__r)*(self.__r))/((mag(self.__r)**2))
            self.__v=self.__newv
            
    def get_patch(self):
        return self.__patch 
    
    def move(self, dt):
         
         self.__r += dt*self.__v 
         self.__patch.center = self.__r

   
class Orbits:
    
    def __init__(self):
        
        self.__ball1 = Ball(0.2,1,[0,0],[np.random.uniform(-1,1),np.random.uniform(-1,1)])
        self.__ball2 = Ball(0.2,1,[0,2.5],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball3 = Ball(0.2,1,[0,5],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball4 = Ball(0.2,1,[0,7.5],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball5 = Ball(0.2,1,[0,-2.5],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball6 = Ball(0.2,1,[0,-5],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball7 = Ball(0.2,1,[0,-7.5],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        """self.__ball8 = Ball(0.2,1,[2.5,0],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball9 = Ball(0.2,1,[5,0],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball10 = Ball(0.2,1,[7.5,-0],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball11 = Ball(0.2,1,[-2.5,0],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball12 = Ball(0.2,1,[-5,0],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball13 = Ball(0.2,1,[-7.5,0],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball14 = Ball(0.2,1,[2,2],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball15 = Ball(0.2,1,[4,4],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball16 = Ball(0.2,1,[6,6],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball17 = Ball(0.2,1,[-2,-2],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball18 = Ball(0.2,1,[-4,-4],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball19 = Ball(0.2,1,[-6,-6],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball20 = Ball(0.2,1,[2,-2],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball21 = Ball(0.2,1,[4,-4],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball22 = Ball(0.2,1,[6,-6],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball23 = Ball(0.2,1,[-2,2],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball24 = Ball(0.2,1,[-4,4],[np.random.uniform(-1,1), np.random.uniform(-1,1)])
        self.__ball25 = Ball(0.2,1,[-6,6],[np.random.uniform(-1,1), np.random.uniform(-1,1)])"""
        self.__container= Ball(masscontainer,10, [0,0],[0,0])
        self.text0=None
        self.L=[self.__ball1, self.__ball2, self.__ball3, self.__ball4, self.__ball5, self.__ball6, self.__ball7 ]   
        
        
        """self.B=[]
        self.L=[]
        self.C=[]
        self.A=[]
        for a in range(0,199):
            x1=np.random.uniform(-9,9)
            x2=np.random.uniform(-9,9)
            v1=np.random.uniform(-1,1)
            v2=np.random.uniform(-1,1)
            MyBall=Ball(1, 0.1 ,[x1, x2], [v1,v2], clr='r')
            
            if np.sqrt(x1**2 + x2**2) <= 9.9:    
                self.B.append(MyBall)
        #print(self.B)
        self.A=list(combinations(self.B, 2))
        #print(self.A)
        
        for i in range(0, len(self.A)-1):
            if (self.A[i][0]).relpos(self.A[i][1])- (self.A[i][0]).radius()+(self.A[i][1]).radius()>=1E-10:
                C1=self.A[i][0]
                C2=self.A[i][1]
                (self.C).append(C1)
                (self.C).append(C2)
                self.L=[item for item, count in collections.Counter(self.C).items() if count ==199]
                
            
        """
       
        self.__container= Ball(masscontainer,10, [0,0],[0,0])
        self.L.append(self.__container)
        self.text0=None
        
        
    def init_figure(self):
        
        ax.add_artist(Ball.container)
        # initialise the text to be animated and add it to the plot
        self.__text0 = ax.text(-9.9,9,"f={:4d}".format(0,fontsize=12))
        patches = [self.__text0]
        # add the patches for the balls to the plot
        for c in self.L:
            pch = c.get_patch()
            ax.add_patch(pch)
            patches.append(pch)
        
        return patches
        
        
    def next_frame(self, framenumber):
        
        self.__text0.set_text("f={:4d}".format(framenumber))
        patches = [self.__text0]
        for d in self.L:
            d.move(0.05)
            d.collidecontainer()
            for d2 in self.L:
                if not d2==d:
                    d.collide(d2)
                patches.append(d.get_patch())
        return patches
if __name__ == "__main__":
    
    fig = plt.figure()
    ax = plt.axes(xlim=(-10, 10), ylim=(-10, 10))
    ax.axes.set_aspect('equal') 
    movie = Orbits()
    
    anim = animation.FuncAnimation( fig,
                                    movie.next_frame, 
                                    init_func = movie.init_figure, 
                                    #frames = None, 
                                    interval = 20,
                                    blit = True)
    plt.show()