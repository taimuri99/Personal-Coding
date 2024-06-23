import math
import numpy as np
import matplotlib.pyplot as plt
from turtle import *
from scipy.signal import correlate2d
from os import system

G =   1 # chose a value between 1-10
dt =  0.1
potential_planets = []
force_planets = []  #arrays for these values are initilised ere and empty so that in every iteration, the values calculated can be added for the four planets
potenergy_planets = []
radius_orbit = []
timeperiod_planets = []

class Planet(object):

    def __init__(self, name, rad, mass, dist, color, vX, vY):
        """ Initialise the variables """
        self.name = name   #name of the planet
        self.rad = rad     #radius of the planet
        self.mass = mass   #its mass
        self.dist = dist   #its distance by the Sun
        self.color = color #colour its trajectory
        self.x = dist     
        self.y = 0.
        self.vX = vX       #initial velocity along x and y
        self.vY = vY
        

        """ Set the turtle """
        self.mw = Turtle()
        self.mw.color(self.color)
        self.mw.shape("circle")

        """ Set the initial position of the sun """
        self.mw.up()
        self.mw.goto(self.x, self.y)
        self.mw.down()

    """ Getters and setters """
    def getName(self):
        return self.name

    def getRad(self):
        return self.rad

    def getMass(self):
        return self.mass

    def getDist(self):
        return self.dist
    
    def getX(self):
        return self.x

    def getY(self):
        return self.y
    
    def movePos(self, nX, nY):
        self.x = nX
        self.y = nY
        self.mw.goto(nX, nY)

    def getVX(self):
        return self.vX

    def getVY(self):
        return self.vY

    def setXV(self, newVX):
        self.vX = newVX

    def setYV(self, newVY):
        self.vY = newVY
        


class Sun(object):

    """ Initialise the sun class """
    def __init__(self):
        self.name = "Sun"
        self.rad = 100.0
        self.m = 15000.0
        self.c = "yellow"
        self.x = 0
        self.y = 0

        """Set the given position of each planet """
        self.mw = Turtle()
        self.mw.up()
        self.mw.goto(self.x, self.y)
        self.mw.down()
        self.mw.dot(self.rad, self.c)

    def getMass(self):
        return self.m

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def __str__(self):
        return self.name

""" Solar system class """
class milkyWay(Planet):
    
    """ Initialise the solar system, keep track of each planet """
    def __init__(self, sun):
        self.star = None
        self.planets = []
        self.mw = Turtle()
        self.mw.ht()

    """ Add planets to the list, initialise the sun """
    def addPlanet(self, planet):
        self.planets.append(planet)

    def addSun(self, planet):
        self.star = planet

    def rotatePlanets(self):
        #verlet velocity method
        # Loop through each planet 
        for p in self.planets:

            # Move the planet 
            
             #Displacement from initial position
            rx = self.star.getX() - p.getX()
            ry = self.star.getY() - p.getY()

            #Pythagoras theorem, to get distance 
            r = math.sqrt(rx ** 2 + ry ** 2)
            
            

            #Calculate the acceleration of each planet 
            #this is independent on the algorithm, why?
            p.accX = G * self.star.getMass() * rx / r ** 3
            p.accY = G * self.star.getMass() * ry / r ** 3
            
            
           
            # New position following the choosen algorithm
            p.movePos(p.getX() + p.vX*dt + (p.accX*dt**2)/2 , p.getY() + p.vY*dt + (p.accY*dt**2)/2)
            
            
            rx = self.star.getX() - p.getX()
            ry = self.star.getY() - p.getY()

            #Pythagoras theorom, to get new distance
            r = math.sqrt(rx ** 2 + ry ** 2)

            #Calculate the  new acceleration of each planet according to the new position which is calculated before to which the planet moves
            #this is independent on the algorithm, why?
            p.accnewX = G * self.star.getMass() * rx / r ** 3
            p.accnewY = G * self.star.getMass() * ry / r ** 3

            # Calculate new velocity of X and Y be coherent with the algorithm choice
            p.setXV(p.getVX() + ((p.accX+p.accnewX)/2)*dt)
            p.setYV(p.getVY() + ((p.accY+p.accnewY)/2)*dt) 
            
            #formulas for energy, potential and force calculate seperate values according to iteration which are appended to the empty arrays initialised in the beginning
            V = -G * self.star.getMass() / r
            potential_planets.append(V)     
            
            E = -G * self.star.getMass() * p.mass / r
            potenergy_planets.append(E)
            
            F = G * self.star.getMass() * p.mass / r**2
            force_planets.append(F)
            
            radius_orbit.append(r)  #radius also added to plots of V,E,F can be made according to distance
          
            
    '''#Euler method has been commented out
       #Loop through each planet 
       for p in self.planets:

            # Move the planet 
            # New position following the choosen algorithm
            p.movePos(p.getX() + p.vX*dt ,p.getY() + p.vY*dt)
                      
            # Displacement from initial position
            rx = self.star.getX() - p.getX()
            ry = self.star.getY() - p.getY()

            # Pythagoras theorem, to get distance 
            r = math.sqrt(rx ** 2 + ry ** 2)

            # Calculate the acceleration of each planet 
            #this is independent on the algorithm, why?
            accX = G * self.star.getMass() * rx / r ** 3
            accY = G * self.star.getMass() * ry / r ** 3

            # Calculate new velocity of X and Y be coherent with the algorithm choice
            p.setXV(p.getVX() + accX*dt)
            p.setYV(p.getVY() + accY*dt)  
            
            V = -G * self.star.getMass() / r
            potential_planets.append(V)
            
            E = -G * self.star.getMass() * p.mass / r
            potenergy_planets.append(E)
            
            F = G * self.star.getMass() * p.mass / r**2
            force_planets.append(F)
            
            radius_orbit.append(r)'''
            

def initiateGalaxy():

    """ Create instances of each planet and the sun """
    sun = Sun()
    mw = milkyWay(sun)
    mw.addSun(sun)

    p1 = Planet("Planet 1", 19, 10, 220, "green", 0.0, 4.0)
    mw.addPlanet(p1)

    p2 = Planet("Planet 2", 30, 60, 280, "blue", 0.0, 5)
    mw.addPlanet(p2)

    p3 = Planet("Planet 3", 40, 50, 320, "red", 0.0, 6.3)
    mw.addPlanet(p3)

    p4 = Planet("Planet 4", 60, 100, 350, "purple", 0.0, 3.)  #all planets are on same y level but different x levels because of starting initial point kept same
    mw.addPlanet(p4)

    """ Initialise how many cycles the simulation should run for """
    
   
    pi=np.pi
    totaltime=((4*pi**2)*60**3/(G*100))**0.5  
    timeFrame = int((2*totaltime)/dt)
    #time loop choose the total number of steps properly 
    #It may need a few periods 
    print timeFrame
    
    for m in range(timeFrame):
        mw.rotatePlanets()
        """ here calculate the period and semi major, semi minor axis"""  
    
    Pmax1 = max(radius_orbit[0::4]) #maximum val of radius from sun is semi major axis
    Pmin1 = min(radius_orbit[0::4]) #min value of radius from sun is semi minor axis
    print 'semi major axis for green planet 1 is', Pmax1
    print 'semi minor axis is green planet 1 is', Pmin1
    c= 0.5*(Pmax1+Pmin1)
    t = c**1.5
    timeperiod_planets.append(t)  #time period calculated for each planet depending on semimajor axis of orbit'''
    
    Pmax2 = max(radius_orbit[1::4]) #maximum val of radius from sun is semi major axis
    Pmin2 = min(radius_orbit[1::4]) #min value of radius from sun is semi minor axis
    print 'semi major axis for blue planet 2 is', Pmax2
    print 'semi minor axis is blue planet 2 is', Pmin2
    c= 0.5*(Pmax2+Pmin2)
    t = c**1.5
    timeperiod_planets.append(t)  #time period calculated for each planet depending on semimajor axis of orbit'''
    
    Pmax3 = max(radius_orbit[2::4]) #maximum val of radius from sun is semi major axis
    Pmin3 = min(radius_orbit[2::4]) #min value of radius from sun is semi minor axis
    print 'semi major axis for red planet 3 is', Pmax3
    print 'semi minor axis is red planet 3 is', Pmin3
    c= 0.5*(Pmax3+Pmin3)
    t = c**1.5
    timeperiod_planets.append(t)  #time period calculated for each planet depending on semimajor axis of orbit'''
    
    Pmax4 = max(radius_orbit[3::4]) #maximum val of radius from sun is semi major axis
    Pmin4 = min(radius_orbit[3::4]) #min value of radius from sun is semi minor axis
    print 'semi major axis for purple planet 4 is', Pmax4
    print 'semi minor axis is purple planet 4 is', Pmin4
    c= 0.5*(Pmax4+Pmin4)
    t = c**1.5
    timeperiod_planets.append(t)  #time period calculated for each planet depending on semimajor axis of orbit'''
    
       
    #print potential_planets
    #print potenergy_planets
    #print force_planets
    #print radius_orbit
    print 'time periods of planets 1,2,3,4 are', timeperiod_planets
    print len(timeperiod_planets)
    

initiateGalaxy()  #graphs calculated of values against distance, second portion of graphs to see fluctuations are calculated by plotting only en, pot and force values agasint no distance'''

plt.plot(radius_orbit[0::4],potential_planets[0::4], color = 'green', label = "Planet 1")
plt.plot(radius_orbit[1::4],potential_planets[1::4], color = 'blue', label = "Planet 2")
plt.plot(radius_orbit[2::4],potential_planets[2::4], color = 'red', label = "Planet 3")
plt.plot(radius_orbit[3::4],potential_planets[3::4], color = 'purple', label = "Planet 4")
plt.xlabel(' Radius (R) ')
plt.ylabel(' Potential (V) ')
plt.legend()
plt.show()  #the list is divided into the first value for every four groupings. these values 0::4 represent planet 1 while 1::4 planet 2, 2::4 planet 3 and 3::4 planet4

plt.plot(radius_orbit[0::4],potenergy_planets[0::4], color = 'green', label = "Planet 1")
plt.plot(radius_orbit[1::4],potenergy_planets[1::4], color = 'blue', label = "Planet 2")
plt.plot(radius_orbit[2::4],potenergy_planets[2::4], color = 'red', label = "Planet 3")
plt.plot(radius_orbit[3::4],potenergy_planets[3::4], color = 'purple', label = "Planet 4")
plt.xlabel(' Radius (R) ')
plt.ylabel(' Potential Energy (E) ')
plt.legend()
plt.show()
    
plt.plot(radius_orbit[0::4],force_planets[0::4], color = 'green', label = "Planet 1")
plt.plot(radius_orbit[1::4],force_planets[1::4], color = 'blue', label = "Planet 2")
plt.plot(radius_orbit[2::4],force_planets[2::4], color = 'red', label = "Planet 3")
plt.plot(radius_orbit[3::4],force_planets[3::4], color = 'purple', label = "Planet 4")
plt.xlabel(' Radius (R) ')
plt.ylabel(' Force (F) ')
plt.legend()
plt.show()

plt.plot(radius_orbit[0::4], color = 'green', label = "Planet 1")
plt.plot(radius_orbit[1::4], color = 'blue', label = "Planet 2")
plt.plot(radius_orbit[2::4], color = 'red', label = "Planet 3")
plt.plot(radius_orbit[3::4], color = 'purple', label = "Planet 4")
plt.xlabel(' Radius (R) ')
plt.legend()
plt.show()

#the list is seperated into the values calculated for each planet an dthhus different significant plots are made for each planets, energy, potential, force acting on it and radial distance from sun
plt.plot(potential_planets[0::4], color = 'green', label = "Planet 1")
plt.plot(potential_planets[1::4], color = 'blue', label = "Planet 2")
plt.plot(potential_planets[2::4], color = 'red', label = "Planet 3")
plt.plot(potential_planets[3::4], color = 'purple', label = "Planet 4")
plt.ylabel(' Potential (V) ')
plt.legend()
plt.show()



