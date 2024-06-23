import math
from turtle import *

G =  # chose a value between 1-10
dt = # work hard to get this right

class Planet(object):

    def __init__(self, name, rad, mass, dist, color, vX, vY):
        """ Initialise the variables """
        self.name = name   #name of the planet
        self.rad = rad	   #radius of the planet
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
        self.rad = 150.0
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

        """ Loop through each planet """
        for p in self.planets:

            """ Move the planet """
            """ New position following the choosen algorithm"""
            p.movePos(p.getX() + #choose the algorithm..........,
                      p.getY() + #choose the algorithm..........)
                      
            """ Displacement from initial position"""
            rx = self.star.getX() - p.getX()
            ry = self.star.getY() - p.getY()

            """ Pythagoras theorom, to get distance """
            r = math.sqrt(rx ** 2 + ry ** 2)

            """ Calculate the acceleration of each planet """
            #this is independent on the algorithm, why?
            accX = G * self.star.getMass() * rx / r ** 3
            accY = G * self.star.getMass() * ry / r ** 3

            """ Calculate new velocity of X and Y be coherent with the algorithm choice"""
            p.setXV(p.getVX() + #choose the algorithm..........)
            p.setYV(p.getVY() + #choose the algorithm..........)
            

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

    p4 = Planet("Planet 4", 60, 100, 350, "purple", 0.0, 3.)
    mw.addPlanet(p4)

    """ Initialise how many cycles the simulation should run for """
    timeFrame = ..... 
    #time loop choose the total number of steps properly 
    #It may need a few periods 

    for m in range(timeFrame):
        mw.rotatePlanets()
         """ here calculate the period and other output quantities"""  
initiateGalaxy()
    

    

