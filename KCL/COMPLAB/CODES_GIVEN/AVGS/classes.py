#"""""
# A brief summary of classes/object
# Using classes to define objects provides a templating facility: class attributes and methods need only be defined once, and you can then instantiate any number of objects, all sharing the same methods. A class is instantiated by calling the class object. e.g.
#               i = myclassname()
#   i is s an instance of the class myclassname.
# There are a few methods available for each class, which can be used to customise your object. They are
# __init__
# you pass parameters to the class when instantiating it, and the __init__ method will be called with these parameters. Usually the method will set various instance variables via self
# __str__(self)
#  exactly like __repr__ except that it is called when the builtin str function is applied to an object. this function is also called for the %s escape of the % operator. For users
# __del__
# called when an object is deleted
# __repr__(self)
# takes exactly one parameter, self, and must return a string. This string is intended to be a representation of the object. For programmers.
# it is conventional in Python to name the first parameter of a method self
#Classes have five attributes stored in the class's name
#__dict__
#__name__
#If the value of __name__  attribute is set to '__main__'  when module run as main program. Otherwise the value of __name__  is set to contain the name of the module.
#__bases__
#__doc__
#__module__
# """"
import numpy as np
# Exampel of classes 
# Projectile's motion: height and distance for a projectile with
# given initial velocity and angle

class projectile:

    # this is called every time a new object is created
    def __init__ (self, v=1.0, theta=45, grav=9.81):

        self.v = v           # velocity m/s
        self.theta = theta   # angle (degrees)
        
        self.thetaRad = theta*np.pi/180.0
        self.vx = v*np.cos(self.thetaRad)
        self.vy = v*np.sin(self.thetaRad)

        self.g = grav


    def maxheight(self):

        # how high does this projectile go?
        # vf_y^2 = 0 = vi_y^2 - 2 g h
        h = self.vy**2/(2.0*self.g)

        return h

    def distance(self):
        
        # time of flight up
        # vf_y = 0 = vi_y - g t
        t = self.vy/self.g

        # total time = up + down
        t = 2.0*t

        d = self.vx*t

        return d

    def __str__(self):
        # a string representation for this class -- so we can print it
        str = " v: %s m/s\n theta: %s (degrees)\n Max-height: %s m\n Max-distance: %s m\n" % \
            (`self.v`, `self.theta`, `self.maxheight()`, `self.distance()`)
        
        return str

# fire a projectile
p1 = projectile()
print p1

# create a list of projectiles
projectiles = []
projectiles.append(p1)

print "For p1"
print projectiles

print "For p2"
print "Max height: in m\n"
p2=projectile(100, 80)
print p2
projectiles.append(p2)
print projectiles

#projectiles.append(projectile(v = 1000, theta = 30))
#print p.maxheight()
     
### "Varying v"

for v in np.arange(1.0,100.,20):
    projectiles.append(projectile(v, 45))

### print "Varying theta"        
for ang in np.arange(10, 90, 10):
    vel = 10.
    projectiles.append(projectile(vel, ang))

print "Vel., Angle, Max height: in m\n, Max distance in m\n"
for p in projectiles:
    print p.v, p.theta, p.maxheight(), p.distance()



        
