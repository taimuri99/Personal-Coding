
# coding: utf-8

# In[1]:


from numpy import *
###from scipy.signal import correlate2d
from numpy.random import randint,choice,uniform
from matplotlib.pyplot import *

class KineticMonteCarlo(object) :

    def __init__( self, Model ) :

        """Initialise a Kinetic Monte Carlo evolution method for
        the given model. The model object requires the following
        vectorised methods:

        Model.evolve : <method>
            evolves model according to given events
        Model.get_rates : <method>
            returns event rates
        """

        # reference state
        self.Model = Model
        self.state = Model.state

        # event rates and evolution rule
        self.evolve = Model.evolve
        self.get_rates = Model.get_rates


    def time_step( self ) :
        """Perform parallel Kinetic Monte-Carlo update on model"""

        # calculate transition rate fractions
        rates = self.get_rates()  ### a k-dimensional array 
        probabilities = cumsum(rates,axis=0)  #### a k-dim array

        total_rate = probabilities[-1]
        probabilities /= total_rate

        # choose events according to rates
        u = uniform(0,1)
        events_index = sum(probabilities < u)
        
        # generate waiting time
        v = uniform(0,1)
        dt = -log(v)/total_rate
        
        # carry out events
        self.evolve(events_index,dt)
    

class LatticeGas(object) :

    def __init__(self, n_sites, density ) :

        # initialisation of lattice gas with given denisty
        self.size = 2*int(ceil(sqrt(n_sites)/2))
        self.state = choice([0,1],size=(self.size,self.size),p=[1-density,density])
        self.time = 0.

        ''' this returna a matrix of 0 and 1; 1 for occupied sites'''
        
    def get_rates(self) :
        pass
         #### to return an array that contains the events_rate
        return array([0.5,0.1,1.])
    

    def evolve(self,events_index,dt) :
        pass
        # pick a site that is compatible with the chosen event and increment according to dt
        # not implemented yet

    def show(self) :
        '''Displays current lattice configuration, with black
        and white representing whether there is a partilce or not
        '''

        figure(figsize=(5,5))
        imshow(self.state,interpolation="none",cmap='bone_r',vmin=0,vmax=1)
        axis('off')
        show()

        
# example of how code should work
model = LatticeGas(50,0.1)
kmc = KineticMonteCarlo(model)

x = kmc.time_step()
model.show()

