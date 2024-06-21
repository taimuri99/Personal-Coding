# coding: utf-8

# In[1]:


from numpy import *
from scipy.signal import correlate2d
from numpy.random import randint,choice,uniform
from matplotlib.pyplot import *
from os import system



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
        print self.state

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

    def __init__(self, n_sites, density, barriers, prefactors, beta=1 ) :

        # initialisation of lattice gas with given denisty
        self.size = 2*int(ceil(sqrt(n_sites)/2))
        self.state = choice([0,1],size=(self.size,self.size),p=[1-density,density])
        self.time = 0.
        self.barriers=array(barriers)
        self.prefactors=array(prefactors)
        self.beta = beta
        

        ''' this return a a matrix of 0 and 1; 1 for occupied sites'''
        
    def get_rates(self) :
        
         #### to return an array that contains the events_rate
        
        self.rate=self.prefactors*exp(-self.barriers*self.beta) #claculatesratesandchoosewhicheventtotakeplace
        print 'In get_rate', self.rate
        return self.rate

    def evolve(self,events_index,dt) :
		print('In evolve'), events_index
		
		
		
		#adsorption/depositioncode
		if events_index == 0 :
			i,j = self.randsite_e()
			self.state[i,j] = 1
			print 'I did a deposition onto site', i,j, 'occ', self.state[i,j]
			
		#diffusion	
		if events_index == 1 :
			i,j = self.randsite_o()
			
			if ((self.state[i,j]==1) and (self.state[i,j+1]==0) and (self.state[i+1,j]==0) and (self.state[i,j-1]==0) and (self.state[i-1,j]==0)) :  #nsew neighbours checked if unoccupied so can diffuse there
					for x in range(4):
						randpos=random.randint(0,3)
						if randpos ==0 :
							self.state[i,j+1] = 1   #northneighbouroccupied
							self.state[i,j] = 0 
							print 'north neighbour is occupied'
						if randpos ==1 :
							self.state[i+1,j] = 1  #westneighbourocc
							self.state[i,j] = 0 
							print 'west neighbour is occupied'
						if randpos ==2 :
							self.state[i,j-1] = 1  #southneighbourocc
							self.state[i,j] = 0 
							print 'south neighbour is occupied'
						if randpos ==3 :
							self.state[i-1,j] = 1  #eastneighbourocc
							self.state[i,j] = 0 
							print 'east neighbour is occupied'
							
			print 'Particle diffused to unoccupied site'
			
							
		#firstdiffusionchosen,and then prerequisistes given that all direct neighbours are unoccupied#4444
							
		
		
		if events_index == 2 :
			i,j = self.randsite_o()
			
			if ((self.state[i,j]==1) and (self.state[i,j+1]==0) and (self.state[i+1,j]==0) and (self.state[i,j-1]==0) and (self.state[i-1,j]==0)) : #nsew neighbours checked if unoccupied so can diffuse there
					for x in range(4):
						randpos=random.randint(0,3)
						if randpos ==0 : #northneighbour
							if (self.state[i-1,j+1]==1 and self.state[i+1,j+1]==1 and self.state[i,j+2]==1):
								self.state[i,j+1] = 1   #northneighboursneighbouroccupied
								self.state[i,j] = 0 
								print 'north neighbour has neighbours'
								
						if randpos ==1 : # west neighbour
							if (self.state[i+1,j+1]==1 and self.state[i+2,j]==1 and self.state[i+1,j-1]==1):
								self.state[i+1,j] = 1  #westneighboursneighbourocc
								self.state[i,j] = 0 
								print 'west neighbour has neighbours'
								
						if randpos ==2 : #south neighbour
							if (self.state[i-1,j-1]==1 and self.state[i+1,j-1]==1 and self.state[i,j-2]==1):
								self.state[i,j-1] = 1  #southneighboursneighbourocc
								self.state[i,j] = 0 
								print 'south neighbour has neighbours'
								
						if randpos ==3 : #east neighbour
							if (self.state[i-1,j-1]==1 and self.state[i-1,j+1]==1 and self.state[i-2,j]==1):
								self.state[i-1,j] = 1  #eastneighboursneighbourocc
								self.state[i,j] = 0
								print 'east neighbour has neighbours'
								
			print 'Particle agglomerated to unoccupied site with unoccupied neighbours'
			
		
		if events_index == 3 :
			i,j = self.randsite_o()
			
			if ((self.state[i,j]==1) and (self.state[i,j+1]==0) and (self.state[i+1,j]==0) and (self.state[i,j-1]==0) and (self.state[i-1,j]==0)) :  #nsew neighbours checked if unoccupied so can diffuse there
					for x in range(4):
						randpos=random.randint(0,3)
						if randpos ==0 : #northneighbour
							if (self.state[i-1,j+1]==1 and self.state[i+1,j+1]==1 and self.state[i,j+2]==1):
								self.state[i,j+1] = 1   #northneighboursneighbouroccupied
								self.state[i,j] = 0 
								print 'north neighbour has neighbours'
								
						if randpos ==1 : # west neighbour
							if (self.state[i+1,j+1]==1 and self.state[i+2,j]==1 and self.state[i+1,j-1]==1):
								self.state[i+1,j] = 1  #westneighboursneighbourocc
								self.state[i,j] = 0 
								print 'west neighbour has neighbours'
						
						if randpos ==2 : #south neighbour
							if (self.state[i-1,j-1]==1 and self.state[i+1,j-1]==1 and self.state[i,j-2]==1):
								self.state[i,j-1] = 1  #southneighboursneighbourocc
								self.state[i,j] = 0 
								print 'south neighbour has neighbours'
								
						if randpos ==3 : #east neighbour
							if (self.state[i-1,j-1]==1 and self.state[i-1,j+1]==1 and self.state[i-2,j]==1):
								self.state[i-1,j] = 1  #eastneighboursneighbourocc
								self.state[i,j] = 0
								print 'east neighbour has neighbours'
								
			print 'Particle aggregated to unoccupied site with unoccupied neighbours'
			
		
		
        
        # pick a site that is compatible with the chosen event and increment according to dt
        # not implemented yet
        # dropping of a new particle on the surface
        
		
		
    def randsite_e(self):
		'''here randomly selected an empty site for deposition/adsorption'''
		mask=self.state == 0
		where(mask)
		ilist,jlist=where(mask)
		emptylist=zip(ilist,jlist)
		k=randint(len(emptylist))
		return emptylist[k]       
        
    def randsite_o(self):
		'''here randomly selected an occ- site for diffusuion'''
		mask=self.state==1
		where(mask)
		ilist,jlist=where(mask)
		occlist=zip(ilist,jlist)
		k=randint(len(occlist))
		return occlist[k]

    def show(self) :
        '''Displays current lattice configuration, with black
        and white representing whether there is a partilce or not
        '''

        figure(figsize=(5,5))
        imshow(self.state,interpolation="none",cmap='bone_r',vmin=0,vmax=1)
        axis('off')
        show()

p=[0.001,0.1,0.1,0.1]
e=[0.01,-0.35,0.05,0.5]
t=0
        
# example of how code should work
model = LatticeGas(50,0.1,barriers=e,prefactors=p)
kmc = KineticMonteCarlo(model)

x = kmc.time_step()


for i in range(4):
	kmc.time_step()
	model.show()

system('convert -delay 4 -loop 0 *.png animation.gif')
system('rm *.png')
