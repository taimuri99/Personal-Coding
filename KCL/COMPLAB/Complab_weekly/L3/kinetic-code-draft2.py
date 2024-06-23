from numpy import *
from scipy.signal import correlate2d
from numpy.random import randint,choice,uniform
from matplotlib.pyplot import *
from os import system

timepass=[]    #arrays inititialsied globally so the values calculated in loop can be added to them as a list. Helps for the plot later on
sumocc=[]
sumunocc=[]
sumislands=[]
sumunbound=[]

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
        print 'timeincrement is' , dt
        
   
        
        # carry out events
        self.evolve(events_index,dt)
    

class LatticeGas(object) :
	
    def __init__(self, n_sites, density, barriers, prefactors, beta=1 ) :

        # initialisation of lattice gas with given denisty
        self.size = 2*int(ceil(sqrt(n_sites)/2))  #length of square lattice formed
        self.state = choice([0,1],size=(self.size,self.size),p=[1-density,density])
        self.time = 0.
        self.barriers=array(barriers)
        self.prefactors=array(prefactors)
        self.beta = beta
        self.totsites=self.size*self.size
        self.dep = 0
        self.sum_occunbound=0

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
			
			self.dep=self.dep+1 #changes a 0site randomly to 1 to mimic a particle dropping into an unoccupied site
			
			i,j = self.randsite_e()
			self.state[i,j] = 1
			print 'I did a deposition onto site', i,j, 'occ', self.state[i,j]
		
		x = self.state
		
		
		#diffusion
		if events_index == 1 :
			
			i,j = self.randsite_o()  #checks initial conditions satisfied of no occupiued neighbours and then chooses a random movement north,south,east or west
					
			if ((self.state[i,j]==1) and (self.state[i,((j+1)%self.size)]==0) and (self.state[((i+1)%self.size),j]==0) and (self.state[i,j-1]==0) and (self.state[i-1,j]==0)) :  #nsew neighbours checked if unoccupied so can diffuse there
					for x in range(4):
						randpos=random.randint(0,4)
						if randpos ==0 :
							self.state[i,((j+1)%self.size)] = 1   #northneighbourshifted
							self.state[i,j] = 0 
							print 'diffused north'
						if randpos ==1 :
							self.state[((i+1)%self.size),j] = 1  #westneighbourshifted
							self.state[i,j] = 0 
							print 'diffused west'
						if randpos ==2 :
							self.state[i,j-1] = 1  #southneighbourshifted
							self.state[i,j] = 0 
							print 'diffused south'
						if randpos ==3 :
							self.state[i-1,j] = 1  #eastneighbourshifted
							self.state[i,j] = 0 
							print 'diffused east'
							
							
		#firstdiffusionchosen,and then prerequisistes given that all direct neighbours are unoccupied#4444
							
		
		#aggregation(1)
		if events_index == 2 :
			
			i,j = self.randsite_o()  #if initial conditions satisfied, secondary conditions for neighbours' neighbours also satisfied then random movement is chosen using randpos (nsew)
				
			if ((self.state[i,j]==1) and (self.state[i,((j+1)%self.size)]==0) and (self.state[((i+1)%self.size),j]==0) and (self.state[i,j-1]==0) and (self.state[i-1,j]==0)) : #nsew neighbours checked if unoccupied so can diffuse there
					for x in range(4):
						randpos=random.randint(0,3)
						if randpos ==0 : #northneighbour
							if (self.state[i-1,((j+1)%self.size)]==0 and self.state[((i+1)%self.size),((j+1)%self.size)]==0 and self.state[i,((j+2)%self.size)]==0):
								self.state[i,((j+1)%self.size)] = 1   #northneighboursneighbourshift
								self.state[i,j] = 0 
								print 'aggregated north'
								
						if randpos ==1 : # west neighbour
							if (self.state[((i+1)%self.size),((j+1)%self.size)]==0 and self.state[((i+2)%self.size),j]==0 and self.state[((i+1)%self.size),j-1]==0):
								self.state[((i+1)%self.size),j] = 1  #westneighboursneighbourshifted
								self.state[i,j] = 0 
								print 'aggregated west'
								
						if randpos ==2 : #south neighbour
							if (self.state[i-1,j-1]==0 and self.state[((i+1)%self.size),j-1]==0 and self.state[i,j-2]==0):
								self.state[i,j-1] = 1  #southneighboursneighbourshifted
								self.state[i,j] = 0 
								print 'aggregated south'
								
						if randpos ==3 : #east neighbour
							if (self.state[i-1,j-1]==0 and self.state[i-1,((j+1)%self.size)]==0 and self.state[i-2,j]==0):
								self.state[i-1,j] = 1  #eastneighboursneighbourshifted
								self.state[i,j] = 0
								print 'aggregated east'
								
								
		
		#aggregation(2)
		if events_index == 3 :
			
			i,j = self.randsite_o()
			
			if ((self.state[i,j]==1) and (self.state[i,((j+1)%self.size)]==0) and (self.state[((i+1)%self.size),j]==0) and (self.state[i,j-1]==0) and (self.state[i-1,j]==0)) :  #nsew neighbours checked if unoccupied so can diffuse there
					for x in range(4):
						randpos=random.randint(0,3)
						if randpos ==0 : #northneighbour
							if (self.state[i-1,((j+1)%self.size)]==0 and self.state[((i+1)%self.size),((j+1)%self.size)]==0 and self.state[i,((j+2)%self.size)]==0):
								self.state[i,((j+1)%self.size)] = 1   #northneighboursneighbourshift
								self.state[i,j] = 0 
								print 'agglomerated north'
								
						if randpos ==1 : # west neighbour
							if (self.state[((i+1)%self.size),((j+1)%self.size)]==0 and self.state[((i+2)%self.size),j]==0 and self.state[((i+1)%self.size),j-1]==0):
								self.state[((i+1)%self.size),j] = 1  #westneighboursneighbourshift
								self.state[i,j] = 0 
								print 'agglomerated west'
						
						if randpos ==2 : #south neighbour
							if (self.state[i-1,j-1]==0 and self.state[((i+1)%self.size),j-1]==0 and self.state[i,j-2]==0):
								self.state[i,j-1] = 1  #southneighboursneighbourshift
								self.state[i,j] = 0 
								print 'agglomerated south'
								
						if randpos ==3 : #east neighbour
							if (self.state[i-1,j-1]==0 and self.state[i-1,((j+1)%self.size)]==0 and self.state[i-2,j]==0):
								self.state[i-1,j] = 1  #eastneighboursneighbourshifted
								self.state[i,j] = 0
								print 'agglomerated east'
								
		
		
		sum_occ = count_nonzero(self.state)
		print 'Total number of occupied sites is' , sum_occ
		
		sum_unocc = self.totsites-sum_occ
		print 'Total number of unoccupied sites is' , sum_unocc
		
		self.time = self.time+dt
		print 'Total time passes is' , self.time
		
		print 'Number of depositions is' , self.dep
		
		if ((self.state[i,j]==1) and (self.state[i,((j+1)%self.size)]==0) and (self.state[((i+1)%self.size),j]==0) and (self.state[i,j-1]==0) and (self.state[i-1,j]==0)) : 
			self.sum_occunbound=self.sum_occunbound+1
		
		sum_island = self.totsites-(sum_unocc + self.sum_occunbound)
		print 'Total number of islands is' , sum_island
		
		print 'Total number of unbound molecules are' , self.sum_occunbound
		
		# outputted values of islands, total time, occupied and unoccupied sites and unbound sites and then added to empty arrays initialised in the beginning 
		
		timepass.append(self.time)
		sumocc.append(sum_occ)
		sumunocc.append(sum_unocc)
		sumislands.append(sum_island)
		sumunbound.append(self.sum_occunbound)
		
        
        # pick a site that is compatible with the chosen event and increment according to dt
        # not implemented yet
        # dropping of a new particle on the surface
        # if any event takes place it will state that and if no event takes place or it does not meet requirements of ifelse statement, it wont do anythinh just give timeval and prereqs
        
	
    
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

        #figure(figsize=(5,5))
        #imshow(self.state,interpolation="none",cmap='bone_r',vmin=0,vmax=1)
        #axis('off')
        #show()
        
		#above code is hashtagged so the graph could be outputted

p=[0.001,0.1,0.1,0.1]
e=[0.01,-0.35,0.05,0.5]



        
# example of how code should work


model = LatticeGas(50,0.1,barriers=e,prefactors=p)

kmc = KineticMonteCarlo(model)

timestep = kmc.time_step()


for i in range(0,500):
	kmc.time_step()
	model.show()
	print i
print (timepass,sumocc,sumunocc,sumislands,sumunbound)


plot(timepass,sumocc,color='skyblue',label="Occupied")
plot(timepass,sumunocc,color='purple', label="Unoccupied")
plot(timepass,sumislands,color='green', label="Islands")
plot(timepass,sumunbound,color='red', label="Unbound")
xlabel('Time (s)')
ylabel('Number of sites')
legend()
show()


