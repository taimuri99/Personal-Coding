from numpy import *
from numpy.random import randint,choice,uniform
from scipy.signal import correlate2d
from scipy.special import ellipk,ellipe
from matplotlib.pyplot import *

class MetropolisMonteCarlo(object) :

    def __init__( self, Model ) :

        """Initialise a Metropolis Monte Carlo sampling method for
        the given model. The model object requires the following
        vectorised methods:

        Model.evolve : <method>
            returns new states for step propositions
        Model.get_probability : <method>
            returns acceptance probabilities
        """

        
        self.Model = Model
        self.state = Model.state

       
        self.evolve = Model.evolve
        self.get_probability = Model.get_probability

    def step( self ) :
       
        state = self.evolve()

        
        alpha = self.get_probability(state)

        
        acceptance_mask = uniform(0,1,size=alpha.shape) < alpha
        self.state[acceptance_mask] = state[acceptance_mask]

        return True

    def equilibrate( self, n_steps=None ) :
       

        n_steps = n_steps if n_steps is not None else self.state.size
        for _ in xrange(n_steps) :
            self.step()

    def get_samples( self, method_name, args=(), n_samples=None ) :
        

        n_samples = n_samples if n_samples is not None else self.state.size
        return array([ getattr(self.Model,method_name)(*args) for _ in xrange(n_samples) if self.step() ])

class isingModel(object) :

    
    def __init__(self, n_sites, interaction_strength=1.0, beta = 1.0,
                 interaction_kernel=None ) :

      
        self.size = 2*int(ceil(sqrt(n_sites)/2))
        self.interaction_strength = interaction_strength
        self.beta = beta

        if interaction_kernel is None: 
            self.interaction_kernel = array([[0,1,0],[1,0,1],[0,1,0]])
        else : 
            self.interaction_kernel = interaction_kernel

     
        self.state = 2*randint(2,size=(self.size,self.size))-1

       
        self.sublattice = array([[ i+j &1 for i in xrange(self.size) ]
                                   for j in xrange(self.size)], dtype=bool)

    def get_energy(self,temperature=None,interaction_strength=None,dphi=0.01) :
       
        if temperature is None and interaction_strength is None : 
            return -self.interaction_strength/2 * mean(self.field*self.state)

        else : 

            
            beta = 1.0/temperature if temperature is not None else self.beta
            interaction_strength = interaction_strength if interaction_strength is not None else self.interaction_strength
            z = 2*beta*interaction_strength

            kappa = 1/sinh(z)**2
            prefactor = interaction_strength/tanh(z)
            scaling = 2/pi*(2*tanh(z)**2-1)

            
            phi = arange(0,pi/2,dphi)
            integral = array([ sum(1/sqrt( 1-4*k*(sin(phi)/(1+k))**2 ))*dphi for k in kappa])

            return -prefactor*(1+scaling*integral)

    def get_magnetisation(self,temperature=None,interaction_strength=None) :
        
        if temperature is None and interaction_strength is None : 
            return mean(self.state)

        else : 

            
            beta = 1.0/temperature if temperature is not None else self.beta
            interaction_strength = interaction_strength if interaction_strength is not None else self.interaction_strength
            z = 2*beta*interaction_strength

            magnetisation = (1.0-sinh(z)**(-4.0)) **(1.0/8.0)
            return nan_to_num(magnetisation)

    def get_heat(self,temperature=None,interaction_strength=None) :
        
        
        beta = 1.0/temperature if temperature is not None else self.beta
        interaction_strength = interaction_strength if interaction_strength is not None else self.interaction_strength
        z = 2*beta*interaction_strength

        k = 2*sinh(z) / cosh(z)**2
        prefactor = (z/tanh(z))**2 / pi

        scaling = 1-tanh(z)**2
        integral = pi/2+ellipk(k)*(2*tanh(z)**2-1)
        return (ellipk(k)-ellipe(k)-scaling*integral)*prefactor-0.3

    def get_probability(self,state) :
        
        energy_change = -self.interaction_strength * self.field*(state-self.state)
        negatives = energy_change < 0

        energy_change[negatives] = 0
        return exp(-self.beta*energy_change)

    def get_field(self) :
        
        self.field = correlate2d( self.state, self.interaction_kernel, mode='same', boundary='wrap')
        return self.field

    def evolve(self) :
        

        
        self.field = self.get_field()

        
        self.sublattice = logical_not(self.sublattice)

        
        state = array(self.state)
        state[self.sublattice] *= -1
        return state

    def show(self) :
        

        figure(figsize=(5,5))
        imshow(self.state,interpolation="none",cmap='bone',vmin=-1,vmax=1)
        axis('off')


def show_results() :
   

    
    fig, ax1 = subplots(figsize=(10,10))
    ax2 = ax1.twinx()

    
    ax1.plot(temperature_numerical/interaction,
             magnetisation_numerical,
             color='red',marker='.',linestyle='')
    ax2.plot(temperature_numerical/interaction,
             energy_numerical,
             color='black',marker='.',linestyle='')
    ax1.plot(temperature_numerical/interaction,
             heat_numerical,
             color='orange',marker='.',linestyle='')

    
    ax1.plot(temperature_theory/interaction,
             magnetisation_theory,
             color='red',label=r'Magnetisation, $M$')
    ax2.plot(temperature_theory/interaction,
             energy_theory,
             color='black',label=r'Energy per site, $U$')
    ax1.plot(temperature_theory/interaction,
             heat_theory,
             color='orange',label=r'Specific Heat Capacity, $c$')

    
    ax1.set_xlabel(r'Temperature $T$, / $Jk_B^{-1}$',fontsize=16)
    ax1.legend(loc=2,fontsize=16)
    ax2.legend(loc=1,fontsize=16)
    show()


def calculate() :
   

    global temperature_theory,temperature_numerical
    global magnetisation_theory,magnetisation_numerical
    global energy_theory,energy_numerical
    global heat_theory,heat_numerical

   
    betas = 1/(linspace(*temperature_range,num=88)*interaction)
    temperature_theory = linspace(*temperature_range,num=200)*interaction
    temperature_numerical = 1/betas

    
    ising = isingModel(n_spins,beta=betas[0],interaction_strength=interaction)
    MonteCarlo = MetropolisMonteCarlo(ising)
    MonteCarlo.equilibrate()

    
    magnetisation_numerical = array([ MonteCarlo.get_samples('get_magnetisation') for ising.beta in betas ])
    energy_numerical = array([ MonteCarlo.get_samples('get_energy') for ising.beta in betas ])

    magnetisation_theory = ising.get_magnetisation(temperature_theory)
    energy_theory = ising.get_energy(temperature_theory)
    heat_theory = ising.get_heat(temperature_theory)

    magnetisation_numerical = abs(mean(magnetisation_numerical,axis=1))
    heat_numerical = n_spins*var(energy_numerical,axis=1)*betas**2
    energy_numerical = mean(energy_numerical,axis=1)


def main() :

    
    global n_spins
    n_spins = 2**8

    global interaction
    interaction = 2.0

    global temperature_range
    temperature_range = (1.5,3.25)

    
    calculate()

    
    show_results()


if __name__ == '__main__':
    main()
