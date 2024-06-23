# importing libraries
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

        # reference state
        self.Model = Model
        self.state = Model.state

        # acceptance probability and evolution functions
        self.evolve = Model.evolve
        self.get_probability = Model.get_probability

    def step( self ) :
        """Perform parallel Metropolis-Hastings update on model

        ---returns---
        flag : <bool>
            Flag that states step has been successful
        """

        # sample proposed state
        state = self.evolve()

        # calculate acceptance probability
        alpha = self.get_probability(state)

        # accept jumps with probability alpha
        acceptance_mask = uniform(0,1,size=alpha.shape) < alpha
        self.state[acceptance_mask] = state[acceptance_mask]

        return True

    def equilibrate( self, n_steps=None ) :
        """Equilibrate model within the default estimate of self.state.size iterations

        ---parameters---
        n_steps : <int>
            Number of steps to take for equilibration,
            if None then defaults to number of state variables
        """

        n_steps = n_steps if n_steps is not None else self.state.size
        for _ in xrange(n_steps) :
            self.step()

    def get_samples( self, method_name, args=(), n_samples=None ) :
        """Sample model attributes given current equilibrated state

        ---parameters---

        method_name : <str>
            Name of methods to call to retrieve sample from self.Model

        args : <tuple>
            Additional arguments for method calls

        n_samples : <int>
            Number of samples to take for attribute,
            if None then defaults to number of state variables

        ---returns---
        samples : <1d-array>
            Returns array with n_samples of method
        """

        n_samples = n_samples if n_samples is not None else self.state.size
        return array([ getattr(self.Model,method_name)(*args) for _ in xrange(n_samples) if self.step() ])

class isingModel(object) :

    '''Two dimensional Ising model, with spins (-1,1) evolving
    according to any given interaction kernel; defaults to
    nearest-neighbour, with periodic boundary conditions

    ---parameters---

    n_sites : <int> or <float>
        Number of sites that approximates to nearest N**2

    interaction_kernel : <2d-array>
        Interaction kernel when convolved with the values
        of the local spins gives the local field around
        the central spin

    '''
    def __init__(self, n_sites, interaction_strength=1.0, beta = 1.0,
                 interaction_kernel=None ) :

        # size, interaction strength and inverse temperature
        self.size = 2*int(ceil(sqrt(n_sites)/2))
        self.interaction_strength = interaction_strength
        self.beta = beta

        if interaction_kernel is None: # defaults to nearest-neighbour
            self.interaction_kernel = array([[0,1,0],[1,0,1],[0,1,0]],dtype=int8)
        else : # custom interaction rule
            self.interaction_kernel = interaction_kernel

        # random initialisation of spins (-1,1)
        self.state = 2*randint(2,size=(self.size,self.size),dtype=int8)-1

        # chessboard sublattice
        self.sublattice = array([[ i+j &1 for i in xrange(self.size) ]
                                   for j in xrange(self.size)], dtype=bool)

    def get_energy(self,temperature=None,interaction_strength=None,dphi=0.01) :
        '''Calculates the mean energy per site of the lattice

        ---parameters---

        temperature : <ndarray>
            return theoretical result at specific temperatures,
            if None then results numerical result

        interaction_strength : <ndarray>
            return theoretical result at specific interaction strengths,
            if None then results numerical result

        ---returns---
        energy : <float>
            mean energy per site
        '''
        if temperature is None and interaction_strength is None : # numerical result
            return -self.interaction_strength/2 * mean(self.field*self.state)

        else : # otherwise return onsagers exact solution

            # optionally provided arguments
            beta = 1.0/temperature if temperature is not None else self.beta
            interaction_strength = interaction_strength if interaction_strength is not None else self.interaction_strength
            z = 2*beta*interaction_strength

            kappa = 1/sinh(z)**2
            prefactor = interaction_strength/tanh(z)
            scaling = 2/pi*(2*tanh(z)**2-1)

            # perform phi integral
            phi = arange(0,pi/2,dphi)
            integral = array([ sum(1/sqrt( 1-4*k*(sin(phi)/(1+k))**2 ))*dphi for k in kappa])

            return -prefactor*(1+scaling*integral)

    def get_magnetisation(self,temperature=None,interaction_strength=None) :
        '''Calculates the magnetisation of the lattice

        ---parameters---

        temperature : <ndarray>
            return theoretical result at specific temperatures,
            if None then results numerical result

        interaction_strength : <ndarray>
            return theoretical result at specific interaction strengths,
            if None then results numerical result

        ---returns---
        magnetisation : <float> or <ndarray>
            magnetisation
        '''
        if temperature is None and interaction_strength is None : # numerical result
            return mean(self.state)

        else : # otherwise return onsagers exact solution

            # optionally provided arguments
            beta = 1.0/temperature if temperature is not None else self.beta
            interaction_strength = interaction_strength if interaction_strength is not None else self.interaction_strength
            z = 2*beta*interaction_strength

            magnetisation = (1.0-sinh(z)**(-4.0)) **(1.0/8.0)
            return nan_to_num(magnetisation,copy=False)

    def get_heat(self,temperature=None,interaction_strength=None) :
        '''Calculates the theoretical specific heat of the lattice

        ---paramters---

        temperature : <ndarray>
            return theoretical result at specific temperatures

        interaction_strength : <ndarray>
            return theoretical result at specific interaction strengths

        ---returns---
        specific heat : <float> or <ndarray>
            specific heat
        '''

        # return onsagers exact solution
        beta = 1.0/temperature if temperature is not None else self.beta
        interaction_strength = interaction_strength if interaction_strength is not None else self.interaction_strength
        z = 2*beta*interaction_strength

        k = 2*sinh(z) / cosh(z)**2
        prefactor = (z/tanh(z))**2 / pi

        scaling = 1-tanh(z)**2
        integral = pi/2+ellipk(k)*(2*tanh(z)**2-1)
        return (ellipk(k)-ellipe(k)-scaling*integral)*prefactor-0.3

    def get_probability(self,state) :
        '''Calculates acceptance probability given current configuration

        ---paramters---
        state : <2d-array>
            lattice of shape (self.size,self.size)
            to calculate the acceptance probability of

        ---returns---
        ratios : <2d-array>
            site probabilities of shape (self.size,self.size)
        '''
        energy_change = -self.interaction_strength * self.field*(state-self.state)
        negatives = energy_change < 0

        energy_change[negatives] = 0
        return exp(-self.beta*energy_change)

    def get_field(self) :
        '''Calculates the local field on each site of the lattice

        ---returns---
        field : <2d-array>
            local field on each site of shape (self.size,self.size)
        '''
        self.field = correlate2d( self.state, self.interaction_kernel, mode='same', boundary='wrap')
        return self.field

    def evolve(self) :
        '''Returns lattice states where the currently selected
        sublattice spins are flipped. Each time this method is
        called the selections changes to the sublattice complement

        ---returns---
        states : <2d-array>
            states with flipped sublattice
        '''

        # calculate local field
        self.field = self.get_field()

        # alternate sublattice
        self.sublattice = logical_not(self.sublattice)

        # return flipped sublattice
        state = array(self.state)
        state[self.sublattice] *= -1
        return state

    def show(self) :
        '''Displays current lattice configuration, with black
        and white representing spin -1 and 1 respectively
        '''

        figure(figsize=(5,5))
        imshow(self.state,interpolation="none",cmap='bone',vmin=-1,vmax=1)
        axis('off')


def show_results() :
    '''produces main figure with results'''

    # dual axis mode
    fig, ax1 = subplots(figsize=(10,10))
    ax2 = ax1.twinx()

    # numerical results
    ax1.plot(temperature_numerical/interaction,
             magnetisation_numerical,
             color='red',marker='.',linestyle='')
    ax2.plot(temperature_numerical/interaction,
             energy_numerical,
             color='black',marker='.',linestyle='')
    ax1.plot(temperature_numerical/interaction,
             heat_numerical,
             color='orange',marker='.',linestyle='')

    # theoretical results
    ax1.plot(temperature_theory/interaction,
             magnetisation_theory,
             color='red',label=r'Magnetisation, $M$')
    ax2.plot(temperature_theory/interaction,
             energy_theory,
             color='black',label=r'Energy per site, $U$')
    ax1.plot(temperature_theory/interaction,
             heat_theory,
             color='orange',label=r'Specific Heat Capacity, $c$')

    # plotting options
    ax1.set_xlabel(r'Temperature $T$, / $Jk_B^{-1}$',fontsize=16)
    ax1.legend(loc=2,fontsize=16)
    ax2.legend(loc=1,fontsize=16)
    show()


def calculate() :
    '''performs all calculations'''

    global temperature_theory,temperature_numerical
    global magnetisation_theory,magnetisation_numerical
    global energy_theory,energy_numerical
    global heat_theory,heat_numerical

    # inverse temperature grid
    betas = 1/(linspace(*temperature_range,num=88)*interaction)
    temperature_theory = linspace(*temperature_range,num=200)*interaction
    temperature_numerical = 1/betas

    # setup and equilibrate model
    ising = isingModel(n_spins,beta=betas[0],interaction_strength=interaction)
    MonteCarlo = MetropolisMonteCarlo(ising)
    MonteCarlo.equilibrate()

    # sample from equilibrium ensemble
    magnetisation_numerical = array([ MonteCarlo.get_samples('get_magnetisation') for ising.beta in betas ])
    energy_numerical = array([ MonteCarlo.get_samples('get_energy') for ising.beta in betas ])

    magnetisation_theory = ising.get_magnetisation(temperature_theory)
    energy_theory = ising.get_energy(temperature_theory)
    heat_theory = ising.get_heat(temperature_theory)

    magnetisation_numerical = abs(mean(magnetisation_numerical,axis=1))
    heat_numerical = n_spins*var(energy_numerical,axis=1)*betas**2
    energy_numerical = mean(energy_numerical,axis=1)


def main() :

    # parameters
    global n_spins
    n_spins = 2**8

    global interaction
    interaction = 2.0

    global temperature_range
    temperature_range = (1.5,3.25)

    # run monte carlo routine
    calculate()

    # plot results
    show_results()


# execute main code
if __name__ == '__main__':
    main()
