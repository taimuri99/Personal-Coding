#=============================================================================
#       Quantum Mechanical Simulation using 
# 		Finite-Difference Time-Domain (FDTD) Method
#	    Report the FDTD method and ensure that there is a strong and clear correspondence between
#	    the theory and the implementation (explicilty report the physical meaning of variables in the code)
#       Say clearly what the program does and what physics want to model.
#       You should describe the code in your report and comment all the relevant part of the code 
#       The time step, dt, has strict requirements or else the simulation becomes unstable.
#       Comment your choice in your report.
#
#       You should tweak the height of the potential as well as the
#       barrier thickness, and report the different behaviour:
#       full reflection with no noticeable transmission, transmission and reflection, mostly transmission with tunneling
#
#       This script requires pylab and numpy. You might consider to change plotting instruction in 
#	matplotlib or turtle if you find it more conformtable. Add a statement on the graphics used
#
#	Please also state here which editor you have used and which if you opt for a Python2.7 or Python3
#============================================================================
# 
# Please enter here your name: Taimur Ahmad Khan
#
#============================================================================
#  
import numpy as np
import pylab
#
# Set pylab to interactive mode so plots update when run outside ipython
pylab.ion()
#
# On your report add a sentence commenting various interactive graphical mode you have used during the module 
#=============================================================================
# Utility functions
#
def Gaussian(x,t,sigma):
    """  A Gaussian curve.
        x = Variable
        t = time shift
        sigma = standard deviation      """
    return np.exp(-(x-t)**2/(2*sigma**2)) #formula to mimic gaussian distribution
def free(npts):
    """This utility function creates a list of zeros of length npts which is defined below. 
    To be used for a free particle and in all the region where the particle is free of any potential obstacles.
    """
    return np.zeros(npts)
def step(npts,v0):
    "Add a function to mimick a step potential of height v0. Return the potential v"
    v = free(npts) 
    v[600:] =v0 #the 'step' will have a potential of V0 which can be changed as to the requirements of the simulations.
    return v   #600 spatial points for the potential step/barrier
def doublebarrier(npts,v0,thickness):
    "Add a function to mimick a barrier of height v0 and thickness thickness"
    v = free(npts)
    v[600:600+thickness]= v0 #first barrier 
    v[900:900+thickness]= v0 #second barrier 
    return v  #both barriers have a thickness of 25 and span from gridpoints 600-625 and 900-925 respectively
def well(npts,v0,thickness):
    v = free(npts)
    v[0:150+thickness] = v0 #left edge
    v[875:] = v0 #right edge
    return v  #both of the wells have set width between the grid points
def doublewell(npts,v0,thickness):
    v=free(npts)
    v[0:thickness] =v0 #left pillar
    v[600:600+thickness] =v0 #middle pillar
    v[1225:1225+thickness] =v0 #right pilllar
    return v #double well shown by three pillars with two pillars showing a well potential
    

#>>>>>>> External Changes
def fillax(x,y,*args,**kw):
    """Fill the space between an array of y values and the x axis.
    All args/kwargs are passed to the pylab.fill function.
    Returns the value of the pylab.fill() call.
    REPORT: In the report section on graphical tools, describe in your report what pylab.fill does
    """
    xx = np.concatenate((x,np.array([x[-1],x[0]],x.dtype)))
    yy = np.concatenate((y,np.zeros(2,y.dtype)))
    return pylab.fill(xx, yy, *args,**kw)
#
#=============================================================================
#
#  Simulation Parameters.  Be sure to include decimal points on appropriate
#  variables so they become floats instead of integers.
#  Comment all the meaning of input variables. 
#  EXPLAIN IN THE CODE AND IN THE REPORT ALL THE CASES YOU CONSIDER:
#  for example 'free: Zero potential, packet propagates freely
#  Add relevant logical variables for all the other cases you must implement
#  Potential parameters: what they mean and what physical process they tune
#  REPORT IN THE CODE AND IN THE REPORT WHAT INPUT PARAMETER YOU NEED AND PASS THEM FROM INPUT
#  ON YOUR REPORT: say clearly which data you collect per each set of input parameters. A table coudl be useful
#

N    =   1250       #number of gridpoints, note that xe^y is x*10^y
T    =   5*N     #  Make sure that the number of total steps avoid to reach the boundary.
Tp   =   50    #  
dx   = 1.0e0    #  
m    = 1.0e0    #  all variables initialised to be used in formulae in functions
hbar = 1.0e0    #  all reduced values in terms of angstroms(position) and hartrees(energy)
#
# 
V0   =  1.5e-02 #the value is changed to 1.0e-02 so that the energy of wave-packet is more than the potential obstacle
THCK =    25    #thickness
POTENTIAL =   'free' #choose the type of potential you want from the functions defined by typing in the chosen potential type      
sigma =   40              

X    = dx*np.linspace(0,N,N)        
x0 = round(N/2) - 5*sigma   
k0 = 2.*np.pi/sigma         #wave-number
E = (hbar**2/2.0/m)*(k0**2+0.5/sigma**2) #equation given
#=============================================================================
# Code begins
if POTENTIAL=='free': #for the free particle with no potential obstacles
    V = free(N)
#
# ADD HERE ALL THE OTHER TYPE OF POTENTIAL YOU WILL IMPLEMENT in agreement with utitlity function deintion
# please make sure that on your report you sketch those potentials
#
elif POTENTIAL=='step':  # for a singular potential region where V of the region is greater than the V in the rest
	V = step(N,V0)
elif POTENTIAL=='doublebarrier':  # for a double barrier potential or 2 potential obstacles
	V= doublebarrier(N,V0,THCK) 
elif POTENTIAL=='well':  # a potential well
	V = well(N,V0,THCK)
elif POTENTIAL == 'doublewell': #a potential doublewell or two potential wells
	V = doublewell(N,V0,THCK)
else:
    raise ValueError("Unrecognized potential type: %s" % POTENTIAL) #function chosen will be executed only and when chosen an output
#of potential, energy of wavepacket, potential obstacle width
#  EXPLAIN IN THE REPORT CLEARLY THE ALGORITHM IMPLEMENTATION. Check the correspondence between
#  equations and variable names in the code. The maximum stable time step is a function of the potential, V
#  COMMENT ON THE choice of time step, true everytime you have a dynamics and the time is discretised
#
Vmax = V.max()
dt   = hbar/(2*hbar**2/(m*dx**2)+Vmax)         #time-step calculated using the equation, maximum value allowed 
c1   = hbar*dt/(m*dx**2)                       # c1 and c2 are the two coefficents in the equations for the predicted future states for both real 
c2   = 2*dt/hbar                               # and imaginary wavefunctions. The equations given below and the report
c2V  = c2*V                                    # CE of the FDTD algorithm which discretises the PDEs
# formula diff from ref paper but introduction of 2 in numerator of both c1 and c2 conserves final simulation
#
#
print 'One-dimensional Schrodinger equation - time evolution'  #the values are calculated and outputted with the simulation of different potential obstacles.
print 'Wavepacket energy:   ',E
print 'Potential type:      ',POTENTIAL
print 'Potential height V0: ',V0
print 'Barrier thickness:   ',THCK
print 'Well thickness:      ',THCK  
#
#============ Algorithm implementation: to be clearly explained in the report
#
psi_r = np.zeros((3,N))     #  Real component
psi_i = np.zeros((3,N))     #  Imaginary component
psi_p = np.zeros(N,)        # Observable probability (magnitude-squared of the complex wave function).
#  Temporal indexing constants, used for accessing rows of the wavefunctions.
PA = 0                 #  Past
PR = 1                 #  Present
FU = 2                 #  Future
#  Initialize wave function.  
xn = range(1,N/2)
x = X[xn]/dx    #  Normalized position coordinate
gg = Gaussian(x,x0,sigma)
cx = np.cos(k0*x)
sx = np.sin(k0*x) 
psi_r[PR,xn] = cx*gg #imaginary and real functions in terms of the wavefunction of free particle, cosx and sin x for real and im respectively
psi_i[PR,xn] = sx*gg
psi_r[PA,xn] = cx*gg
psi_i[PA,xn] = sx*gg
#
#
psi_p = psi_r[PR]**2 + psi_i[PR]**2 #square of mag of wavefunc
#
#  Normalization: what you should do?
#  REPORT THE PROBABILITY ID CALCULATED
P   = dx * psi_p.sum()                      
nrm = np.sqrt(P)
psi_r /= nrm
psi_i /= nrm
psi_p /= P
#
#================================================================
#
#===Plotting Initialize the figure and axes.
pylab.figure()
xmin = X.min()
xmax = X.max()
ymax = 1.5*(psi_r[PR]).max()
pylab.axis([xmin,xmax,-ymax,ymax])
#  Initialize the plots with their own line objects.  The figures plot MUCH
#  faster if you simply update the lines as opposed to redrawing the entire
#  figure.  For reference, include the potential function as well.
lineR, = pylab.plot(X,psi_r[PR],'b',alpha=0.7,label='Real')
lineI, = pylab.plot(X,psi_i[PR],'r',alpha=0.7,label='Imaginary')
lineP, = pylab.plot(X,6*psi_p,'k',label='Probability')
pylab.title('Potential height: %.2e' % V0)
pylab.xlabel('Spatial Points')  #xdomain is position
pylab.ylabel('Wave Function (hartree)') #yrange is energy ibn hartrees
#
#
if Vmax !=0 :
    # Scaling factor for energies, so they fit in the same plot as the wavefunctions
    Efac = ymax/2.0/Vmax
    print 'Factor to printing energy/wfc on the same plot',Efac
    V_plot = V*Efac
#    pylab.plot(X,V_plot,':k',zorder=0)   #  Potential line.
    fillax(X,V_plot, facecolor='y', alpha=0.2,zorder=0)
    # Plot the wavefunction energy, in the same scale as the potential
    pylab.axhline(E*Efac,color='g',label='Energy',zorder=1)
pylab.legend(loc='lower right')
pylab.draw()
#  pylab.xlim(xmin,xmax)
#
#  Direct index assignment is MUCH faster than using a spatial FOR loop, so
#  these constants are used in the update equations.  Remember that Python uses
#  zero-based indexing.
#
IDX1 = range(1,N-1)                            #  psi [ k ]
IDX2 = range(2,N)                              #  psi [ k + 1 ]
IDX3 = range(0,N-2)                            #  psi [ k - 1 ]
#the before and after points of differential steps 
#=============================================================Time evolution
#
for t in range(T+1):
    # Precompute a couple of indexing constants, this speeds up the computation
    psi_rPR = psi_r[PR]
    psi_iPR = psi_i[PR]
    #  Apply the update equations calculated using the fDTD method, calculating predicted future states
    psi_i[FU,IDX1] = psi_i[PA,IDX1] + \
                      c1*(psi_rPR[IDX2] - 2*psi_rPR[IDX1] +
                          psi_rPR[IDX3])
    psi_i[FU] -= c2V*psi_r[PR]

    psi_r[FU,IDX1] = psi_r[PA,IDX1] - \
                      c1*(psi_iPR[IDX2] - 2*psi_iPR[IDX1] +
                          psi_iPR[IDX3])
    psi_r[FU] += c2V*psi_i[PR]
    #  Increment the time steps.  PR -> PA and FU -> PR as present becomes past, future becomes present, new future is calculated and process is repeated and then next future is calculated
    psi_r[PA] = psi_rPR
    psi_r[PR] = psi_r[FU]
    psi_i[PA] = psi_iPR
    psi_i[PR] = psi_i[FU]
    #
     #
    if t % Tp == 0: #keeps t in range of timeperiod and npot exceed it. if t reaches desired amount, the simulation stops
        #  Compute observable probability for the plot
        #  REPORT: Explain what you ar plotting and what is the probability
        #  This is mandatory to understand what data you collect
        psi_p = psi_r[PR]**2 + psi_i[PR]**2 #probability
        #  Update the plots.
        lineR.set_ydata(psi_r[PR])
        lineI.set_ydata(psi_i[PR])
        # Note: we plot the probability density amplified by a factor so it's a bit easier to see.
        lineP.set_ydata(6*psi_p)
#
        pylab.draw() #moving simulation produced instead of only an initial and final plot
        pylab.pause(0.00001)
        
    if t == 0: #if the time variable is satisfied or the time period finishes.
        pylab.savefig("initialplot.png")
        ### you might want to save more of them
        
pylab.savefig("finalplot.eps")
# So the windows don't auto-close at the end if run outside ipython
pylab.ioff()
pylab.show()


