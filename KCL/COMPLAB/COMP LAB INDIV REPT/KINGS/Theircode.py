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
# Please enter here your name: STUDENT NAME
#
#============================================================================
#  
import numpy as np
import pylab
#
# Set pylab to interactive mode so plots update when run outside ipython
pylab.ion()
#
# On your report add a sentence commenting various interactive graphilcal mode you have used during the module 
#=============================================================================
# Utility functions
#
def Gaussian(x,t,sigma):
    """  A Gaussian curve.
        x = Variable
        t = time shift
        sigma = standard deviation      """
    return np.exp(-(x-t)**2/(2*sigma**2))
def free(npts):
    """This utility function creates a list of zeros of length npts. 
    To be used for a free particle and in all the region where the particle is free
    """
    return np.zeros(npts)
def step(npts,v0):
    "Add a function to mimick a step potential of height v0. Return the potential v"
    return v
def barrier(npts,v0,thickness):
    "Add a function to mimick a barrier of height v0 and thickness thickness"
    return v
    

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
#  REPORT IN THE CODE AND IN THE REPORT WHAT INOUT PARAMETER YOU NEED AND PASS THEM FROM INPUT
#  ON YOUR REPORT: say clearly which data you collect per each set of input parameters. A table coudl be useful
#
''''
N    =          #
T    =          #  Make sure that the number of total steps avoid to reach the boundary. 
Tp   =          #  
dx   = 1.0e0    #  
m    = 1.0e0    #  
hbar = 1.0e0    #  
#
# 
V0   =   
THCK =        
POTENTIAL =          
sigma =                 
''''    
X    = dx*np.linspace(0,N,N)        
x0 = round(N/2) - 5*sigma   
k0 = 2.*np.pi/sigma         
E = (hbar**2/2.0/m)*(k0**2+0.5/sigma**2) 
#=============================================================================
# Code begins
if POTENTIAL=='free':
    V = free(N)
#
# ADD HERE ALL THE OTHER TYPE OF POTENTIAL YOU WILL IMPLEMENT in agreement with utitlity function deintion
# please make sure that on your report you sketch those potentials
#
elif POTENTIAL==''  # for a step potential
elif POTENTIAL==''  # for a double barrier 
elif POTENTIAL==''  # for a well
else:
    raise ValueError("Unrecognized potential type: %s" % POTENTIAL)
#
#  EXPLAIN IN THE REPORT CLEARLY THE ALGORITHM IMPLEMENTATION. Check the correspondence between
#  equations and variable names in the code. The maximum stable time step is a function of the potential, V
#  COMMENT ON THE choice of time step, true everytime you have a dynamics and the time is discretised
#
Vmax = V.max()                                 #
dt   = hbar/(2*hbar**2/(m*dx**2)+Vmax)         #  
c1   = hbar*dt/(m*dx**2)                       #  
c2   = 2*dt/hbar                               #  
c2V  = c2*V                                    # p
#
#
#
print 'One-dimensional Schrodinger equation - time evolution'
print 'Wavepacket energy:   ',E
print 'Potential type:      ',POTENTIAL
print 'Potential height V0: ',V0
print 'Barrier thickness:   ',THCK
print 'Well thickness:      ',THCK  
#
#============ Algorithm implementation: to be clearly explained in the report
#
psi_r = np.zeros((3,N))     #  Real
psi_i = np.zeros((3,N))     #  Imaginary
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
psi_r[PR,xn] = cx*gg
psi_i[PR,xn] = sx*gg
psi_r[PA,xn] = cx*gg
psi_i[PA,xn] = sx*gg
#
#
psi_p = psi_r[PR]**2 + psi_i[PR]**2
#
#  Normalization: what you shoudl do?
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
lineI, = pylab.plot(X,psi_i[PR],'r',alpha=0.7,label='Imag')
lineP, = pylab.plot(X,6*psi_p,'k',label='Prob')
pylab.title('Potential height: %.2e' % V0)
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
#
#=============================================================Time evolution
#
for t in range(T+1):
    # Precompute a couple of indexing constants, this speeds up the computation
    psi_rPR = psi_r[PR]
    psi_iPR = psi_i[PR]
    #  Apply the update equations.
    psi_i[FU,IDX1] = psi_i[PA,IDX1] + \
                      c1*(psi_rPR[IDX2] - 2*psi_rPR[IDX1] +
                          psi_rPR[IDX3])
    psi_i[FU] -= c2V*psi_r[PR]

    psi_r[FU,IDX1] = psi_r[PA,IDX1] - \
                      c1*(psi_iPR[IDX2] - 2*psi_iPR[IDX1] +
                          psi_iPR[IDX3])
    psi_r[FU] += c2V*psi_i[PR]
    #  Increment the time steps.  PR -> PA and FU -> PR
    psi_r[PA] = psi_rPR
    psi_r[PR] = psi_r[FU]
    psi_i[PA] = psi_iPR
    psi_i[PR] = psi_i[FU]
    #
     #
    if t % Tp == 0:
        #  Compute observable probability for the plot
        #  REPORT: Explain what you ar plotting and what is the probability
        #  This is mandatory to understand what data you collect
        psi_p = psi_r[PR]**2 + psi_i[PR]**2
        #  Update the plots.
        lineR.set_ydata(psi_r[PR])
        lineI.set_ydata(psi_i[PR])
        # Note: we plot the probability density amplified by a factor so it's a bit easier to see.
        lineP.set_ydata(6*psi_p)
#
        pylab.draw()
    if t == 0:
        pylab.savefig("graph_ini.png")
        ### you might want to save more of them
        
pylab.savefig("graph.eps")
# So the windows don't auto-close at the end if run outside ipython
pylab.ioff()
pylab.show()

