import numpy as np
from pylab import *

# RC circuit following Kirchoff's law
# Write here the values of resistor/capacitor/applied voltage
Res = 20. # resistance (Ohm)
Cap = 0.1 # capacitance (Farad)
Volt = 10. # applied voltage (Volt)
charge = Volt*Cap #at long t for charging; at t=0 for discharging
curr = Volt/Res #Ohm's law
print 'Reference charge [C]=',charge, 'Resistence [Ohm]=',Res, 'Capacity [F]=', Cap
print 'Max current [Ampere]=',curr

#Integrating factor/integrating ODEs :define the characteristic time of the RC-circuit
#meaning how quickly the capacitor charges/discharges (RC in seconds)
tau = Res*Cap #resistance times capacitance
print 'Characteristic time [sec]', tau

t=np.linspace(0,20,num=10)  # time between 0-10 sec time step 0.1
#t=np.arange(0,10,0.1)  # time between 0-10 sec time step 0.1 -equivalent Python instruction

def expo(char_time):
    expo = exp(-t/char_time)
    return expo
        
def CHARGE(charge,current,Vappl,tau):
    ""   
#    RC-charging When a voltage source is applied to an RC circuit, the capacitor, C charges up through the resistance, R
    ""
    IntFact = expo(tau)
    qc=charge*(1.-IntFact)         # charge versus the time
    currc=current*tau*IntFact       # current i = dq/dt
    biasc=Vappl*(1.-IntFact)       # voltage across the capacitor
    return qc,currc,biasc
        
def DISCHARGE(charge,current,Vappl,tau):
    IntFact = expo(tau)
    qdis=charge*IntFact    # charge versus the time
    currdis= -current*tau*IntFact # current i = dq/dt
    biasdis=Vappl*IntFact  #voltage 
    return qdis,currdis,biasdis

        
q=DISCHARGE(charge,curr,Volt,tau)[0]
curr=DISCHARGE(charge,curr,Volt,tau)[1]
bias=DISCHARGE(charge,curr,Volt,tau)[2]

plot(t,q,'ro')
plot(t,curr,'bo')
plot(t,bias,'go')
xlabel("time (seconds)")
ylabel("Capacitor Charge (red), Current (blue), Voltage (green)")
#caption("Red line for charge; Blue line for current; gree stands for voltage")
title("RC circuit")
show()
