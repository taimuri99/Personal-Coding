import numpy 

class Atom:
        def __init__(self,name,atnum,wght) :
            self.n=name
            self.z=atnum #atomic number
            self.w=wght
            
        def get_neutron_content(self):  
            proton_content=self.z*1.008
            self.neutcount=self.w - proton_content
            
### == Actual program starts here == ###
### read the file
fin=open('atomtab.dat','r')
lines=fin.readlines()
fin.close()

#create list of atoms
atoms=[]
for line in lines[1:]:
    l=line.split()
    name=l[0] ; atnum=int(l[1]) ; wght=float(l[2])
    at=Atom(name,atnum,wght)
    at.get_neutron_content()
    atoms.append(at)
    
#print neutron content per atoms in a file
fout=open('neutron_cont.dat','w')
for at in atoms:    
    print>>fout,at.n,at.neutcount
fout.close()
