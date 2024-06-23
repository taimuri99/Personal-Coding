#!/usr/bin/env python
# coding: utf-8

# In[10]:


from numpy import*
import numpy as np
import matplotlib.pyplot as pl

from matplotlib.pyplot import *


# In[11]:


class Students:
    def __init__(self,name,marks,month):
        self.n=name
        self.m=marks 
        self.d=month


# In[12]:


fin=open('firstdata.dat','r')
lines=fin.readlines()  #datafile used for the values required in the data presentation
fin.close()


# In[15]:


markslist=[] #arrays initialised empty to which values appended and then plotted
monthlist=[]
for line in lines[0:41]: #python starts from 0
    l=line.split()
    name=l[0] ; marks=int(l[1]) ; month=int(l[2])
    markslist.append(marks)
    monthlist.append(month) #values appended to the initialised arrays
aver_marks=np.mean(markslist)
aver_month=np.mean(monthlist)

print(aver_marks,aver_month) 


# In[19]:


bins=[]
for b in range(max(monthlist)+1):bins.append(b+0.5)
pl.hist(monthlist,bins)
pl.savefig('months.png')
pl.suptitle('Histogram')
pl.xlabel('Months')
pl.ylabel('Frequency')
pl.show()


# In[ ]:




