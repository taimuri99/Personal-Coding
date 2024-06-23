#!/usr/bin/env python
# coding: utf-8

# In[15]:


from numpy import*
import numpy as np

from matplotlib.pyplot import *


# In[16]:


from numpy.random import randint #for outputting random numbers(random.randint)


# In[17]:


J=0.0


# In[18]:


for i in range(10000): #iterated 10000 times 
    x = randint (1,365,size=(30)) #30 diff numbers if repeat, 1 added to J
    seen = set()
    y = list(set(n for n in x if n in seen or seen.add(n)))
    if len(y): 
        J=J+1.0
    else: 
        J=J  #if the birthday is repeated more than once, it will add to the variable J which is initialised as 0


# In[19]:


l=(J/i)*100
print l


# In[ ]:




