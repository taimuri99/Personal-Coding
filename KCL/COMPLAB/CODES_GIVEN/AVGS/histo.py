#!/usr/bin/python
# python/examples/histo1.py
#
# Simple example of histogram plot 
#           by Tchilabalo Abozou Kpanzou
#          and David MacKay

import random
def points(n):
    # Make points by squaring x
    # where x is random in (0,1)
    Y = []
    for i in range(n):
        x = random.random()
        Y.append(x*x)  # x is a random number produced and y is the square of that random outputted x
    return Y

Y = points(10000) # Make a load of points
from pylab import *
figure()
hist(Y, 100) # The second argument is the number of bins
show()



