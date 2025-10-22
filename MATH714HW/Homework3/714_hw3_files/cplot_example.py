#!/usr/bin/python3
import numpy as np
import custom_plot as cplt
from math import sin

# Load in the wall matrix
q=np.loadtxt("van_vleck.txt",dtype=np.int8)

# Create an empty array and fill it with some made-up data
w=np.zeros((73,160))
for i in range(73):
    x=0.1*i
    for j in range(160):
        y=0.1*j
        w[i,j]=sin(x+0.5*y)

# Call the first custom plotting routine
cplt.plot1("test1.png",w,q,0,1,4)

# Call the second custom plotting routine
cplt.plot2("test2.png",w,q,-1.1,1.1,4)
