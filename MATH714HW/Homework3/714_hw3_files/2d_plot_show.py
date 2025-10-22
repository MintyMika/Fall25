#!/usr/bin/python3
import numpy as np
from skimage import io
import matplotlib.pyplot as plt

# Load in the test image 
a=io.imread("van_vleck_map.png",as_gray=True)

# Check size
(y,x)=a.shape
print(x,"by",y,"pixels")

# Plot the image
plt.imshow(a,cmap=plt.cm.gray)
plt.show()
