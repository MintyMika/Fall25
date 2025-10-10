#!/usr/bin/python3
from math import exp,tan,cos
from random import random
import sys

def f(z):
    return exp(-z)*tan(z)

def ddf(z):
    tanz=tan(z)
    secz=1./cos(z)
    return exp(-z)*(2*secz*(-1+tanz)+tanz)

def deriv(x1,x2,x3,x4):
    h1=x2-x1
    h2=x3-x2
    h3=x4-x3

    a=-2*(2*h2+h3-h1)/(h2*(h2+h3)*h1)
    b=2*(h2+h3-h1)/(h2*h3*(h2+h1))
    c=-2*(h2-h1)/(h3*(h2+h3)*(h2+h3+h1))
    d=2*(2*h2+h3)/(h1*(h2+h1)*(h2+h3+h1))

    return f(x2)*a+f(x3)*b+f(x4)*c+f(x1)*d


for i in range(100,301):
    H=10**(-i/100)
    x2=random()*H
    x3=random()*H
    if x2>x3: (x2,x3)=(x3,x2)

    print(H,abs(deriv(0,x2,x3,H)-ddf(x2)))

