#!/usr/bin/python

def f(y2,y3):
    h3=y3-y2
    h4=1-y2
    h1=y2
    return abs(-2*h3*h4+2*h1*(h3+h4))

q=40

while(q<10000):
    d=1./q

    s=0
    for i in range(q+1):
        for j in range(i,q+1):
            s+=f(i*d,j*d)

    print(d,s*d*d)
    q+=q//2
