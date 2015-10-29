#!/usr/bin/env python
import numpy as np
import sys
from scipy.integrate import quad
f,df=float(sys.argv[1]),float(sys.argv[2])
fL=f-df/2.
fH=f+df/2.
fmLow=88e6
fmHigh=108e6
def g(x):
    if(x>=fmLow and x<fmHigh):
        return 1.
    else:
        return 0. 

def f(x):
    if(x>=fL and x<fH):
        return 1.
    else:
        return 0.
print quad(lambda x: f(x)*g(x),fLtemp,fHtemp)[0]/(fH-fL)
