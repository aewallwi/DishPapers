#!/usr/bin/env python
import numpy as np
import sys
z=float(sys.argv[1])
h=.67
c=3e8
littleh=0.67
bigH=littleh*100./3.086e19#in inverse seconds
DH=c/bigH
Om0=0.28
OmL=0.32
f21=1.42040575177 #21 cm frequency in GHz
minDelay=float(sys.argv[2])#minimum delay in ns
print '%.5f'%(np.sqrt(OmL+Om0*(1.+z)**3.)(2.*np.pi*f21*minDelay)/(DH*(1.+z)**2.))


