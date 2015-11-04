#!/usr/bin/env python
import numpy as np
import sys
z,dz=float(sys.argv[1]),float(sys.argv[2])
f21=1.42040575177 #21 cm frequency in GHz
print '%.5f'%(f21*dz/(1+z)**2.)

