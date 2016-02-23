#!/usr/bin/env python
#bw in GHz
#df in MHz
import sys
import numpy as np
bw,df=float(sys.argv[1]),float(sys.argv[2])
nchan=int(np.round(bw*1e3/df))
if(np.mod(nchan,2)==1):
    nchan+=1
print '%d'%(nchan)
