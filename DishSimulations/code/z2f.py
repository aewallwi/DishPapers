#!/usr/bin/env python
#convert z to frequency in GHz. 
import sys
z=float(sys.argv[1])
f=1.42040575177/(z+1.)
print '%.5f'%(f)

