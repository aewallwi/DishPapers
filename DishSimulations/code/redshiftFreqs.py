#!/usr/bin/env python
import numpy as np
import sys
from scipy.integrate import quad

f21=1.42040575177e9 #21 cm frequency
def dz2df(dz,z):
    return f21*dz/(1+z)**2.
def z2f(z):
    return f21/(1+z)
zmin,zmax,dz,odir=sys.argv[1:]
zmin=float(zmin)
zmax=float(zmax)
dz=float(dz)
#now divy up redshift into dz intervals
zList=np.arange(zmin,zmax+dz,dz)
fmFrac=np.zeros(len(zList))
fCenters=np.zeros(len(zList))
dFs=np.zeros(len(zList))
#compute center frequencies
fCenters=z2f(zList)
dFs=dz2df(dz,zList)
#compute fractional coverage within FM
fmLow=88e6
fmHigh=108e6
def g(x):
    if(x>=fmLow and x<fmHigh):
        return 1.
    else:
        return 0. 
for mm in range(len(fCenters)):
    fLtemp=fCenters[mm]-dFs[mm]/2.
    fHtemp=fCenters[mm]+dFs[mm]/2.
    #print fLtemp/1e6
    #print fHtemp/1e6
    def f(x):
        if(x>=fLtemp and x<fHtemp):
            return 1.
        else:
            return 0.
#    print f(fCenters[mm])*g(fCenters[mm])
    print quad(lambda x: f(x)*g(x),fLtemp,fHtemp)
    fmFrac[mm]=quad(lambda x: f(x)*g(x),fLtemp,fHtemp)[0]/(fHtemp-fLtemp)
    print fmFrac[mm]
output=np.round(np.vstack([zList,fCenters/1e9,dFs/1e9,fmFrac]).T,decimals=5)
#save to text filed
np.savetxt(odir+'/bandInfo.txt',output,fmt='%.5f')


        


