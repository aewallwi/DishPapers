#!/usr/bin/env python
import numpy as np
import glob,re,os
import scipy.interpolate as interp
import cosmology
zre=re.compile('z[0-9]{3}.[0-9]{2}')
zlist=[]
plist=[]
klist=[]

#now get lowres power
flist=glob.glob('Deldel_T_power_spec/ps_*')
cdir=os.path.dirname(os.path.realpath(__file__))
for fname in flist:
    z=float(zre.findall(fname)[0][1:])
    ps=np.loadtxt(cdir+'/'+fname)
    plist=np.concatenate((plist,ps[:,1]))
    klist=np.concatenate((klist,ps[:,0]))
    zlist=np.concatenate((zlist,z*np.ones(ps.shape[0])))
print klist
#maxk=10.
#plist=plist[klist<maxk]
#zlist=zlist[klist<maxk]
#klist=klist[klist<maxk]
'''
flist=glob.glob('pspec/ps_*')
for fname in flist:
    z=float(zre.findall(fname)[0][1:])
    ps=np.loadtxt(cdir+'/'+fname)
    kfilter=ps[:,0]>maxk
    plist=np.concatenate((plist,ps[kfilter,1]))
    klist=np.concatenate((klist,ps[kfilter,0]))
    zlist=np.concatenate((zlist,z*np.ones(ps.shape[0])[kfilter]))
maxk=klist.max()/2.
plist=plist[klist<maxk]
zlist=zlist[klist<maxk]
klist=klist[klist<maxk]
flist=glob.glob('pspec_vhighres/ps_*')
for fname in flist:
    z=float(zre.findall(fname)[0][1:])
    ps=np.loadtxt(cdir+'/'+fname)
    kfilter=ps[:,0]>maxk
    plist=np.concatenate((plist,ps[kfilter,1]))
    klist=np.concatenate((klist,ps[kfilter,0]))
    zlist=np.concatenate((zlist,z*np.ones(ps.shape[0])[kfilter]))

'''

print zlist.shape
print klist.shape
print plist.shape
print zlist
klist=np.log10(klist)
print 'min klist'
print klist.min()
plist=np.log10(plist)
maxz=zlist.max()
minz=zlist.min()
points=np.vstack([klist,zlist]).T

def delta21(k,z):
    k10=np.log10(k)
    return 10.**(interp.griddata(points,plist,(k10,z),fill_value=0.,method='linear'))

#beampp is the integrated beam squared
def ps21(k,z,flux=False,beampp=1.,kmin=10**(klist.min())):
    if(flux):
        cfactor=cosmology.eta2kpara(1.,z)*(cosmology.u2kperp(1.,z)**2.)*(cosmology.t2i(1420e6/(z+1.),1e-3)**2.)*beampp
    else:
        cfactor=1.
    output=cfactor*2.*np.pi**2.*delta21(k,z)/k**3.
    throwouts=k<=kmin
    if(np.any(throwouts)):
        output[throwouts]=0.
    return output

def delta21_c(kperp,kpara,z):
    k=np.sqrt(kperp**2.+kpara**2.)
    return delta21(k,z)

def ps21_c(kperp,kpara,z,flux=False,beampp=1.):
    k=np.sqrt(kperp**2.+kpara**2.)
    return ps21(k,z,flux,beampp)

