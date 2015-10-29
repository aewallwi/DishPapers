import numpy as np
import scipy.interpolate as interp
import scipy.integrate as integrate
import cosmology as cosmo
sigmaT=6.6524587e-25
nfList=np.loadtxt('nfList.txt')
zList=np.loadtxt('zList.txt')
nfFunc=interp.interp1d(zList,nfList)
def intFunc(x):
    if(x<zList.min()):
        return cosmo.baryon_n(x)*sigmaT*cosmo.DHcm*cosmo.de(x)/(1.+x)
    elif(x>=zList.min() and x<=zList.max()):
        return cosmo.baryon_n(x)*sigmaT*cosmo.DHcm*cosmo.de(x)/(1.+x)*(1.-nfFunc(x))
    else:
        return 0.
        
tauSum=quad(intFunc,0,zList.min())[0]+quad(intFunc,zList.min(),zList.max())[0]
print tauSum
