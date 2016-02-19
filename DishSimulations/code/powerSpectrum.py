#************************************************************
#Aaron Ewall-Wice
#Power spectrum tools
#aaronew@mit.edu
#May 14th 2015
#************************************************************
import numpy as np
import numpy.fft as fft

def fft3(data):
    return fft.fftshift(fft.fftn(fft.fftshift(data)))


#compute binned 1D power spectrum
def ps1D(data,dx,kMin,kMax,nk,logbins=True):
    nSide=data.shape[0]
    print 'nSide=%d'%(nSide)
    print 'dx=%.2f Mpc'%(dx)
    dk=2.*np.pi/(nSide*dx)
    print 'dk=%.2f/Mpc'%(dk)
    psCube=np.abs(fft3(data))**2.*dx**3./(nSide**3.)
    #generate kaxis
    if(logbins):
        kMin=np.log10(kMin)
        kMax=np.log10(kMax)
    dkBin=(kMax-kMin)/nk
    print 'dkb=%f'%(dkBin)
    print 'kMax=%f'%(kMax)
    print 'kMin=%f'%(kMin)
    kAxis=dk*np.arange(-nSide/2.,nSide/2.)
    kCubeX,kCubeY,kCubeZ=np.meshgrid(kAxis,kAxis,kAxis)
    kCube=np.sqrt(kCubeX**2.+kCubeY**2.+kCubeZ**2.)
    del kCubeX
    del kCubeY
    del kCubeZ
    if(logbins):
        binCube=np.round((np.log10(kCube)-kMin)/dkBin).astype(int)
    else:
        binCube=np.round((kCube-kMin)/dkBin).astype(int)
    #now bin
    kCube=kCube.flatten()
    binCube=binCube.flatten()
    psCube=psCube.flatten()
    ps=np.zeros(nk)
    kaxis=np.zeros(nk)
    counts=np.zeros(nk)
    for mm in range(nk):
        bins=binCube==mm
        counts[mm]=len(binCube[bins])
        ps[mm]=np.mean(psCube[bins])
        kaxis[mm]=np.mean(kCube[bins])
    return kaxis,counts,ps


    
        

    
    
    
    
