#!/usr/bin/env python
import numpy as np
from numpy import fft as fft
import matplotlib.pyplot as plt
import scipy.interpolate as interp
#************************************************************
#1d fft convolution
#************************************************************
def convFFT(x,y):
    nx=len(x)
    xf=np.pad(x,(nx/2,nx/2),mode='constant')
    yf=np.pad(y,(nx/2,nx/2),mode='constant')
    xf=(fft.fft(fft.fftshift(xf)))
    yf=(fft.fft(fft.fftshift(yf)))
    return fft.fftshift(np.real((fft.ifft(xf*yf))))[nx/2:3*nx/2]
def dconvFFT(x,y):
    nx=len(x)
    xf=np.pad(x,(nx/2,nx/2),mode='constant')
    yf=np.pad(y,(nx/2,nx/2),mode='constant')
    xf=(fft.fft(fft.fftshift(xf)))
    yf=(fft.fft(fft.fftshift(yf)))
    return np.real(fft.fftshift(fft.ifft(xf/yf)))[nx/2:3*nx/2]

#************************************************************
#bandpass derived from Rich Bradley's simulation of inital 
# feed design. The delay response is well modeled by a gaussian centered
#at 40 nanoseconds with a power law tail. 
#************************************************************
def heraFeedDelay(tau):
    cutoff=(80.-30.)*1e-9
    mu=(72.-30.)*1e-9
    vals1=np.abs(tau)<=cutoff
    print vals1
    vals2=np.abs(tau)>cutoff
    pwise=np.zeros(len(tau))
    func1=lambda x: np.exp(-(np.abs(x)-mu)**2./(2.*(7e-9)**2.))
    func2=lambda x: func1(tau[vals1].max())*((np.abs(x)-cutoff)*1e9+1.)**-3.
    pwise[vals1]=func1(tau[vals1])
    pwise[vals2]=func2(tau[vals2])
    pwise[tau<=0.]=0.
    #plt.plot(tau*1e9,pwise)
    #plt.xlim(1e-2,1e3)
    #plt.yscale('log')
    #plt.show()
    return pwise
          




#def initHeraFeedBand(f,params,sRate=5e9):
#    df=params[0]
#    bw=params[1]
#    fc=params[2]
#    print df,bw,fc
#    nf=int(np.round(bw/df))
#    if(np.mod(nf,2)==1):
#        nf+=1
#    #find oversampling factor by dividing bandwidth by sampling rate
#    sFactor=int(sRate/bw)
#    dt=1./sRate
#    nt=nf*sFactor
#    tau=np.arange(-nt/2,nt/2)*dt
#    kernel=(np.sinc(np.pi*tau*bw)*np.exp(-1j*2*np.pi*tau*fc))
##    kernel=kernel/np.sqrt(np.mean(np.abs(kernel)**2.))
#    delaySpectrum=convFFT(heraFeedDelay(tau),kernel)
#    band=np.fft.fftshift(np.fft.fft(np.fft.fftshift(delaySpectrum)))
#    find=(np.round(f/df)).astype(int)+nt/2
#    breturn=band[find]
#    breturn=breturn/breturn[-1]
#    plt.plot(f,np.abs(breturn))
#    plt.show()
#    return f,breturn

#params are a list of frequencies and complex gains to be interpolated
def initHeraFeedBand(params):
    f=params[0]
    b=params[1]
    fc=params[2]
    bfReal=interp.interp1d(f,np.real(b))
    bfImag=interp.interp1d(f,np.imag(b))
    nFactor=bfReal(fc)+1j*bfImag(fc)
    bfunc=lambda(x): (bfReal(x)+1j*bfImag(x))/nFactor
    return bfunc

def heraFeedBand(f,params):
    return params[0](f)


#def heraFeedBand(f,params):
 #   fc=params[0]
 #   df=params[1]
 #   breturn=params[2]
 #   find=np.array(np.round((f-fc)/df)).astype(int).flatten()+len(breturn)/2
 #   return breturn[find]
    

#************************************************************
#flat bandpass with a constant arbitrary gain given by params[0]
#************************************************************
def flatband(f,params):
    g=params[0]
    return g

#************************************************************
#band pass with a single reflection
#************************************************************
def reflectionsimple(f,params):
    r=params[0]
    phi=params[1]
    tau=params[2]
    return 1./(1.-r*np.exp(1j*(2.*np.pi*tau*f+phi)))

#************************************************************
#use discrete fourier series to define bandpass
#params containes complex fourier coefficient and delay for 
#series
#params=[r0,tau0,r1,tau1,...]
#dc term must also be specified!
#************************************************************
def fourierSeries(f,params):
    nref=len(params)/2
    output=np.zeros(len(f))
    for mm in range(nref):
        output+=params[mm*2]*np.exp(1j*params[mm*2+1]*f)
    return output

#************************************************************
#represents bandpass for a tile. wraps a bandpass function 
#(bfunc) which accepts f in Hz and paramters
#************************************************************
class Bandpass():
    def __init__(self,bfunc,params):
        if(bfunc.__name__=='heraFeedBand'):
#            df=params[0]
#            bw=params[1]
#            fc=params[2]
#            nf=int(bw/df)
#            if(np.mod(nf,2)==1):
#                nf+=1
#            faxis=np.arange(-nf/2,nf/2)*df+fc
#            faxis,breturn=initHeraFeedBand(faxis,params,sRate=5e9)
#            params=[fc,df,breturn]
#            self.bfunc=lambda f: bfunc(f,params)
#            self.params=params
            self.params=[initHeraFeedBand(params)]
            self.bfunc=lambda f: bfunc(f,self.params)
        else:
            self.bfunc=lambda f: bfunc(f,params)
            self.params=params
        
    def set_params(self,params):
        self.params=params
        self.bfunc=lambda f: bfunc(f,params)
    def gain(self,f):
        return self.bfunc(f)
