#!/usr/bin/env python
#************************************************************
#Aaron Ewall-Wice
#aaronew@mit.edu
#November 25th 2015
#************************************************************
# A script that compares many different estimates of the
# HERA dish power kernel.
# 1) Autocorrelations
# 2) Reflectometry measurements from UCB
# 3) Reflectometry measurements from Rich Bradley at Green Bank
# 4) Simulations
#************************************************************
import numpy as n
import numpy.fft as fft
import aipy as a
import capo.arp as arp
import matplotlib.pyplot as p
import scipy.signal as signal
import glob
import csv
#************************************************************
#hard code directories of files for now
#sorry :P
#************************************************************
autocorrDir='/Users/aaronew/Dropbox/HERA_autocorrs/' #parent directory of autocorrelations
simDir='/Users/aaronew/Dropbox/DishPapers/DishSimulations/' #parent directory for simulations
reflectometryDir='/Users/aaronew/Dropbox/DishPapers/DishReflectometry/'#parent directory for UCB reflectometery measurements
reflectometryRBDir='/Users/aaronew/Dropbox/DishPapers/Reflectometry_RB/'#parent directory for Rich Bradley's reflectometery measurements

#************************************************************
#first load up and plot autocorrelations. This code is
#primarily taken from Aaron Parsons arp/scripts/auto_delays.py
#************************************************************
thresh=.5
tol=1e-9
ants=[0,62,104,96]
fileList=glob.glob(autocorrDir+'*AR')
c0,c1=140,930
antstr=','.join(['%d_%d'%(i,i) for i in ants])
times,dat,flg=arp.get_dict_of_uv_data(fileList,antstr=antstr,polstr='xx',verbose=True)
print dat.keys()
div=False
colors=['']*10
g0,g1,g2={},{},{}
w0,w1,w2={},{},{}
for i,ant in enumerate(ants):
    print ant
    bl=a.miriad.ij2bl(ant,ant)
    print bl
    fqs=n.linspace(.1,.2,dat[bl]['xx'].shape[1])[c0:c1]
    tau=fft.fftfreq(fqs.size,fqs[1]-fqs[0])
    d,f=dat[bl]['xx'][:,c0:c1],flg[bl]['xx'][:,c0:c1]
    ntimes,nfreqs=d.shape
    w=n.logical_not(f).astype(n.float)
    d*=w
    wt=w.sum(axis=1);wt.shape=(-1,1)
    dt=n.sum(d,axis=1); dt.shape=(-1,1)
    d,w=d/dt*wt,w*wt
    d,w=d.sum(axis=0),w.sum(axis=0);d/=w.max();w/=w.max()
    ok=w>thresh
    g0[ant],w0[ant]=n.where(ok,d/w,0),n.where(ok,1,0)
    window=a.dsp.gen_window(d.size,'blackman-harris')
    dw,ww=d*window,w*window
    _dw,_ww=fft.ifft(dw),fft.ifft(ww)
    sync=(fqs/.150)**-4.5
    w2[ant]=w0[ant]
    g2[ant]=n.where(w2[ant]>0,g0[ant]/sync,0)
    _dw,_ww=fft.ifft(dw),fft.ifft(ww)
    gain=a.img.beam_gain(_ww)
    mdl,info=a.deconv.clean(_dw,_ww,tol=tol,stop_if_div=div)
    mdl+=info['res']/gain
    mdl/=mdl[0]
    p.plot(tau[:tau.size/2],10*n.log10(n.abs(mdl[:tau.size/2])),colors[i],label='%d/S'%(ant))
    mdl[:16]=0
    mdl[-15:]=0
    ns=2./n.sqrt(100e6/1024/2*10.7*ntimes)



#************************************************************
#plot simulations
#************************************************************

pwFeedAndDish=n.loadtxt(simDir+'data/feedSim_500ns_pW.txt',skiprows=2) #2xntimes array with col 0 = times (ns) and col 1 = plane wave e field.
vFeedAndDish=n.loadtxt(simDir+'data/feedSim_500ns_o1.txt',skiprows=2) #2xntimes array with col0 = times (ns) and col 1 = feed voltage output.
ntFeedAndDish=len(pwFeedAndDish)
dtFeedAndDish=vFeedAndDish[1,0]-vFeedAndDish[0,0]
#zero pad to avoid wrap around
pwFeedAndDish_pad=n.pad(pwFeedAndDish[:,1],(ntFeedAndDish/2,ntFeedAndDish/2),mode='constant')
vFeedAndDish_pad=n.pad(vFeedAndDish[:,1],(ntFeedAndDish/2,ntFeedAndDish/2),mode='constant')
#frequency axis,padded
fAxisFeedAndDish_pad=fft.fftfreq(len(vFeedAndDish_pad),dtFeedAndDish)
pwFeedAndDish_f=fft.fft(pwFeedAndDish_pad)
vFeedAndDish_f=fft.fft(vFeedAndDish_pad)
dfFeedAndDish=1./(len(vFeedAndDish_pad)*dtFeedAndDish)
#select frequencies between 100 and 200 MHz
flow=.1
fhigh=.2
selection=n.logical_and(fAxisFeedAndDish_pad>=flow,fAxisFeedAndDish_pad<fhigh)
kernelF=(n.abs(vFeedAndDish_f/pwFeedAndDish_f)[selection])**2.
tAxisFeedAndDishSelection=n.arange(-len(kernelF)/2.,len(kernelF)/2.)/(fhigh-flow)
windowFeedAndDish=signal.blackmanharris(kernelF.size)
windowFeedAndDish/=windowFeedAndDish.mean()
kernelFeedAndDish=fft.fftshift(fft.ifft(fft.fftshift(windowFeedAndDish*kernelF)));kernelFeedAndDish/=n.abs(kernelFeedAndDish).max()
p.plot(tAxisFeedAndDishSelection,10*n.log10(n.abs(kernelFeedAndDish)),color='k',ls='-',lw=2,label='CST Simulation')

#************************************************************
#load and plot measurements
#taken by Rich Bradley of the HERA dish and feed along with
#the feed only
#************************************************************
#************************************************************
#loads S11 data from HP8753D 06.12 VNA
#************************************************************
fAxisRB=n.arange(1601)*100e6/1e9/1601.+100e-3 #frequencies measured by VNA are 1601 channels evenly spaced
#between 100 and 200 MHz

def fromVNAHP(filename):
    data=n.loadtxt(filename,skiprows=9,delimiter=',')
    return data[:,0]+1j*data[:,1]

s11FeedOnlyRB=fromVNAHP(reflectometryRBDir+'/FeedOnly/DATA02.d1') #S11 of the feed only pointing straight up, effectively measures Gamma_f (the feed reflection coefficient)
s11FeedAndDishRB=fromVNAHP(reflectometryRBDir+'/DishAndFeed/HERA_Dish1_S11_feed5_0m_A.d1') #S11 of the feed suspended over the dish. 

#generate delay response from Rich's measurements using S11 of the feed and dish and the feed only
#and the equations in Nipanjana's paper and memo drafts 
responseFeedAndDishFreqRB=(s11FeedAndDishRB-s11FeedOnlyRB)*s11FeedOnlyRB/(1.-s11FeedOnlyRB)+(1.-s11FeedOnlyRB) #This is the frequency dependent gain. 
windowRB=signal.blackmanharris(len(responseFeedAndDishFreqRB))#windowing function 
windowRB=windowRB/windowRB.mean()#normalize window
tAxisFeedAndDishRB=fft.fftshift(fft.fftfreq(len(responseFeedAndDishFreqRB),fAxisRB[1]-fAxisRB[0]))
delayKernelFeedAndDishRB=fft.fftshift(fft.ifft(fft.fftshift(n.abs(responseFeedAndDishFreqRB)**2.*windowRB)))#power kernel, the DFT of the square of the frequency dependent gain


p.plot(tAxisFeedAndDishRB,10.*n.log10(n.abs(delayKernelFeedAndDishRB)/n.abs(delayKernelFeedAndDishRB).max()),color='k',lw=2,ls=':',label='RB Reflectometry')

#************************************************************
#taken by UCB Green Bank Expidition. Most of this
#Load and plot inferred delay response from S11 Data
#since no feed only measurement exists for this data,
#the delay response assumes a frequency independent
#gamma_a which is estimated from the zero delay bin
#of the dft. 
#************************************************************
dfile_amp=reflectometryDir+'/alldata/FC_NC41_AB_DB.csv'
dfile_phs=reflectometryDir+'/alldata/FC_NC41_AB_P.csv'
#************************************************************
#This method is copied from
#from the DishReflectometry/plot_delay3.py
#************************************************************
def fromcsv(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = n.array(list(d)[18:-1], dtype=n.float)
    return x[:,0]/1e9, x[:,1]

dfreq,dfreqdb=fromcsv(dfile_amp)
dfreq,dphs=fromcsv(dfile_phs)

dSelection=n.logical_and(dfreq>=.1,dfreq<=.2)
dMeasured=(10.**(dfreqdb/20.)*n.exp(1j*dphs*n.pi/180.))[dSelection]
fMeasured=dfreq[dSelection]
tAxisMeasured=fft.fftshift(fft.fftfreq(len(dMeasured),fMeasured[1]-fMeasured[0]))
#windowMeasured=n.hamming(len(fMeasured))
windowMeasured=signal.blackmanharris(len(fMeasured))
windowMeasured=windowMeasured/windowMeasured.mean()
gAMeasured=fft.fftshift(fft.ifft(fft.fftshift(dMeasured*windowMeasured)))[len(dMeasured)/2]
responseMeasuredFreq=(dMeasured-gAMeasured)*(gAMeasured/(1-gAMeasured))+(1.-gAMeasured)
responseMeasured=fft.fftshift(fft.ifft(fft.fftshift(responseMeasuredFreq*windowMeasured)))
responseMeasuredKernel=signal.fftconvolve(responseMeasured,n.conj(responseMeasured)[::-1],mode='same')
p.plot(tAxisMeasured,10*n.log10(n.abs(responseMeasuredKernel/n.abs(responseMeasuredKernel).max())),color='k',ls='--',label='UCB Reflectometry')
#p.axvline(75.*1./3e8*1e9,label='negative delay of \nreflections originating at\n base of 75 meter cable',color='r',ls='--')
p.legend()
p.grid()
p.xlabel('$\\tau$ (ns)')
p.ylabel('Power (dB)')
p.xlim(-30,1000)
p.savefig('compare_data.pdf')
p.show()
