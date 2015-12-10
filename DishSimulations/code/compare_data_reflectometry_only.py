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
fAxisRB_PAPER=n.arange(1601)*200e6/1e9/1601.+150e-3
#between 100 and 200 MHz

def fromVNAHP(filename):
    data=n.loadtxt(filename,skiprows=9,delimiter=',')
    return data[:,0]+1j*data[:,1]

s11FeedOnlyRB=fromVNAHP(reflectometryRBDir+'/FeedOnly/DATA02.d1') #S11 of the feed only pointing straight up, effectively measures Gamma_f (the feed reflection coefficient)
s11FeedAndDishRB=fromVNAHP(reflectometryRBDir+'/DishAndFeed/HERA_Dish1_S11_feed5_0m_A.d1') #S11 of the feed suspended over the dish.
PAPERSelection=n.logical_and(fAxisRB_PAPER>=.1,fAxisRB_PAPER<.2)

s11PAPER=fromVNAHP(reflectometryRBDir+'/PAPER/DATA01.d1')
s11PAPERDeep=fromVNAHP(reflectometryRBDir+'/PAPER/DATA03.d1')


fAxisRB_PAPER=fAxisRB_PAPER[PAPERSelection]
s11PAPER=s11PAPER[PAPERSelection]
s11PAPERDeep=s11PAPERDeep[PAPERSelection]

windowRB_PAPER=signal.blackmanharris(len(s11PAPER))
windowRB_PAPER=windowRB_PAPER/windowRB_PAPER.mean()

#generate PAPER delay
tAxisPAPER=fft.fftshift(fft.fftfreq(len(s11PAPER),fAxisRB_PAPER[1]-fAxisRB_PAPER[0]))
responsePAPERFreq=(1.-s11PAPER)
responsePAPERFreqDeep=(1.-s11PAPERDeep)


delayResponsePAPER=fft.fftshift(fft.ifft(fft.fftshift(responsePAPERFreq*windowRB_PAPER)))
delayResponsePAPER[tAxisPAPER<0]=0.


delayKernelPAPER=signal.fftconvolve(delayResponsePAPER,n.conj(delayResponsePAPER)[::-1],mode='same')
#delayKernelPAPER=fft.fftshift(fft.ifft(fft.fftshift(n.abs(responsePAPERFreq)**2.*windowRB_PAPER)))
#delayKernelPAPERDeep=fft.fftshift(fft.ifft(fft.fftshift(n.abs(responsePAPERFreqDeep)**2.*windowRB_PAPER)))
#delayKernelPAPERCross=fft.fftshift(fft.ifft(fft.fftshift(n.conj(responsePAPERFreqDeep)*responsePAPERFreq*windowRB_PAPER)))


p.plot(tAxisPAPER,10.*n.log10(n.abs(delayKernelPAPER)/n.abs(delayKernelPAPER).max()),label='PAPER Reflectometry RB')
#p.plot(tAxisPAPER,10.*n.log10(n.abs(delayKernelPAPER)/n.abs(delayKernelPAPERDeep).max()),label='PAPER Reflectometry RB Deep')
#p.plot(tAxisPAPER,10.*n.log10(n.abs(delayKernelPAPER)/n.abs(delayKernelPAPERCross).max()),label='PAPER Reflectometry RB Cross')


#generate delay response from Rich's measurements using S11 of the feed and dish and the feed only
#and the equations in Nipanjana's paper and memo drafts 
responseFeedAndDishFreqRB=(s11FeedAndDishRB-s11FeedOnlyRB)*s11FeedOnlyRB/(1.-s11FeedOnlyRB)+(1.-s11FeedOnlyRB) #This is the frequency dependent gain. 
windowRB=signal.blackmanharris(len(responseFeedAndDishFreqRB))#windowing function 
windowRB=windowRB/windowRB.mean()#normalize window
tAxisFeedAndDishRB=fft.fftshift(fft.fftfreq(len(responseFeedAndDishFreqRB),fAxisRB[1]-fAxisRB[0]))
delayResponseFeedAndDishRB=fft.fftshift(fft.ifft(fft.fftshift(responseFeedAndDishFreqRB*windowRB)))

delayKernelFeedAndDishRB=signal.fftconvolve(delayResponseFeedAndDishRB,n.conj(delayResponseFeedAndDishRB)[::-1],mode='same')
#delayKernelFeedAndDishRB=fft.fftshift(fft.ifft(fft.fftshift(n.abs(responseFeedAndDishFreqRB)**2.*windowRB)))#power kernel, the DFT of the square of the frequency dependent gain

p.plot(tAxisFeedAndDishRB,10.*n.log10(n.abs(delayKernelFeedAndDishRB)/n.abs(delayKernelFeedAndDishRB).max()),color=[.5,.5,.5],lw=2,ls='-',label='RB Reflectometry')

##load and plot a closed circuit measurement that Rich took.
#responseClosed=fromVNAHP(reflectometryRBDir+'/FeedOnly/DATA00.d1')#S11 of closed circuit
#delayKernelClosed=fft.fftshift(fft.ifft(fft.fftshift(n.abs(responseClosed)**2.*windowRB)))
#p.plot(tAxisFeedAndDishRB,10.*n.log10(n.abs(delayKernelClosed)),label='closed circuit')



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

#load simulated autocorrelation
#autocorr3=n.load('threeZoneAutocorr.npz')
#p.plot(autocorr3['delays']*1e9,10.*n.log10(autocorr3['autocorrDelay']),color='k',lw=1,label='three zone\nbeam model')


#************************************************************
#plot Rich's PAPER Reflectometry
#************************************************************


p.legend(ncol=2)
p.grid()
p.xlabel('$\\tau$ (ns)')
p.ylabel('Power (dB)')
p.xlim(-30,1000)
p.gcf().set_size_inches(10,10)
p.savefig('compare_data_reflectometry.pdf')

