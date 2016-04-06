#! /usr/bin/env python

import numpy as np, pylab as plt, aipy as a
import sys, csv


def fromcsv(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = np.array(list(d)[18:-1], dtype=np.float)
    return x[:,0]/1e9, x[:,1]

def lin(db,ph): # Converting the magnitude and phase of the measurements into complex returnloss in voltage ratio.
    return 10**(db/20) * np.exp(2j*np.pi*ph/360)

def pwr(d): # Computing the absolute magnitude of power in db. 
    return 10*np.log10(np.abs(d)**2)
    
def vna_open(db_open,ph_open,fq):
    cable_phse = (2*np.pi*fq*30.48)/(0.81*0.3)
    ph_open_cable = ph_open+cable_phse
    d_open = lin(db_open,ph_open)
    meas = lin(db_open,ph_open_cable)
    gamma2 = (d_open-1)/(meas+1)
    return gamma2
    
def vna_short(db_open,ph_open,fq):
    cable_phse = (2*np.pi*fq*30.48)/(0.85*0.3)
    ph_open_cable = ph_open+cable_phse
    d_open = lin(db_open,ph_open)
    meas = lin(db_open,ph_open_cable)
    gamma2 = (d_open+1)/(meas-1)
    return gamma2

for filename in sys.argv[1:]:
    BASE = filename[:-len('.csv')]
    db_file = BASE + '_DB.csv'
    ph_file = BASE + '_P.csv'
    
    fq,db = fromcsv(db_file) # Reading the magnitude and the phase of the datafile to be processed
    fq,ph = fromcsv(ph_file) 
    
    fq_f,db_f = fromcsv('HERA_FEED_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_f,ph_f = fromcsv('HERA_FEED_P.csv')
    
    fq_open,db_open = fromcsv('Calibration/OPENNN_DB.csv')# Reading the open file to calibrate for the multiple signal reflections in the 50 feet long cable
    fq_open,ph_open = fromcsv('Calibration/OPENNN_P.csv')
    
    fq_short,db_short = fromcsv('Calibration/SHORTT_DB.csv')# Reading the open file to calibrate for the multiple signal reflections in the 50 feet long cable
    fq_short,ph_short = fromcsv('Calibration/SHORTT_P.csv')

        
    d = lin(db,ph) 
    d_f = lin(db_f,ph_f) 
    d_open = lin(db_open,ph_open)  
    d_short = lin(db_short,ph_short)
    
    
    valid = np.where(np.logical_and(fq < .200, fq > .100)) # restrict to PAPER band
    
    fq, d, db, ph = fq[valid], d[valid], db[valid], ph[valid]
    
    fq_f,d_f, db_f, ph_f = fq_f[valid], d_f[valid], db_f[valid],ph_f[valid]
    
    fq_open,d_open, db_open, ph_open = fq_open[valid], d_open[valid], db_open[valid],ph_open[valid]
    
    fq_short,d_short, db_short, ph_short = fq_short[valid], d_short[valid], db_short[valid],ph_short[valid]
    

    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    
    
    WINDOW = 'blackman-harris'
    window = a.dsp.gen_window(fq.size, WINDOW)
    
   

    if True: 
    
        #Computing VNA input return loss from the open data. 
        
        gamma_vna_open = vna_open(db_open,ph_open,fq)
        #plt.plot(fq,pwr(gamma_vna_open),linewidth=2.5,label='Feed only')
        #plt.plot(fq,pwr((1+gamma_vna_open)),linewidth=2.5,label='Feed only')
        #gamma_vna_short = vna_short(db_short,ph_short,fq)
         
        
        
        d_t =1+d_f # Computing the transmission coefficient
        
        d1 = ((d-d_f)*d_f/d_t ) # Reflections corrected to produce system bandpass in the receiving mode 
        
        d2 = d1+d_t # Zero term added for final correction

        _d = np.fft.ifft(np.abs(d)**2*window)/window.mean()
        _df = np.fft.ifft(np.abs(d_f)**2*window)/window.mean() #Take the Fourier transform of the SQUARE of the absolute MAGNITUDE of the voltage return loss of the feed.
        _dt = np.fft.ifft(np.abs(d_t)**2*window)/window.mean() #Take the Fourier transform of the SQUARE of the absolute MAGNITUDE of the voltage transmission of the feed.
        _dopen = np.fft.ifft(np.abs(d_open)**2*window)/window.mean() #Take the Fourier transform of the SQUARE of the absolute magnitude of the voltage return loss of the open.
        _dgamma2 = np.fft.ifft(np.abs(gamma_vna_open)**2*window)/window.mean() #Take the Fourier transform of the SQUARE of the absolute magnitude of the voltage return loss of the vna input.
        _dH = np.fft.ifft(np.abs(d2)**2*window) /window.mean()#Take the Fourier transform of the SQUARE of the absolute magnitude of the voltage return loss  of the HERA element.
        
        #_dopen[[np.where((np.fft.fftshift(tau) <= -90) & (np.fft.fftshift(tau) >= -170))]]=0.0
        
    	#plt.plot(np.fft.fftshift(tau),pwr(np.fft.fftshift(np.abs(_d))),linewidth=2.5,label='feed')
    	plt.plot(np.fft.fftshift(tau),pwr(np.fft.fftshift(np.abs(_d))),linewidth=2.5,label='gammaopen')
    	plt.plot(np.fft.fftshift(tau),pwr(np.fft.fftshift(np.abs(_dH))),linewidth=2.5,label='gamma2')
    	
    	#filtered_x = np.array([np.where((np.fft.fftshift(tau) <= -90) & (np.fft.fftshift(tau) >= -170))])
 #       plt.plot((np.fft.fftshift(tau))[filtered_x], pwr(np.fft.fftshift(np.abs(_dopen)))[filtered_x], linewidth=2.0)
        #print len(plt.plot((np.fft.fftshift(tau))[filtered_x]))
    	#plt.plot(tau,pwr(_dH),linewidth=2.5,label='FEED transmission(mag)')
    	
    	#plt.plot(fq,ph,linewidth=2.5,label='Dish+feed return loss in transmission mode(ph)')
    	#plt.plot(fq,ph_f,linewidth=2.5,label='Feed return loss in transmission mode(ph)')
    	#plt.plot(fq,180*np.arctan(np.imag(d2)/np.real(d2))/np.pi,linewidth=2.5,label='Irec/Isky')
    	#plt.plot(fq,pwr(d),linewidth=2.5,label='Dish+feed return loss in transmission mode(mag)')
    	#plt.plot(fq,pwr(d_f),linewidth=2.5,label='Feed return loss in transmission mode(mag)')
    	#plt.plot(fq,pwr(d2),linewidth=2.5,label='Irec/Isky')
    	#plt.plot(fq,pwr(d2),linewidth=2.5,label='FEED transmission in receving mode(mag)')
    	
    	
    	
    
plt.xlim(-350,350)
#plt.ylim(-60,0)
#plt.xlabel('Delay (nS)')
#plt.ylabel('Delay Spectrum')
plt.xlabel('Frequency (GHz)')
plt.ylabel('Return loss magnitude (dB)')
#plt.ylabel('Return loss phase (deg)')
plt.grid()
plt.legend(loc='lower right')
plt.show()

#fig = plt.figure()
#fig.set_size_inches(8.5,8.5) 
#ax = fig.add_subplot(1,1,1) 
#ax.plot(fq,pwr(d),linewidth=2.5,label='Dish+feed return loss in transmission mode(mag)')
#ax.plot(fq,pwr(d_f),linewidth=2.5,label='Feed return loss in transmission mode(mag)')
#ax.plot(fq,pwr(d2),linewidth=2.5,label='Irec/Isky')
#ax.set_xlabel("Frequency(GHz)")
#ax.set_ylabel("Return loss magnitude (dB)")
#fig.savefig("RL_mag1.jpeg",dpi=100, bbox_inches='tight')


	







