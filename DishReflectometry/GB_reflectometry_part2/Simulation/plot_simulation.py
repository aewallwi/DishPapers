#! /usr/bin/env python

import numpy as np, pylab as plt, aipy as a
import sys, csv


def fromcsv(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = np.array(list(d)[18:-1], dtype=np.float)
    #print d[0]
    return x[:,0]/1e3, x[:,1]

for filename in sys.argv[1:]:
    BASE = filename[:-len('.csv')]
    db_file = BASE + '_DB.csv'
    ph_file = BASE + '_P.csv'

    WINDOW = 'blackman-harris'
    #WINDOW = 'hamming'
    
    fq,db = fromcsv(db_file) # Reading the magnitude and the phase of the datafile to be processed
    fq,ph = fromcsv(ph_file) 
    
    fq_f,db_f = fromcsv('simulation_feeds11_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_f,ph_f = fromcsv('simulation_feeds11_P.csv')
    
    d = 10**(db/20) * np.exp(2j*np.pi*ph/360) # 20 to put into voltage amplitude, not power
    d_f = 10**(db_f/20) * np.exp(2j*np.pi*ph_f/360) # 20 to put into voltage amplitude, not power
    
    valid = np.where(np.logical_and(fq < .200, fq > .090)) # restrict to PAPER band
    
    fq, d, db, ph, fq_f,d_f, db_f, ph_f = fq[valid], d[valid], db[valid], ph[valid],fq_f[valid], d_f[valid], db_f[valid],ph_f[valid]
    
    
    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    window = a.dsp.gen_window(fq.size, WINDOW)
   

    if True: # approx account for 1st reflection of sky signal off of feed
        
        d_t =1+d_f # Transmission coefficient
        d1 = ((d-d_f)*d_f/d_t ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2 = d1+d_t # Zero term added for final correction

        _d = np.fft.ifft((d_f)**2*window)/window.mean()
        _dw = np.fft.ifft((d)**2*window)/window.mean()
        _dww = np.fft.ifft((d2)**2*window) /window.mean()# compensate for window amplitude
        
        #plt.plot(fq_f,10.0*np.log10(np.abs((d_f))**2),linewidth=2.5,label='return loss feed (simulation)' )
        plt.plot(fq_f,10.0*np.log10(np.abs((d))**2),linewidth=2.5,label='return loss HERA element (simulations)' )
        
        #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_d))), linewidth=2.5, label= 'Feed Only (Simulation)')
        #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dw))), linewidth=2.5, label= 'HERA element in transmission (Simulation)')
        #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dww))), linewidth=2.5, label='HERA element in reception (Simulation)')
        #_dww_rb = np.fft.ifft((d2_rb)**2*window2) / window2.mean() # compensate for window amplitude
    

    #import IPython; IPython.embed()
    
    
    
   

#-----------------plotting returnloss magnitude--------------
#plt.plot(fq, 10.0*np.log10((np.abs(d)**2)), label='Feed on dish')
#plt.xlim(-350,350)
#plt.ylim(-60,0)
#plt.xlabel('Delay (nS)')
#plt.ylabel('Delay Spectrum')
plt.xlabel('Frequency (GHz)')
plt.ylabel('Return loss magnitude (dB)')
plt.grid()
plt.legend(loc='lower right')
plt.show()

#-----------------plotting returnloss phase--------------
#plt.plot(fq, ph, label='Feed on dish')
#plt.xlabel('Frequency (GHz)')
#plt.ylabel('Phase (Degree)')
#plt.grid()
#plt.legend(loc='lower right')
#plt.show()
#-----------------plotting delay spectrum--------------
#plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dw))), linewidth=2.5, label=BASE.replace('_', ' '))
#plt.xlim(-350,350)
#plt.xlabel('Delay(ns)')
#plt.ylabel('Delay spectrum (dB)')
#plt.grid()
#plt.legend(loc='lower right')
#plt.show()
	







