#! /usr/bin/env python

import numpy as np, pylab as plt, aipy as a, scipy as sp
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
    
def k_parallel(tau):
    f21 = 1.42 #In GHz
    z = (1.42/0.15 -1)
    Omega_M = 0.27
    Omega_L = 0.73
    Omega_K = 1-Omega_M-Omega_L
    H0 = 100*1000 #m per sec h per Mpc
    E = (Omega_M*(1+z)**3+Omega_K*(1+z)**2+Omega_L)**0.5
    k_ll = 2*np.pi*tau*H0*E/(3*10**8)*(1+z)**2    
    return k_ll
    
def multiy(x,y1,y2):
    fig, ax1 = plt.subplots()
    ax1.plot(x,y1, color = 'k', linewidth=2.5)
    ax1.set_xlabel('Freq (GHz)')
    ax1.set_ylabel('Magnitude (dB)')
    ax2 = ax1.twinx()
    ax2.plot(x,y2,linewidth=2.5)
    ax2.set_ylabel('Phase (Deg)')
    ax2.set_ylabel('exp', color='b')
    for tl in ax2.get_yticklabels():
    	tl.set_color('b')
    

for filename in sys.argv[1:]:
    BASE = filename[:-len('.csv')]
    db_file = BASE + '_DB.csv'
    ph_file = BASE + '_P.csv'
    
    fq,db = fromcsv(db_file) # Reading the magnitude and the phase of the datafile to be processed
    fq,ph = fromcsv(ph_file) 
    
    fq_f,db_f = fromcsv('HERA_FEED_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_f,ph_f = fromcsv('HERA_FEED_P.csv')
    
    
    fq_p,db_p = fromcsv('PAPER1_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_p,ph_p = fromcsv('PAPER1_P.csv')
    
    fq_open,db_open = fromcsv('Calibration/OPENNN_DB.csv')# Reading the open file to calibrate for the multiple signal reflections in the 50 feet long cable
    fq_open,ph_open = fromcsv('Calibration/OPENNN_P.csv')
    
    fq_short,db_short = fromcsv('Calibration/SHORTT_DB.csv')# Reading the open file to calibrate for the multiple signal reflections in the 50 feet long cable
    fq_short,ph_short = fromcsv('Calibration/SHORTT_P.csv')
    

        
    d = lin(db,ph) 
    d_f = lin(db_f,ph_f) 
    d_open = lin(db_open,ph_open)  
    d_short = lin(db_short,ph_short)
    d_p = lin(db_p, ph_p)
    
    
    
    valid = np.where(np.logical_and(fq < .175, fq > .140)) # restrict to PAPER band
    
    valid1 = np.where(np.logical_and(fq < .175, fq > .140)) # restrict to PAPER band
    
    fq, d, db, ph = fq[valid], d[valid], db[valid], ph[valid]
    
    fq_f,d_f, db_f, ph_f = fq_f[valid], d_f[valid], db_f[valid],ph_f[valid]
    
    fq_open,d_open, db_open, ph_open = fq_open[valid], d_open[valid], db_open[valid],ph_open[valid]
    
    fq_short,d_short, db_short, ph_short = fq_short[valid], d_short[valid], db_short[valid],ph_short[valid]
    
    fq_p,d_p, db_p, ph_p = fq_p[valid], d_p[valid], db_p[valid],ph_p[valid]
    
    
    
    
    

    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])

    
    
    WINDOW = 'blackman-harris'
    window = a.dsp.gen_window(fq.size, WINDOW)
    
    
   

    if True: 
    
        #Computing VNA input return loss from the open data. 
        k_ll = k_parallel(tau) 
        gamma_vna_open = vna_open(db_open,ph_open,fq)
        #fig, ax1 = plt.subplots()
        #ax1.plot(fq,pwr(gamma_vna_open), color = 'k', linewidth=2.5)
        #ax1.set_xlabel('Freq (GHz)')
        #ax1.set_ylabel('Magnitude (dB)')
        #ax2 = ax1.twinx()
        #ax2.plot(fq,180*np.arctan(np.imag(gamma_vna_open)/np.real(gamma_vna_open))/np.pi,linewidth=2.5)
        #ax2.set_ylabel('Phase (Deg)')
        #ax2.set_ylabel('exp', color='b')
        #for tl in ax2.get_yticklabels():
    	#	tl.set_color('b')

        
        
        
        d_t =1+d_f # Computing the transmission coefficient
        
        d1 = ((d-d_f)*d_f/d_t ) # Reflections corrected to produce system bandpass in the receiving mode 
        
        d2 = d1+d_t # Zero term added for final correction
        
        
        
        
        

        _dp = np.fft.ifft(np.abs(d_p)**2*window)/window.mean()
        plt.plot(np.fft.fftshift(tau),0.5*pwr(np.fft.fftshift(np.abs(_dp))),linewidth=2.5,label='PAPER')
        
        _df = np.fft.ifft(np.abs(d_f)**2*window)/window.mean() #Take the Fourier transform of the SQUARE of the absolute MAGNITUDE of the voltage return loss of the feed.
        plt.plot(np.fft.fftshift(tau),0.5*pwr(np.fft.fftshift(np.abs(_df))),linewidth=2.5,label='HERA feed')
        
        _dt = np.fft.ifft(np.abs(d_t)**2*window)/window.mean() #Take the Fourier transform of the SQUARE of the absolute MAGNITUDE of the voltage transmission of the feed.
        #plt.plot(np.fft.fftshift(tau),pwr(np.fft.fftshift(np.abs(_dt))),linewidth=1.5,label='Transmission')
        
        _dopen = np.fft.ifft(np.abs(d_open)**2*window)/window.mean() #Take the Fourier transform of the SQUARE of the absolute magnitude of the voltage return loss of the open.
        #plt.plot(np.fft.fftshift(tau),0.5*pwr(np.fft.fftshift(np.abs(_dopen))),linewidth=3.5,label='Open')
        
        _dvna = np.fft.ifft(np.abs(gamma_vna_open)**2*window)/window.mean() #Take the Fourier transform of the SQUARE of the absolute magnitude of the voltage return loss of the vna input.
        #plt.plot(np.fft.fftshift(tau),pwr(np.fft.fftshift(np.abs(_dt))),linewidth=1.5,label=' VNA input return loss')
        
        _dH = np.fft.ifft(np.abs(d2)**2*window) /window.mean()#Take the Fourier transform of the SQUARE of the absolute magnitude of the voltage return loss  of the HERA element.
        plt.plot(np.fft.fftshift(tau),0.5*pwr(np.fft.fftshift(np.abs(_dH))),linewidth=2.5,label='HERA feed on dish in reception')
        
        
        #_dopen[[np.where((np.fft.fftshift(tau) <= -90) & (np.fft.fftshift(tau) >= -170))]]=0.0
        
    	
    	#plt.plot(np.fft.fftshift(k_ll),pwr(np.fft.fftshift(np.abs(_df))),linewidth=1.5,label='Feed')
    	
    	
    	
    	
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
plt.ylabel('Delay Spectrum(dB)')
plt.xlabel('Frequency (GHz)')
#plt.ylabel('Return loss magnitude (dB)')
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


	







