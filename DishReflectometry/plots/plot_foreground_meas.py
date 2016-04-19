#! /usr/bin/env python

import numpy as np, pylab as plt, aipy as a
import sys, csv


#def fromcsv(filename):
 #   print 'Reading', filename
  #  d = csv.reader(open(filename,'r'), delimiter=',')
   # x = np.array(list(d)[18:-1], dtype=np.float)
    #print d[0]
    #return x[:,0]/1e3, x[:,1]

def fromcsv1(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = np.array(list(d)[18:-1], dtype=np.float)
    #print d[0]
    return x[:,0]/1e9, x[:,1]
    
def lin(db,ph): # Converting the magnitude and phase of the measurements into complex returnloss in voltage ratio
    return 10**(db/20) * np.exp(2j*np.pi*ph/360)

def pwr(d): # Computing the absolute magnitude of power in db. 
    return 10*np.log10(np.abs(d)**2)
def pwr1(d): # Computing the absolute magnitude of power in db. 
    return 10*np.log10(np.abs(d)**1)
    
def k_parallel(tau):
    f21 = 1.42 #In GHz
    z = (1.42/0.15 -1) 
    Omega_M = 0.27
    Omega_L = 0.73
    Omega_K = 1-Omega_M-Omega_L
    H0 = 100*1000 #m per sec h per Mpc
    
    E = (Omega_M*(1+z)**3+Omega_K*(1+z)**2+Omega_L)**0.5
    
    k_ll = 2*np.pi*tau*f21*H0*E/((3*10**8)*(1+z)**2)    
    
    return k_ll
if True:
#for filename in sys.argv[1:]:
 #   BASE = filename[:-len('.csv')]
  #  db_file = BASE + '_DB.csv'
   # ph_file = BASE + '_P.csv'

    WINDOW = 'blackman-harris'
    #WINDOW = 'hamming'
    
    #fq,db = fromcsv(db_file) # Reading the magnitude and the phase of the datafile to be processed
    #fq,ph = fromcsv(ph_file) 
    
    
    fq_feed,db_feed = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_FEED_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_feed,ph_feed = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_FEED_P.csv')
    
    fq_H,db_H = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_H,ph_H = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_P.csv')
    
    
    
    
    
    
    #d = lin(db,ph)
    d_feed = lin(db_feed,ph_feed)
    d_H = lin(db_H,ph_H)
    
    
    
    valid = np.where(np.logical_and(fq_H < .200, fq_H > .100)) # restrict to PAPER band
    
    #fq, db, ph, d = fq[valid], db[valid], ph[valid], d[valid]
    
    fq_feed, db_feed, ph_feed, d_feed = fq_feed[valid], db_feed[valid], ph_feed[valid], d_feed[valid]
    
    fq_H,db_H, ph_H, d_H = fq_H[valid],db_H[valid], ph_H[valid], d_H[valid]
    
    
    tau = np.fft.fftfreq(fq_feed.size, fq_feed[1]-fq_feed[0])
    
    
    window = a.dsp.gen_window(fq_feed.size, WINDOW)
   

    if True: # approx account for 1st reflection of sky signal off of feed
    
        k_ll = k_parallel(tau)
        #print k_ll
        
        
        d_t_H =1+d_feed # Transmission coefficient
        d1_H = ((d_H-d_feed)*d_feed/d_t_H ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2_H = d1_H+d_t_H # Zero term added for final correction
        
        _d2H = np.fft.ifft(np.abs(d2_H)**2*window) /window.mean()# compensate for window amplitude
        
        
        #np.savetxt("delay_spectrum_paper_tau.txt",np.fft.fftshift(tau1))
        
    
        
        #plt.plot(fq_f,180*np.arctan2(np.imag(d_f),np.real(d_f))/np.pi,linewidth=2.5,label='return loss HERA element (simulation)' )
        #plt.plot(fq_feed,ph_feed,linewidth=2.5,label='return loss HERA element (measured)' )
        
        #plt.plot(fq,pwr(d2),linewidth=2.5,label='return loss feed (simulation)' )
        #plt.plot(fq_H,pwr(d2_H),linewidth=2.5,label='return loss feed (measured)' )
        
        #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_d))), linewidth=2.5, label= 'Feed Only (Simulation)')
        #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_d2))), linewidth=2.5, label= 'HERA (Simulation)')
        
        #plt.plot(np.fft.fftshift(tau), 0.5*pwr(np.fft.fftshift(np.abs(_d2))), linewidth=2.5, label='HERA feed')
        
        #ax2.set_ylabel('exp', color='b')
        #for tl in ax2.get_yticklabels():
    	#	tl.set_color('b')
    	
    	
    	
        
        
        #plt.plot(1000*delay, delay_spectrum, linewidth=2.5, label='Simulated power kernel of HERA element')
        
        #_dww_rb = np.fft.ifft((d2_rb)**2*window2) / window2.mean() # compensate for window amplitude
    
        d_t_10 =1+d_feed/10 # Transmission coefficient
        d1_H = ((d_H/10-d_feed/10)*(d_feed/10)/d_t_H ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2_H = d1_H+d_t_10 # Zero term added for final correction
        
        _d3H = np.fft.ifft(np.abs(d2_H)**2*window) /window.mean()# compensate for window amplitude
    #import IPython; IPython.embed()
    
    
    
        npzdata = np.load('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/data/spec_on_foreground_reflected_power_21cmfast_14.6m_150.0_MHz_subband_v2.npz')
#There will be 4 sets of values under these keys: tau, achrmbeam, chrmbeam, funcbeam
        tau2 = npzdata['tau'] # 1000 values
        achrmbeam_spec = npzdata['achrmbeam'] # 3x1000 array
        chrmbeam_spec = npzdata['chrmbeam'] # 3x1000 array
        funcbeam_spec = npzdata['funcbeam'] # 3x1000 array
        
        
        #ax1.plot(k_ll, -achrmbeam_spec[0], linewidth=1.5, color ='k', label='k||>0.1 ')
        plt.plot(tau2, -achrmbeam_spec[0], linewidth=1.5, color ='r', label='k||>0.1 ')
        #ax1.plot(k_ll, -achrmbeam_spec[1], linewidth=1.5, color ='b', label='k||>0.15')
        plt.plot(tau2, -achrmbeam_spec[1], linewidth=1.5, color ='b', label='k||>0.15')
        #ax1.plot(k_ll, -achrmbeam_spec[2], linewidth=1.5, color ='r', label='k||>0.2')
        plt.plot(tau2, -achrmbeam_spec[2], linewidth=1.5, color ='k', label='k||>0.2')
        
        #ax1.fill_between(k_ll, -achrmbeam_spec[0],color='0.5', label=' ')
        plt.fill_between(tau2, -achrmbeam_spec[0],color='0.5', label=' ')
        #ax1.fill_between(k_ll, -achrmbeam_spec[1],color='0.75', label=' ')
        plt.fill_between(tau2, -achrmbeam_spec[1],color='0.75', label=' ')
        #ax1.fill_between(k_ll, -achrmbeam_spec[2],color='0.9', label=' ')
        plt.fill_between(tau2, -achrmbeam_spec[2],color='0.9', label=' ')
        
        #ax2.plot(np.fft.fftshift(tau1), 0.5*pwr(np.fft.fftshift(np.abs(_d2H))), color='m', linewidth=4.5, label='HERA element (Measured)')
        plt.plot(np.fft.fftshift(tau), 0.5*pwr(np.fft.fftshift(np.abs(_d2H))),  color='m', linewidth=2.5, label='HERA (Measurement)')
        #ax2.plot(np.fft.fftshift(tau1), 0.5*pwr(np.fft.fftshift(np.abs(_d3H))), linewidth=2.0, color ='0.1', label='HERA element in with 20dB better return loss')
        plt.plot(np.fft.fftshift(tau), 0.5*pwr(np.fft.fftshift(np.abs(_d3H))), color='0.75',linewidth=2.0, label='20dB improvement')
        
        d_t_10 =1+d_feed/30 # Transmission coefficient
        d1_H = ((d_H/30-d_feed/30)*(d_feed/30)/d_t_H ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2_H = d1_H+d_t_10 # Zero term added for final correction
        
        _d3H = np.fft.ifft(np.abs(d2_H)**2*window) /window.mean()# compensate for window amplitude
        
        
        plt.plot(np.fft.fftshift(tau), 0.5*pwr(np.fft.fftshift(np.abs(_d3H))), linewidth=2.0, color ='c', label='30dB improvement')
        
#-----------------plotting returnloss magnitude--------------
#plt.plot(fq, 10.0*np.log10((np.abs(d)**2)), label='Feed on dish')
#ax1.set_xlim(0,max(k_ll))
#ax2.set_xlim(0,500)
plt.xlim(-0,450)
#plt.ylim(-60,0)

#tau2,achrmbeam_spec[0],achrmbeam_spec[1],achrmbeam_spec[2] = tau2[valid2],achrmbeam_spec[0][valid2],achrmbeam_spec[1][valid2],achrmbeam_spec[2][valid2]
plt.xlabel('Delay (nS)')
plt.ylabel('Delay Spectrum (dB)')
#plt.xlabel('Frequency (GHz)')
#plt.ylabel('Return loss magnitude (dB)')
#plt.ylabel('Return loss phase (deg)')
plt.grid()
plt.legend(loc='upper right')
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
	







