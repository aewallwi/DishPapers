#! /usr/bin/env python

import numpy as np, pylab as plt, aipy as a
import sys, csv


def fromcsv(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = np.array(list(d)[18:-1], dtype=np.float)
    #print d[0]
    return x[:,0]/1e3, x[:,1]

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

    WINDOW = 'blackman-harris'
    #WINDOW = 'hamming'
     
    fq,db = fromcsv('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/Simulation/simulation_HERAs11_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq,ph = fromcsv('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/Simulation/simulation_HERAs11_P.csv')
    
    fq_f,db_f = fromcsv('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/Simulation/simulation_feeds11_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_f,ph_f = fromcsv('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/Simulation/simulation_feeds11_P.csv')
    
    
    fq_feed,db_feed = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_FEED_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_feed,ph_feed = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_FEED_P.csv')
    
    fq_H,db_H = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_H,ph_H = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_P.csv')
    
    delay,delay_spectrum = fromcsv('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/Simulation/simulation_aaron.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    
    
    d = lin(db,ph) 
    d_f = lin(db_f,ph_f) 
    d_feed = lin(db_feed,ph_feed)
    d_H = lin(db_H,ph_H)
    
    
    valid = np.where(np.logical_and(fq_f < .200, fq_f > .100)) # restrict to PAPER band
    valid1 = np.where(np.logical_and(fq_H < .200, fq_H > .100)) # restrict to PAPER band
    
    fq, db, ph, d = fq[valid], db[valid], ph[valid],d[valid] 
    fq_f, db_f, ph_f, d_f = fq_f[valid], db_f[valid], ph_f[valid],d_f[valid]
    
    fq_feed, db_feed, ph_feed, d_feed = fq_feed[valid1], db_feed[valid1], ph_feed[valid1], d_feed[valid1]
    
    fq_H,db_H, ph_H, d_H = fq_H[valid1],db_H[valid1], ph_H[valid1], d_H[valid1]
    
    
    tau = np.fft.fftfreq(fq_f.size, fq_f[1]-fq_f[0])
    tau1 = np.fft.fftfreq(fq_H.size, fq_H[1]-fq_H[0])
    
    window = a.dsp.gen_window(fq_f.size, WINDOW)
    window1 = a.dsp.gen_window(fq_H.size, WINDOW)
   

    if True: # approx account for 1st reflection of sky signal off of feed

        
        d_t =1+d_f # Transmission coefficient
        d1 = ((d-d_f)*d_f/d_t ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2 = d1+d_t # Zero term added for final correction
        
        d_t_H =1+d_feed # Transmission coefficient
        d1_H = ((d_H-d_feed)*d_feed/d_t_H ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2_H = d1_H+d_t_H # Zero term added for final correction
        
        _d2 = np.fft.ifft(np.abs(d2)**2*window)/window.mean()# compensate for window amplitude
        _d2H = np.fft.ifft(np.abs(d2_H)**2*window1) /window1.mean()# compensate for window amplitude
        

        plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_d2))), linewidth=2.5, label= 'HERA (Simulation)')

        d_t_10 =1+d_feed/10 # Transmission coefficient
        d1_H = ((d_H/10-d_feed/10)*(d_feed/10)/d_t_H ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2_H = d1_H+d_t_10 # Zero term added for final correction
        
        _d3H = np.fft.ifft(np.abs(d2_H)**2*window1) /window1.mean()# compensate for window amplitude
   
        plt.plot(np.fft.fftshift(tau1), 0.5*pwr(np.fft.fftshift(np.abs(_d2H))),   linewidth=2.5, label='HERA (Measurement)')
        #ax2.plot(np.fft.fftshift(tau1), 0.5*pwr(np.fft.fftshift(np.abs(_d3H))), linewidth=2.0, color ='0.1', label='HERA element in with 20dB better return loss')
        plt.plot(np.fft.fftshift(tau1), 0.5*pwr(np.fft.fftshift(np.abs(_d3H))), color='0.75',linewidth=2.0, label='20dB improvement')
        
        d_t_10 =1+d_feed/30 # Transmission coefficient
        d1_H = ((d_H/30-d_feed/30)*(d_feed/30)/d_t_H ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2_H = d1_H+d_t_10 # Zero term added for final correction
        
        _d3H = np.fft.ifft(np.abs(d2_H)**2*window1) /window1.mean()# compensate for window amplitude
        
        
plt.xlim(-450,450)
plt.xlabel('Delay (nS)')
plt.ylabel('Delay Spectrum (dB)')
plt.grid()
plt.legend(loc='upper right')
plt.show()


	







