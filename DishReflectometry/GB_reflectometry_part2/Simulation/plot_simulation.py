#! /usr/bin/env python

import numpy as np, pylab as plt, aipy as a
import sys, csv
from matplotlib import colors

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
    
    
    fq_feed,db_feed = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_FEED_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_feed,ph_feed = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_FEED_P.csv')
    
    fq_H,db_H = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_H,ph_H = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_P.csv')
    
    #fq_feed =fq_f
    #db_feed =db_f
    #ph_feed =ph_f
    #fq_H = fq
    #db_H = db
    #ph_H = ph
    
    delay,delay_spectrum = fromcsv('simulation_aaron.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    
    
    
    
    
    d = lin(db,ph) 
    d_f = lin(db_f,ph_f) 
    d_feed = lin(db_feed,ph_feed)
    d_H = lin(db_H,ph_H)
    
    
    valid = np.where(np.logical_and(fq < .200, fq > .10)) # restrict to PAPER band
    valid1 = np.where(np.logical_and(fq_H < .200, fq_H > .100)) # restrict to PAPER band
    
    fq, db, ph, d = fq[valid], db[valid], ph[valid], d[valid]
     
    fq_f, db_f, ph_f, d_f = fq_f[valid], db_f[valid], ph_f[valid],d_f[valid]
    
    fq_feed, db_feed, ph_feed, d_feed = fq_feed[valid1], db_feed[valid1], ph_feed[valid1], d_feed[valid1]
    
    fq_H,db_H, ph_H, d_H = fq_H[valid1],db_H[valid1], ph_H[valid1], d_H[valid1]
    
    
    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    tau1 = np.fft.fftfreq(fq_H.size, fq_H[1]-fq_H[0])
    
    window = a.dsp.gen_window(fq.size, WINDOW)
    window1 = a.dsp.gen_window(fq_H.size, WINDOW)
   

    if True: # approx account for 1st reflection of sky signal off of feed
    
        k_ll = k_parallel(tau)
        #print k_ll
        
        d_t =1+d_f # Transmission coefficient
        d1 = ((d-d_f)*d_f/d_t ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2 = d1+d_t # Zero term added for final correction
        
        d_t_H =1+d_feed # Transmission coefficient
        d1_H = ((d_H-d_feed)*d_feed/d_t_H ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2_H = d1_H+d_t_H # Zero term added for final correction
        
        _d2 = np.fft.ifft(np.abs(d2)**2*window)/window.mean()# compensate for window amplitude
        _d2H = np.fft.ifft(np.abs(d2_H)**2*window1) /window1.mean()# compensate for window amplitude
        
        
        #np.savetxt("delay_spectrum_paper_tau.txt",np.fft.fftshift(tau1))
        
    
        
        #plt.plot(fq_f,180*np.arctan2(np.imag(d_f),np.real(d_f))/np.pi,linewidth=2.5,label='return loss HERA element (simulation)' )
        #plt.plot(fq_feed,ph_feed,linewidth=2.5,label='return loss HERA element (measured)' )
        
        #plt.plot(fq,db,linewidth=2.5,label='return loss feed (simulation)' )
        #plt.show()
        #plt.plot(fq_H,pwr(d2_H),linewidth=2.5,label='return loss feed (measured)' )
        
        #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_d2))), linewidth=2.5, color = '#ff8c00', label= 'EM Simulation')
        #plt.plot(1000*delay, delay_spectrum, linewidth=2.5, color = 'black', label='Time domain simulation')
        plt.plot(np.fft.fftshift(tau1), 10.0*np.log10(np.fft.fftshift(np.abs(_d2H))), linewidth=2.5, color = 'blue' ,  label= 'Measurement')
        
        
        #plt.plot(np.fft.fftshift(tau), 0.5*pwr(np.fft.fftshift(np.abs(_d2))), linewidth=2.5, label='HERA feed')
        
        #ax2.set_ylabel('exp', color='b')
        #for tl in ax2.get_yticklabels():
    	#	tl.set_color('b')
    	
    	
    	
        
        
        #plt.plot(1000*delay, delay_spectrum, linewidth=2.5, label='Time domain simulation')
        
        #_dww_rb = np.fft.ifft((d2_rb)**2*window2) / window2.mean() # compensate for window amplitude
    
        d_t_10 =1+d_feed/10 # Transmission coefficient
        d1_H = ((d_H/10-d_feed/10)*(d_feed/10)/d_t_H ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2_H = d1_H+d_t_10 # Zero term added for final correction
        
        _d3H = np.fft.ifft(np.abs(d2_H)**2*window1) /window1.mean()# compensate for window amplitude
    #import IPython; IPython.embed()
    
    
    
        npzdata = np.load('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/data/spec_on_foreground_reflected_power_21cmfast_14.6m_150.0_MHz_subband_v2.npz')
#There will be 4 sets of values under these keys: tau, achrmbeam, chrmbeam, funcbeam
        tau2 = npzdata['tau'] # 1000 values
        achrmbeam_spec = npzdata['achrmbeam'] # 3x1000 array
        chrmbeam_spec = npzdata['chrmbeam'] # 3x1000 array
        funcbeam_spec = npzdata['funcbeam'] # 3x1000 array

        #valid2 = np.where(np.logical_and(k_ll <0.251530831794, k_ll > 0.00)) # restrict to PAPER band
        
        k_ll = k_parallel(tau2)
        #k_ll = k_ll[valid2] 
        #print k_ll
        
        #fig, ax1 = plt.subplots()
        
        #ax1.set_xlabel('k||')
        #ax1.set_ylabel('delay spectrum (dB)')
        
        #ax2 = ax1.twiny()
        #ax2.set_xlabel('delay (ns)')
        
        
        #ax1.plot(k_ll, -achrmbeam_spec[0], linewidth=1.5, color ='k', label='k||>0.1 ')
        plt.plot(tau2, -achrmbeam_spec[0], linewidth=1.5, color ='r', label='FG Simulation k||>0.1 ')
        #ax1.plot(k_ll, -achrmbeam_spec[1], linewidth=1.5, color ='b', label='k||>0.15')
        plt.plot(tau2, -achrmbeam_spec[1], linewidth=1.5, color ='b', label='FG Simulation k||>0.15')
        #ax1.plot(k_ll, -achrmbeam_spec[2], linewidth=1.5, color ='r', label='k||>0.2')
        plt.plot(tau2, -achrmbeam_spec[2], linewidth=1.5, color ='k', label='FG Simulation k||>0.2')
        
        #ax1.fill_between(k_ll, -achrmbeam_spec[0],color='0.5', label=' ')
        plt.fill_between(tau2, -achrmbeam_spec[0],color='0.5', label=' ')
        #ax1.fill_between(k_ll, -achrmbeam_spec[1],color='0.75', label=' ')
        plt.fill_between(tau2, -achrmbeam_spec[1],color='0.75', label=' ')
        #ax1.fill_between(k_ll, -achrmbeam_spec[2],color='0.9', label=' ')
        plt.fill_between(tau2, -achrmbeam_spec[2],color='0.82', label=' ')
        
        #ax2.plot(np.fft.fftshift(tau1), 0.5*pwr(np.fft.fftshift(np.abs(_d2H))), color='m', linewidth=4.5, label='HERA element (Measured)')
        #plt.plot(np.fft.fftshift(tau1), 0.5*pwr(np.fft.fftshift(np.abs(_d2H))),   linewidth=2.5, label='HERA (Measurement)')
        #ax2.plot(np.fft.fftshift(tau1), 0.5*pwr(np.fft.fftshift(np.abs(_d3H))), linewidth=2.0, color ='0.1', label='HERA element in with 20dB better return loss')
        #plt.plot(np.fft.fftshift(tau1), 0.5*pwr(np.fft.fftshift(np.abs(_d3H))), color='c',linewidth=2.0, label='20dB improvement')
        
        
        npzdata2 = np.load('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/Simulation/delayplot/delay_spectrum.npz')
        tau3 = npzdata2['tau_hera']
        ddave = npzdata2['dhera']
        #plt.plot(tau3, -ddave,linewidth=1.5 , label='Dave simulation')
        
        d_t_10 =1+d_feed/10 # Transmission coefficient
        d1_H = ((d_H/10-d_feed/10)*(d_feed/10)/d_t_H ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2_H = d1_H+d_t_10 # Zero term added for final correction
        
        _d3H = np.fft.ifft(np.abs(d2_H)**2*window1) /window1.mean()# compensate for window amplitude
        
        #plt.plot(np.fft.fftshift(tau1), 10.0*np.log10(np.fft.fftshift(np.abs(_d3H))), linewidth=2.5, color = 'gray' ,  label= '')
        
        np.savez('delayspectrum', delay = np.fft.fftshift(tau1), amp = 10.0*np.log10(np.fft.fftshift(np.abs(_d2H))) )
        
        npzdata3 = np.load('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/refl_cstr.npz')
        tau4 = npzdata3['tau_cstr']
        tau5 = npzdata3['tau_meas']
       
        meas_aar = npzdata3['refl_meas']
        meas_aar_010 = npzdata3['refl_cstr_dspec_k0.10']
        meas_aar_015 = npzdata3['refl_cstr_dspec_k0.15']
        meas_aar_020 = npzdata3['refl_cstr_dspec_k0.20']
        
        meas_aar_oqe_010 = npzdata3['refl_cstr_oqe_k0.10']
        meas_aar_oqe_015 = npzdata3['refl_cstr_oqe_k0.15']
        meas_aar_oqe_020 = npzdata3['refl_cstr_oqe_k0.20']
        
        #plt.plot(tau5, meas_aar, linewidth=1.5, color ='k', label='FG Simulation k||>0.2')
        #plt.plot(tau4, meas_aar_010, linewidth=1.5, color ='r', label='FG Simulation k||>0.1')
        #plt.plot(tau4, meas_aar_015, linewidth=1.5, color ='b', label='FG Simulation k||>0.15')
        #plt.plot(tau4, meas_aar_020, linewidth=1.5, color ='k', label='FG Simulation k||>0.2')
        
        plt.plot(tau4, meas_aar_oqe_010, linewidth=1.5, color ='m', label='FG Simulation OQE k||>0.1')
        plt.plot(tau4, meas_aar_oqe_015, linewidth=1.5, color ='g', label='FG Simulation OQE k||>0.15')
        plt.plot(tau4, meas_aar_oqe_020, linewidth=1.5, color ='c', label='FG Simulation OQE k||>0.2')
        #plt.plot(np.fft.fftshift(tau1), 0.5*pwr(np.fft.fftshift(np.abs(_d3H))), linewidth=2.0, color ='m', label='30dB improvement')
        
#-----------------plotting returnloss magnitude--------------
#plt.plot(fq, 10.0*np.log10((np.abs(d)**2)), label='Feed on dish')
#ax1.set_xlim(0,max(k_ll))
#ax2.set_xlim(0,500)
plt.xlim(-500,500)
plt.ylim(-60,0)

#tau2,achrmbeam_spec[0],achrmbeam_spec[1],achrmbeam_spec[2] = tau2[valid2],achrmbeam_spec[0][valid2],achrmbeam_spec[1][valid2],achrmbeam_spec[2][valid2]
plt.xlabel('Delay (nS)')
plt.ylabel('Delay Spectrum (dB)')
#plt.xlabel('Frequency (GHz)')
#plt.ylabel('Return loss magnitude (dB)')
#plt.ylabel('Return loss phase (deg)')
plt.grid()
plt.legend(loc='upper left')
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
	







