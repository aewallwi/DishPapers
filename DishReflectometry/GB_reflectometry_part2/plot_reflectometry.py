#! /usr/bin/env python

import numpy as np, pylab as plt, aipy as a
import sys, csv

def fromcsv2(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = np.array(list(d)[18:-1], dtype=np.float)
    return x[:,0]/1e3, x[:,1]
def fromcsv1(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = np.array(list(d)[18:-1], dtype=np.float)
    return x[:,0]/1e3, x[:,1]
def fromcsv(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = np.array(list(d)[18:-1], dtype=np.float)
    return x[:,0]/1e9, x[:,1]

for filename in sys.argv[1:]:
    BASE = filename[:-len('.csv')]
    db_file = BASE + '_DB.csv'
    ph_file = BASE + '_P.csv'

    WINDOW = 'blackman-harris'
    #WINDOW = 'hamming'
    
    fq,db = fromcsv(db_file) # Reading the magnitude and the phase of the datafile to be processed
    fq,ph = fromcsv(ph_file) 
    
    fq_f,db_f = fromcsv('HERA_FEED_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_f,ph_f = fromcsv('HERA_FEED_P.csv')
    
    d = 10**(db/20) * np.exp(2j*np.pi*ph/360) # 20 to put into voltage amplitude, not power
    d_f = 10**(db_f/20) * np.exp(2j*np.pi*ph_f/360) # 20 to put into voltage amplitude, not power
    
    fq_rb,db_rb = fromcsv1('measurement_dishAndFeed_DB.csv')# Reading Rich's dish+feed measurements magnitude and phase file.
    fq_rb,ph_rb = fromcsv1('measurement_dishAndFeed_P.csv')
    
    fq_rb_f,db_rb_f = fromcsv1('measurement_FeedOnly_DB.csv')#Reading Rich's feed only file to do the zero delay corrections.
    fq_rb_f,ph_rb_f = fromcsv1('measurement_FeedOnly_P.csv')
    
    #fq_rb_f,db_rb_f = fromcsv('Calibration/OPEN_DB.csv')#Reading Rich's feed only file to do the zero delay corrections.
    #fq_rb_f,ph_rb_f = fromcsv('Calibration/OPEN_P.csv')
    
    d_rb = 10**(db_rb/20) * np.exp(2j*np.pi*ph_rb/360) # 20 to put into voltage amplitude, not power 
    d_rb_f = 10**(db_rb_f/20) * np.exp(2j*np.pi*ph_rb_f/360) # 20 to put into voltage amplitude, not power 
    

    valid = np.where(np.logical_and(fq < .200, fq > .100)) # restrict to PAPER band
    valid1 = np.where(np.logical_and(fq_rb_f < .200, fq_rb_f > .100)) # restrict to PAPER band
    
    fq, d, db, ph, fq_f,d_f, db_f, ph_f = fq[valid], d[valid], db[valid], ph[valid],fq_f[valid], d_f[valid], db_f[valid],ph_f[valid]
    
    fq_rb, d_rb, db_rb, ph_rb, fq_rb_f,d_rb_f,db_rb_f,ph_rb_f = fq_rb[valid1], d_rb[valid1], db_rb[valid1], ph_rb[valid1], fq_rb_f[valid1],d_rb_f[valid1],db_rb_f[valid1], ph_rb_f[valid1]
    #elif v == '130-160 MHz':
     # fq, d, db, ph = fq[valids[v]], d[valids[v]], db[valids[v]], ph[valids[v]]
      #print np.abs(d[0])
    #elif v == '160-190 MHz':    
     # fq, d, db, ph = fq[valids[v]], d[valids[v]], db[valids[v]], ph[valids[v]]
    
    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    tau1 = np.fft.fftfreq(fq_rb.size, fq_rb[1]-fq_rb[0])
    window = a.dsp.gen_window(fq.size, WINDOW)
    window2 = a.dsp.gen_window(fq_rb_f.size, WINDOW)
   

    if True: # approx account for 1st reflection of sky signal off of feed
        #_dw *= np.abs(_d[0])
        #_d *= np.abs(_d[0])
        
        #d -= np.abs(_d[0])
        #d1 = (np.abs(d)-np.abs(d_f))*np.abs(d_f)/(1-np.abs(d_f))+ (1-np.abs(d_f))
        #d_rb_t = 1+d_rb_f # Transmission coefficient
        #d1_rb = ((d_rb-d_rb_f)*d_rb_f/d_rb_t ) # Reflections corrected to produce system bandpass in the receiving mode
        #d2_rb = d1_rb+d_rb_t # Zero term added for final correction
        
        _dw1 = np.fft.ifft((d_rb_f)**2*window2) /window2.mean()
        
        d_t =1+d_f # Transmission coefficient
        d1 = ((d-d_f)*d_f/d_t ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2 = d1+d_t # Zero term added for final correction

        _d = np.fft.ifft((d2)**2*window)/window.mean()
        _dww = np.fft.ifft(np.abs(d2)**2*window) /window.mean()# compensate for window amplitude
        #print np.fft.fftshift(tau)
        
        #_dww_rb = np.fft.ifft((d2_rb)**2*window2) / window2.mean() # compensate for window amplitude

    #import IPython; IPython.embed()
    #plt.plot(fq_f,np.abs(d_f),linewidth=2.5,label='Transmission UCB(real)')
    #plt.plot(fq_f,np.abs(d2),linewidth=2.5,label='Transmission UCB(real)')
    #plt.show()
    #plt.plot(fq,d_t.imag,linewidth=2.5,label='Transmission UCB (imag)')
    #plt.plot(fq,10.0*np.log10(np.abs(d2)),linewidth=2.5,label='Corrected measurements UCB (abs mag)')
    #plt.plot(fq_rb_f,d_rb_t.real,linewidth=2.5,label='Transmission RB (real)')
    #plt.plot(fq_rb_f,d_rb_t.imag,linewidth=2.5,label='Transmission RB (imag)')
    #plt.plot(fq_rb_f,10.0*np.log10(np.abs(d2_rb)),linewidth=2.5,label='Corrected measurements RB (abs mag)')
    
    
   #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_d))), linewidth=2.5, label= 'Square Window')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_d))), linewidth=2.5, label='Feed only (measurement)')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dww))), linewidth=2.5, label='current balun')
    
    
    db = db-30.0
    d = 10**(db/20) * np.exp(2j*np.pi*ph/360) # 20 to put into voltage amplitude, not power
    db_f = db_f-30.0
    d_f = 10**(db_f/20) * np.exp(2j*np.pi*ph_f/360) # 20 to put into voltage amplitude, not power
    
    fq_rb,db_rb = fromcsv1('measurement_dishAndFeed_DB.csv')# Reading Rich's dish+feed measurements magnitude and phase file.
    fq_rb,ph_rb = fromcsv1('measurement_dishAndFeed_P.csv')
    
    fq_rb_f,db_rb_f = fromcsv1('measurement_FeedOnly_DB.csv')#Reading Rich's feed only file to do the zero delay corrections.
    fq_rb_f,ph_rb_f = fromcsv1('measurement_FeedOnly_P.csv')
    
    #fq_rb_f,db_rb_f = fromcsv('Calibration/OPEN_DB.csv')#Reading Rich's feed only file to do the zero delay corrections.
    #fq_rb_f,ph_rb_f = fromcsv('Calibration/OPEN_P.csv')
    
    d_rb = 10**(db_rb/20) * np.exp(2j*np.pi*ph_rb/360) # 20 to put into voltage amplitude, not power 
    d_rb_f = 10**(db_rb_f/20) * np.exp(2j*np.pi*ph_rb_f/360) # 20 to put into voltage amplitude, not power 
    

    valid = np.where(np.logical_and(fq < .200, fq > .100)) # restrict to PAPER band
    valid1 = np.where(np.logical_and(fq_rb_f < .200, fq_rb_f > .100)) # restrict to PAPER band
    
    fq, d, db, ph, fq_f,d_f, db_f, ph_f = fq[valid], d[valid], db[valid], ph[valid],fq_f[valid], d_f[valid], db_f[valid],ph_f[valid]
    
    fq_rb, d_rb, db_rb, ph_rb, fq_rb_f,d_rb_f,db_rb_f,ph_rb_f = fq_rb[valid1], d_rb[valid1], db_rb[valid1], ph_rb[valid1], fq_rb_f[valid1],d_rb_f[valid1],db_rb_f[valid1], ph_rb_f[valid1]
    #elif v == '130-160 MHz':
     # fq, d, db, ph = fq[valids[v]], d[valids[v]], db[valids[v]], ph[valids[v]]
      #print np.abs(d[0])
    #elif v == '160-190 MHz':    
     # fq, d, db, ph = fq[valids[v]], d[valids[v]], db[valids[v]], ph[valids[v]]
    
    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    tau1 = np.fft.fftfreq(fq_rb.size, fq_rb[1]-fq_rb[0])
    window = a.dsp.gen_window(fq.size, WINDOW)
    window2 = a.dsp.gen_window(fq_rb_f.size, WINDOW)
   

    if True: # approx account for 1st reflection of sky signal off of feed
        #_dw *= np.abs(_d[0])
        #_d *= np.abs(_d[0])
        
        #d -= np.abs(_d[0])
        #d1 = (np.abs(d)-np.abs(d_f))*np.abs(d_f)/(1-np.abs(d_f))+ (1-np.abs(d_f))
        #d_rb_t = 1+d_rb_f # Transmission coefficient
        #d1_rb = ((d_rb-d_rb_f)*d_rb_f/d_rb_t ) # Reflections corrected to produce system bandpass in the receiving mode
        #d2_rb = d1_rb+d_rb_t # Zero term added for final correction
        
        _dw1 = np.fft.ifft((d_rb_f)**2*window2) /window2.mean()
        
        d_t =1+d_f # Transmission coefficient
        d1 = ((d-d_f)*d_f/d_t ) # Reflections corrected to produce system bandpass in the receiving mode 
        d2 = d1+d_t # Zero term added for final correction

        _d = np.fft.ifft(np.abs(d_f)**2*window)/window.mean()
        _dww = np.fft.ifft(np.abs(d2)**2*window) /window.mean()# compensate for window amplitude
        #print np.fft.fftshift(tau)
        
        #_dww_rb = np.fft.ifft((d2_rb)**2*window2) / window2.mean() # compensate for window amplitude

    #import IPython; IPython.embed()
    #plt.plot(fq_rb_f,(np.abs(d_rb_t))**2,linewidth=2.5,label='Transmission UCB(real)')
    #plt.plot(fq,d_t.imag,linewidth=2.5,label='Transmission UCB (imag)')
    #plt.plot(fq,10.0*np.log10(np.abs(d2)),linewidth=2.5,label='Corrected measurements UCB (abs mag)')
    #plt.plot(fq_rb_f,d_rb_t.real,linewidth=2.5,label='Transmission RB (real)')
    #plt.plot(fq_rb_f,d_rb_t.imag,linewidth=2.5,label='Transmission RB (imag)')
    #plt.plot(fq_rb_f,10.0*np.log10(np.abs(d2_rb)),linewidth=2.5,label='Corrected measurements RB (abs mag)')
    
    
   #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_d))), linewidth=2.5, label= 'Square Window')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_d))), linewidth=2.5, label='Feed only (measurement)')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dww))), linewidth=2.5, label='Better_balun')
    
    
   

#-----------------plotting returnloss magnitude--------------
#plt.plot(fq, 10.0*np.log10((np.abs(d)**2)), label='Feed on dish')
plt.xlim(-350,350)
#plt.ylim(-60,0)
plt.xlabel('Delay (nS)')
plt.ylabel('Delay Spectrum')
#plt.xlabel('Frequency (GHz)')
#plt.ylabel('Return loss magnitude (dB)')
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
	







