#! /usr/bin/env python

import numpy as np, pylab as plt, aipy as a
import sys, csv

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
    
    d_rb = 10**(db_rb/20) * np.exp(2j*np.pi*ph_rb/360) # 20 to put into voltage amplitude, not power 
    d_rb_f = 10**(db_rb_f/20) * np.exp(2j*np.pi*ph_rb_f/360) # 20 to put into voltage amplitude, not power 
    
    
    
#    valids = {
#          '100- 130 MHz'  : np.where(np.logical_and(fq>.100 ,fq<.130)), 
#          '130 - 160 MHz' : np.where(np.logical_and(fq>.130 ,fq<.160)), 
#          '160 - 190 MHz' : np.where(np.logical_and(fq>.160 ,fq<.190)),
#          }
          

#for i,v in enumerate(valids.keys()):

    #valid = np.ones(fq.size, dtype=np.bool) # use entire sampled band
    #valid = np.where(fq < .250) # restrict to HERA band
#    if v == '100-130 MHz': 
    valid = np.where(np.logical_and(fq < .200, fq > .100)) # restrict to PAPER band
    
    fq, d, db, ph, d_f, db_f, ph_f = fq[valid], d[valid], db[valid], ph[valid], d_f[valid], db_f[valid],ph_f[valid]
    fq_rb, d_rb, db_rb, ph_rb, fq_rb_f,d_rb_f,db_rb_f,ph_rb_f = fq_rb[valid], d_rb[valid], db_rb[valid], ph_rb[valid], fq_rb_f[valid],d_rb_f[valid],db_rb_f, ph_rb_f[valid]
    #elif v == '130-160 MHz':
     # fq, d, db, ph = fq[valids[v]], d[valids[v]], db[valids[v]], ph[valids[v]]
      #print np.abs(d[0])
    #elif v == '160-190 MHz':    
     # fq, d, db, ph = fq[valids[v]], d[valids[v]], db[valids[v]], ph[valids[v]]
    
    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    window = a.dsp.gen_window(fq.size, WINDOW)
   

    if True: # approx account for 1st reflection of sky signal off of feed
        #_dw *= np.abs(_d[0])
        #_d *= np.abs(_d[0])
        
        #d -= np.abs(_d[0])
        #d1 = (np.abs(d)-np.abs(d_f))*np.abs(d_f)/(1-np.abs(d_f))+ (1-np.abs(d_f))
        d_t_rb = 1+d_rb_f
        d2 = (d_rb-d_rb_f)/d_rb_f+d_t_rb
        d_t =1+d_f
        d1 = ((d-d_f)*d_f/d_t )+ d_t
#        d1 = (1+d_f)
        
        #d1 *= np.abs(_dw[0])/(1-np.abs(_dw[0]))
        _d = np.fft.ifft((np.abs(d1))**2)
        _dww = np.fft.ifft((np.abs(d1))**2*window) /window.mean()# compensate for window amplitude
        _dww_rb = np.fft.ifft((np.abs(d2))**2*window) / window.mean() # compensate for window amplitude
    
    print np.abs(d[0])
    H0=100#h Km/sec/Mpc
    c = 3*10**5# Km/sec
    f21 = 1420#GHz
    Omega_m = 0.27
    Omega_lambda = 0.73
    Omega_k = (1-Omega_m-Omega_lambda)
    z = 8.0
    Ez = Omega_m*(1+z)**3+Omega_k*(1+z)**2+Omega_lambda*(1+z)
    k_ll = (f21*H0*Ez)/(c*(1+z)**2)
    
    #x1 = np.fft.fftshift(tau)
    #x2 = 2*np.pi*x1*k_ll
    #y1 = 10.0*np.log10(np.fft.fftshift(np.abs(_dww)))
    #fig = plt.figure()
    #ax2 = fig.add_subplot(111)
    #ax1 = ax2.twinx()
    
   # line1, = ax2.plot(x1, y1, linewidth=2.5, label=BASE.replace('_', ' '))
   # line2, = ax1.plot(x1, 10.0*np.log10(np.fft.fftshift(np.abs(_dww))), linewidth=2.5, label=BASE.replace('_', ' '))
    
    
   #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_d))), linewidth=2.5, label= 'Square Window')
    plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dww))), linewidth=2.5, label=BASE.replace('_', ' '))
    plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dww_rb))), linewidth=2.5, label='DishandFeed_RB'.replace('_', ' '))
    
    
    #plt.plot(fq, 10.0*np.log10((np.abs(_dww)**2)), label=BASE.replace('_', ' '))
    #plt.plot(fq, db,linewidth=2.5, label=BASE.replace('_', ' '))
    #plt.plot(fq, ph, linewidth=2.5, label=BASE.replace('_', ' '))
   

#-----------------plotting returnloss magnitude--------------
#plt.plot(fq, 10.0*np.log10((np.abs(d)**2)), label='Feed on dish')
#plt.xlim(-350,350)
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
	







