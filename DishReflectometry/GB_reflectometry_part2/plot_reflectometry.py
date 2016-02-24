#! /usr/bin/env python

import numpy as np, pylab as plt, aipy as a
import sys, csv

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

    fq,db = fromcsv(db_file)
    fq,ph = fromcsv(ph_file)
    #d = 10**(db/10) * np.exp(2j*np.pi*ph/360) # power
    d = 10**(db/20) * np.exp(2j*np.pi*ph/360) # 20 to put into voltage amplitude, not power
    valids = {
          '100- 130 MHz'  : np.where(np.logical_and(fq>.100 ,fq<.130)), 
          '130 - 160 MHz' : np.where(np.logical_and(fq>.130 ,fq<.160)), 
          '160 - 190 MHz' : np.where(np.logical_and(fq>.160 ,fq<.190)),
          }
          

for i,v in enumerate(valids.keys()):

    #valid = np.ones(fq.size, dtype=np.bool) # use entire sampled band
    #valid = np.where(fq < .250) # restrict to HERA band
    if v == '100-130 MHz': 
    #valid1 = np.where(np.logical_and(fq < .130, fq > .100)) # restrict to PAPER band
    #valid2 = np.where(np.logical_and(fq < .160, fq > .130)) # restrict to PAPER band
    #valid3 = np.where(np.logical_and(fq < .190, fq > .160))
#if i==0:
      fq, d, db, ph = fq[valid[v]], d[valid[v]], db[valid[v]], ph[valid[v]]
#elif i==1:
 #   fq, d, db, ph = fq[valid2], d[valid2], db[valid2], ph[valid2]
#elif i==2:    
  #  fq, d, db, ph = fq[valid3], d[valid3], db[valid3], ph[valid3]
    
    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    window = a.dsp.gen_window(fq.size, WINDOW)
    
    _d = np.fft.ifft((np.abs(d))**2)
    _dw = np.fft.ifft((np.abs(d))**2*window) / window.mean() # compensate for window amplitude

    if True: # approx account for 1st reflection of sky signal off of feed
        #_dw *= np.abs(_d[0])
        #_d *= np.abs(_d[0])
        
        #d -= np.abs(_d[0])
        d *= np.abs(_dw[0])/(1-np.abs(_dw[0]))
        _d = np.fft.ifft((np.abs(d))**2)
        _dww = np.fft.ifft((np.abs(d))**2*window) /window.mean()# compensate for window amplitude
    
    
    print np.abs(_d[0])
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
    
    #plt.plot(fq, 10.0*np.log10((np.abs(_dww)**2)), label=BASE.replace('_', ' '))
    #plt.plot(fq, db,linewidth=2.5, label=BASE.replace('_', ' '))
    #plt.plot(fq, ph, linewidth=2.5, label=BASE.replace('_', ' '))
   

#-----------------plotting returnloss magnitude--------------
#plt.plot(fq, 10.0*np.log10((np.abs(d)**2)), label='Feed on dish')
plt.xlim(-350,350)
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
	







