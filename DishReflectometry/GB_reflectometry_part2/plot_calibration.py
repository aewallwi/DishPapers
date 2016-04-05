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
    
    fq,db = fromcsv(db_file) # Reading the magnitude and the phase of the datafile to be processed
    fq,ph = fromcsv(ph_file) 
    d = 10**(db/20) * np.exp(2j*np.pi*ph/360) # Converting the data to linear date. 20 to put into voltage amplitude, not power


    valid = np.where(np.logical_and(fq < .200, fq > .100)) # restrict to PAPER band
    
    fq, d, db, ph = fq[valid],d[valid], db[valid], ph[valid]

    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    window = a.dsp.gen_window(fq.size, WINDOW)
   

    if True: # approx account for 1st reflection of sky signal off of feed
        #_dw *= np.abs(_d[0])
        #_d *= np.abs(_d[0])
        
        #d -= np.abs(_d[0])
        #d1 = (np.abs(d)-np.abs(d_f))*np.abs(d_f)/(1-np.abs(d_f))+ (1-np.abs(d_f))
        #d_t = 1+d_f
        #d1 = (d-d_f)
        #d2 = d1*d_f
        #d3 = d2/d_t
        
        _dw1 = np.fft.ifft((d)**2*window) /window.mean()# compensate for window amplitude
        #_dw2 = np.fft.ifft((d_f)**2*window) /window.mean()# compensate for window amplitude
        #_dw3 = np.fft.ifft((d1)**2*window) /window.mean()# compensate for window amplitude
        #_dw4 = np.fft.ifft((d2)**2*window) /window.mean()# compensate for window amplitude
        #_dw5 = np.fft.ifft((d3)**2*window) /window.mean()# compensate for window amplitude
        #_dw6 = np.fft.ifft((d_t)**2*window) /window.mean()# compensate for window amplitude
       
    
    #import IPython; IPython.embed()
    print np.abs(d[0])
    
    #plt.plot(fq_f,((1+d_f)**2).real,linewidth=2.5,label='Gamma_t(real)')
    #plt.plot(fq_f,((1+d_f)**2).imag,linewidth=2.5,label='Gamma_t (imag)')
    #plt.plot(fq,10.0*np.log10(np.abs((d)**2)),linewidth=2.5,label='Feed Abs amp')
    #plt.plot(fq,np.abs(d_t),linewidth=2.5,label='Gamma_t Abs amp')
    #plt.plot(fq,np.abs(d_f),linewidth=2.5,label='Gamma_f Abs amp')
   
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(_dw1.real)), linewidth=2.5, label='Dish+Feed')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(_dw2.real)), linewidth=2.5, label='Feed ')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(_dw3.real)), linewidth=2.5, label='Dish+Feed-Feed')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(_dw4.real)), linewidth=2.5, label='(Dish+Feed-Feed)*Feed')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(_dw5.real)), linewidth=2.5, label='(Dish+Feed-Feed)*Feed/(1+Feed)')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(_dw6.real)), linewidth=2.5, label='[(Dish+Feed-Feed)*Feed/(1+Feed)]+(1+FEED)')
    
    plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dw1)**1)), linewidth=2.5, label=BASE)
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dw2)**1)), linewidth=2.5, label='Feed ')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dw3)**1)), linewidth=2.5, label='Dish+Feed-Feed')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dw4)**1)), linewidth=2.5, label='(Dish+Feed-Feed)*Feed')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dw5)**1)), linewidth=2.5, label='(Dish+Feed-Feed)*Feed/(1+Feed)')
    #plt.plot(np.fft.fftshift(tau), 10.0*np.log10(np.fft.fftshift(np.abs(_dw6)**1)), linewidth=2.5, label='[(Dish+Feed-Feed)*Feed/(1+Feed)]+(1+FEED)')
    
    
    
plt.xlim(-350,350)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Delay spectrum magnitude (dB)')
#plt.ylabel('Return loss magnitude (dB)')
plt.grid()
plt.legend(loc='lower right')
plt.show()


	







