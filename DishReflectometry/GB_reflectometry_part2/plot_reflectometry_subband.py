#! /usr/bin/env python

import numpy as np, pylab as plt, aipy as a
import sys, csv

def fromcsv(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = np.array(list(d)[18:-1], dtype=np.float)
    return x[:,0]/1e9, x[:,1]
    #d = 10**(db/10) * np.exp(2j*np.pi*ph/360) # power
    

def take_delay(db, db_f, ph, ph_f, fq, window='blackman-harris'):
    '''Take reflectometry data in dB and phase to return delay transform.
       Returns windowed and non windowed delay spectra for input spectrum.'''    
    d = 10**(db/20) * np.exp(2j*np.pi*ph/360) # 20 to put into voltage amplitude, not power
    d_f = 10**(db_f/20) * np.exp(2j*np.pi*ph_f/360) # 20 to put into voltage amplitude, not power
    d_t =1+d_f # Computing the transmission coefficient
    d1 = ((d-d_f)*d_f/d_t ) # Reflections corrected to produce system bandpass in the receiving mode 
    d2 = d1+d_t # Zero term added for final correction
    #d1 = (np.abs(d)-np.abs(d_f))
    #d1 = (np.abs(d)-np.abs(d_f))*np.abs(d_f)/(1-np.abs(d_f))+ (1-np.abs(d_f))
    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    window = a.dsp.gen_window(fq.size, window)
    _d = np.fft.ifft(np.abs(d2)**2.0)
    _dw = np.fft.ifft(np.abs(d2)**2.0*window) / window.mean() #compensate for window amplitude
    
    if True:
    #if False:
        
        #_dw *= ( np.abs(_dw[0])/ (1- np.abs(_dw[0])))  # these should be changed to the dc bin of the windowed data.
        #_d *= ( np.abs(_d[0])/ (1- np.abs(_d[0])))  # these should be changed to the dc bin of the windowed data.

     return np.fft.fftshift(_dw), np.fft.fftshift(_d), np.fft.fftshift(tau)

colors = np.array([(31,119,180), (255,127,14), (44,160,44), (214,39,40), (148,103,189)])/255.

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


valids = {
	          '100 - 130 MHz'  : np.where(np.logical_and(fq>.100 ,fq<.130)), 
    	      '130 - 160 MHz' : np.where(np.logical_and(fq>.130 ,fq<.160)), 
        	  '160 - 190 MHz' : np.where(np.logical_and(fq>.160 ,fq<.190)),
         }
plots = []
names = []
for i,v in enumerate(valids.keys()):
    
    if v == 'no cage: 50 - 1000 MHz':
        plots.append(plt.plot(dns, ddb, linewidth=2, label='%s'%v))
        names.append(v)
        #print dns,ddb
    elif v == 'no cage: 100 - 200 MHz':
        dw, d, tau = take_delay(dfreqdb[valids[v]], dphs[valids[v]], dfreq[valids[v]], window='hamming')
        plots.append(plt.plot(tau, 10*n.log10(n.abs(dw)), linewidth=2, label='%s'%v))
        names.append(v)
    else:
        dw, d, tau = take_delay(db[valids[v]], db_f[valids[v]], ph[valids[v]], ph_f[valids[v]], fq[valids[v]], window='blackman-harris')
        plots.append(plt.plot(tau, 10*np.log10(np.abs(dw)), linewidth=2, label='%s'%v))
        names.append(v)
plots = np.array(plots)
names = np.array(names)
plt.xlim(-350,350) 
#plt.ylim(-10, 1)
plt.xlabel('delay (ns)')
plt.ylabel('Delay spectrum (dB)')
plt.grid(1)
#plots = plots[[2,1,3,0,4]]
#print names
#names = names[[2,1,3,0,4]]
#print names
names = [i for i in names]
plots = [i[0] for i in plots]
print names, plots
#p.legend()
plt.legend(plots, names) 

plt.show()

          	 
          

#for i,v in enumerate(valids.keys()):

    #valid = np.ones(fq.size, dtype=np.bool) # use entire sampled band
    #valid = np.where(fq < .250) # restrict to HERA band
#    if v == '100-130 MHz': 
    #valid = np.where(np.logical_and(fq < .200, fq > .100)) # restrict to PAPER band
    
    #fq, d, db, ph, d_f, db_f, ph_f = fq[valid], d[valid], db[valid], ph[valid], d_f[valid], db_f[valid],ph_f[valid]
    #elif v == '130-160 MHz':
     # fq, d, db, ph = fq[valids[v]], d[valids[v]], db[valids[v]], ph[valids[v]]
      #print np.abs(d[0])
    #elif v == '160-190 MHz':    
     # fq, d, db, ph = fq[valids[v]], d[valids[v]], db[valids[v]], ph[valids[v]]
    






