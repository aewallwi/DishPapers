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

    valid = np.ones(fq.size, dtype=np.bool) # use entire sampled band
    #valid = np.where(fq < .250) # restrict to HERA band
    #valid = np.where(np.logical_and(fq < .2, fq > .1)) # restrict to PAPER band
    fq, d = fq[valid], d[valid]
    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    window = a.dsp.gen_window(fq.size, WINDOW)
    _d = np.fft.ifft(d)
    _dw = np.fft.ifft(d*window) / window.mean() # compensate for window amplitude

    if True: # approx account for 1st reflection of sky signal off of feed
        _dw *= np.abs(_d[0])
        _d *= np.abs(_d[0])

    print np.abs(_d[0])
    #plt.figure(1); plt.semilogy(fq, np.abs(d)**2, label=BASE)

    #plt.figure(2)
    #plt.semilogy(np.fft.fftshift(tau), np.fft.fftshift(np.abs(_d)**2))
    #plt.semilogy(np.fft.fftshift(tau), np.fft.fftshift(np.abs(_dw)), label=BASE)
    plt.semilogy(np.fft.fftshift(tau), np.fft.fftshift(np.abs(_dw)**2), label=BASE)

#plt.figure(1); plt.legend(loc='best')
#plt.figure(2)
plt.xlim(-30,350)
plt.ylim(1e-10,1)
plt.grid()
plt.legend(loc='best')
plt.show()
