#! /usr/bin/env python

import numpy as np, pylab as plt, aipy as ap
import sys, csv

def fromcsv(filename):
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = np.array(list(d)[18:-1], dtype=np.float)
    return x[:,0]/1e9, x[:,1]

def load_data(filename):
    assert(filename.endswith('_DB.csv'))
    BASE = filename[:-len('_DB.csv')]
    fq,db = fromcsv(BASE + '_DB.csv')
    fq,ph = fromcsv(BASE + '_P.csv')
    #d = 10**(db/10) * np.exp(2j*np.pi*ph/360) # power
    d = 10**(db/20) * np.exp(2j*np.pi*ph/360) # 20 to put into voltage amplitude, not power
    return fq, d

def load_rb_data(filename):
    d = np.loadtxt(filename, skiprows=9, delimiter=',')[:-1]
    d = d[:,0] + 1j*d[:,1]
    fq = np.linspace(.1,.2,1600)
    return fq, d

def dly_coeff_transform(fq, dish_volt, feed_volt):
    window = ap.dsp.gen_window(fq.size, 'blackman-harris')
    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    refl_coeff = np.abs(feed_volt)**2
    tran_coeff = 1 - refl_coeff
    coeff = refl_coeff / tran_coeff
    _dw = np.fft.ifft(window * np.abs(dish_volt)**2 * coeff) / window.mean() # compensate for window amplitude
    return np.fft.fftshift(tau), np.fft.fftshift(_dw)

def dly_transform(fq, d_volt):
    window = ap.dsp.gen_window(fq.size, 'blackman-harris')
    tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
    _d = np.fft.ifft(np.abs(d_volt)**2)
    _dw = np.fft.ifft(window * np.abs(d_volt)**2) / window.mean() # compensate for window amplitude
    if True: # approx account for 1st reflection of sky signal off of feed
        _dw *= np.abs(_d[0])
    return np.fft.fftshift(tau), np.fft.fftshift(_dw)

def valid_fqs(fq, d):
    #valid = np.ones(fq.size, dtype=np.bool) # use entire sampled band
    #valid = np.where(fq < .250) # restrict to HERA band
    valid = np.where(np.logical_and(fq < .2, fq > .1)) # restrict to PAPER band
    #valid = np.where(np.logical_and(fq < .19, fq > .11)) # restrict to PAPER band
    return fq[valid], d[valid]

fq,dish = load_rb_data('DishAndFeed/HERA_Dish1_S11_Feed5_3m_A.d1')
fq,feed = load_rb_data('FeedOnly/DATA02.d1')
fq,paper = load_rb_data('PAPER/DATA03.d1')
fq,balun = load_rb_data('PAPER/DATA02.d1')
tau,_dish = dly_coeff_transform(fq, dish, feed)
#tau,_paper = dly_transform(fq, paper)
tau,_paper = dly_coeff_transform(fq, paper, feed)
tau,_balun = dly_transform(fq, balun)
fig = plt.figure(figsize=(6,5))
fig.subplots_adjust(left=.15, top=.95, bottom=.15, right=.95)
plt.plot(tau, 10*np.log10(np.abs(_dish)), label='HERA reflect')
plt.plot(tau, 10*np.log10(np.abs(_paper)), label='PAPER reflect')
plt.plot(tau[800:805], 10*np.log10(np.abs(_balun[800:805])), label='Received')
plt.ylim(-65,0)
plt.xlim(0,500)
plt.ylabel('Power [dB]')
plt.xlabel('Delay [ns]')
plt.legend(loc='best')
plt.grid()
plt.show()
