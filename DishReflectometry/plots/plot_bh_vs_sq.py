#! /usr/bin/env python
import numpy as n 
import pylab as p
import sys, csv
import aipy as a

def fromcsv(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = n.array(list(d)[18:-1], dtype=n.float)
    return x[:,0]/1e9, x[:,1]

def take_delay(db, ph, fq, window='blackman-harris'):
    '''Take reflectometry data in dB and phase to return delay transform.
       Returns windowed and non windowed delay spectra for input spectrum.'''
    d = 10**(db/20) * n.exp(2j*n.pi*ph/360)
    tau = n.fft.fftfreq(fq.size, fq[1]-fq[0])
    window = a.dsp.gen_window(fq.size, window)
    _d = n.fft.ifft(d)
    _dw = n.fft.ifft(d*window) / window.mean() #compensate for window amplitude
    
    if True:
        _dw *= ( n.abs(_dw[0])/ (1- n.abs(_dw[0])))  # these should be changed to the dc bin of the windowed data.
        _d *= ( n.abs(_d[0])/ (1- n.abs(_d[0])))  # these should be changed to the dc bin of the windowed data.

    return n.fft.fftshift(_dw), n.fft.fftshift(_d), n.fft.fftshift(tau)

colors = n.array([(31,119,180), (255,127,14), (44,160,44), (214,39,40), (127,127,127), (148,103,189)])/255.


#file_base = sys.argv[1]
file_base = '../alldata/NC41_12'
amp = '_DB.csv'
phs = '_P.csv'
fq, amps = fromcsv(file_base + amp)
fq, phs= fromcsv(file_base + phs)

valid = n.where(n.logical_and(fq>.1, fq<.2))
dw, d, tau = take_delay(amps[valid], phs[valid], fq[valid])
dwhm, d, tau = take_delay(amps[valid], phs[valid], fq[valid], window='hamming')
dwhn, d, tau = take_delay(amps[valid], phs[valid], fq[valid], window='hanning')
p.plot(tau, 10*n.log10(n.abs(dw)**2), linewidth=2, label='blackman-harris', color=colors[0])
p.plot(tau, 10*n.log10(n.abs(d)**2), linewidth=2, label='square', color = colors[1])
p.plot(tau, 10*n.log10(n.abs(dwhm)**2), linewidth=2, label='hamming', color=colors[2])
p.plot(tau, 10*n.log10(n.abs(dwhn)**2), linewidth=2, label='hanning', color=colors[3])

p.xlim(-30,350) 
p.ylim(-100, 1)
p.vlines(60, -100,100, linestyle='--', linewidth=2)
p.hlines(-60,-100 ,500, linestyle='--', linewidth=2)
p.xlabel('delay (ns)')
p.ylabel('return loss (dB)')
p.grid(1)
p.legend() 

p.show()
