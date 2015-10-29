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
#file_base = '../alldata/NC41_'
file_base = '../alldata/F'
amp_end = '_DB.csv'
phs_end = '_P.csv'

#valids = [6,9,12,15]
#valids = ['A_NC41_3','A_NC41_AB_3','C_NC41','C_NC41_AB','C_NC41_UD']
#labels = ['feed alone 3ft above ground','feed alone 3ft above absorber','feed in cage on ground','feed in cage on absorber','feed in cage upsidedown']
valids = ['C_NC41','C_NC41_AB']
labels = ['feed in cage on ground','feed in cage on absorber']

fig, ax = p.subplots(nrows=1,ncols=1)
for i,v in enumerate(valids):
    fq,amps = fromcsv(file_base + str(v) + amp_end)
    fq,phs = fromcsv(file_base + str(v) + phs_end)
    #bandwidth = n.where(n.logical_and(fq>.05 ,fq<.25)) #HERA bandwidth
    #bandwidth = n.where(n.logical_and(fq>.1,fq<.2)) #PAPER bandwidth
    bandwidth = n.where(n.logical_and(fq>.05,fq<.5)) #full bandwidth
    dw, d, tau = take_delay(amps[bandwidth], phs[bandwidth], fq[bandwidth])
#    p.plot(tau, 10*n.log10(n.abs(dw)**2), linewidth=2, label=('%s'%v)+'ft', color=colors[i])
    ax.plot(tau, 10*n.log10(n.abs(dw)**2), linewidth=2, label=labels[i], color=colors[i])
p.xlim(-30,350)
p.ylim(-100, 1)
p.vlines(60, -100,100, linestyle='--', linewidth=2)
p.hlines(-60,-100 ,500, linestyle='--', linewidth=2)
p.xlabel('delay (ns)')
p.ylabel('return loss (dB)')
p.grid(1)
p.legend()

fig.subplots_adjust(left=.08, top=.95, bottom=.10,right=0.92)

p.show()
