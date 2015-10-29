#! /usr/bin/env python
import numpy as n 
import pylab as p
import sys, csv

def fromcsv(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = n.array(list(d)[18:-1], dtype=n.float)
    return x[:,0]/1e9, x[:,1]

colors = n.array([(31,119,180), (255,127,14), (44,160,44), (214,39,40), (127,127,127), (148,103,189)])/255.

#file_base = sys.argv[1]
file_base = '../alldata/NC41_12'
amp = '_DB.csv'
phs = '_P.csv'
freqs, amps = fromcsv(file_base + amp)
freqs, phase= fromcsv(file_base + phs)

#rloss = amps * n.exp(2j*n.pi*phase/360.)
#p.figure(figsize=(10,5))
#p.plot(freqs, rloss.real, linewidth=3, label='real')
#p.plot(freqs, rloss.imag, linewidth=3, label='imag')
#p.xlabel('Frequency (GHz)' )
#p.ylabel('Return Loss (dB)')
#p.grid(1)
#p.xlim(.05,.250)
#p.ylim(-30,25)
#p.vlines(.1, -100, 100, linestyle='--', linewidth=2)
#p.vlines(.2, -100, 100, linestyle='--', linewidth=2)
#p.legend()


fig, (ax1,ax2) = p.subplots(2,1, figsize=(10,5), sharex=True)
#ax1.plot(freqs, amps, linewidth=3, label='amplitude', color=(0,107/255.,164/255.))
ax1.plot(freqs, amps, linewidth=3, label='amplitude', color=colors[0])
ax1.set_ylabel('Return Loss (dB)')
ax1.grid(1)
#ax1.set_xlim(.05,.250)
ax1.set_ylim(-25,25)
#p.setp(axes[0].get_xticklabels(), visible=False)
ax1.vlines(.1, -100, 100, linestyle='--', linewidth=2)
ax1.vlines(.2, -100, 100, linestyle='--', linewidth=2)
ax1.vlines(.05, -100, 100, linestyle='--', linewidth=2, color='c')
ax1.vlines(.25, -100, 100, linestyle='--', linewidth=2, color='c')
ax1.vlines(.14, -100, 100, linestyle='--', linewidth=2, color='m')
ax1.vlines(.16, -100, 100, linestyle='--', linewidth=2, color='m')
ax1.legend()

#ax2.plot(freqs, phase*n.pi/180., linewidth=3, label='phase', color=(1,128/255.,14/255.))
ax2.plot(freqs, n.unwrap(phase*n.pi/180), linewidth=3, label='phase', color=colors[1])
ax2.set_xlabel('Frequency (GHz)' )
ax2.set_ylabel('Phase (radians)')
ax2.grid(1)
#ax2.set_xlim(.05,.250)
ax2.set_ylim(-14, 0)
ax2.vlines(.1, -100, 100, linestyle='--', linewidth=2)
ax2.vlines(.2, -100, 100, linestyle='--', linewidth=2)
ax2.vlines(.05, -100, 100, linestyle='--', linewidth=2, color='c')
ax2.vlines(.25, -100, 100, linestyle='--', linewidth=2, color='c')
ax2.vlines(.14, -100, 100, linestyle='--', linewidth=2, color='m')
ax2.vlines(.16, -100, 100, linestyle='--', linewidth=2, color='m')
ax2.legend()

fig.subplots_adjust(left=.08, top=.95, bottom=.10,right=0.88)

p.show()
