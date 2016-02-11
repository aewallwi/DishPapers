#! /usr/bin/env python
import numpy as n
import pylab as p

#In plots dashed curves are delay power spectra with a blackman harris window across the subband. solid curves are power spectra attained via oqe formalism.

names = ['pspecs_achrmbeam_sbinfo.npz', 'pspecs_chrmbeam_sbinfo.npz', 'pspecs_funcbeam_sbinfo.npz']
files = [ n.load(k) for k in names ]

#shape of pC and pI is (3,30,49) = (subband, bl, ks)

sband=0

colors = ['k', 'c', 'm']
for b,bl in enumerate(files[0]['bls']):
    for k,f in enumerate(files):
        p.subplot(111)
        p.semilogy(n.abs(f['pC'][sband,b,...]), colors[k], label=names[k].split('_')[1])
        p.semilogy(n.fft.fftshift(n.abs(f['pI'][sband,b,...])),'--'+colors[k])
    
    p.suptitle('[ '+', '.join(map(str,bl)) + ' ]')
    p.legend()
    p.show()    
         
