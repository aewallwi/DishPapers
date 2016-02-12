#! /usr/bin/env python
import numpy as n
import pylab as p

#In plots dashed curves are delay power spectra with a blackman harris window across the subband. solid curves are power spectra attained via oqe formalism.

#names = ['pspecs_achrmbeam_sbinfo.npz', 'pspecs_chrmbeam_sbinfo.npz', 'pspecs_funcbeam_sbinfo.npz']
names = ['pspecs_chrmbeam_sbinfo.npz']
files = [ n.load(k) for k in names ]

#shape of pC and pI is (3,30,49) = (subband, bl, ks)

sband=0
lst=12

colors = ['k', 'c', 'm']
for b,bl in enumerate(files[0]['bls']):
    print b, bl
    for k,f in enumerate(files):
        p.subplot(111)
        p.semilogy(f['hetas'][k],n.abs(f['pCnorm'][sband,b,:,12]), colors[k], label=names[k].split('_')[1])
        p.semilogy(f['wetas'][k],n.abs(f['pInorm'][sband,b,:,12]),'--'+colors[k])
    
    p.suptitle('[ '+', '.join(map(str,bl)) + ' ]')
    p.legend()
    p.show()    
         
