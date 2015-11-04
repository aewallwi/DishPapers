#! /usr/bin/env python
import numpy as n 
import pylab as p
import sys, csv
import aipy as a
import pd

def fromcsv(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = n.array(list(d)[18:-1], dtype=n.float)
    return x[:,0]/1e9, x[:,1]

def fromcsvdaisy_time(filename):
    ''' Returns delay(ps), power(dB) '''
    print 'Reading Daisy file ', filename
    d = csv.reader(open(filename, 'r'), delimiter=',')
    x = n.array(list(d)[12:-1], dtype=n.float)
    print x
    return x[:,0]*1e-3, 20.0*n.log10(x[:,1])

def fromcsvdaisy_freq(filename):
    ''' Returns delay(ps), power(dB) '''
    print 'Reading Daisy file ', filename
    d = csv.reader(open(filename, 'r'), delimiter=',')
    x = n.array(list(d)[12:-1], dtype=n.float)
    print x
    return x[:,0]/1e9, 20.0*n.log10(x[:,1]), x[:,2]

def take_delay(db, ph, fq, window='blackman-harris'):
    '''Take reflectometry data in dB and phase to return delay transform.
       Returns windowed and non windowed delay spectra for input spectrum.'''
    d = 10**(db/20) * n.exp(2j*n.pi*ph/360)
    tau = n.fft.fftfreq(fq.size, fq[1]-fq[0])
    window = a.dsp.gen_window(fq.size, window)
    _d = n.fft.ifft(d)
    _dw = n.fft.ifft(d*window) / window.mean() #compensate for window amplitude
    
    if True:
    #if False:
        
        _dw *= ( n.abs(_dw[0])/ (1- n.abs(_dw[0])))  # these should be changed to the dc bin of the windowed data.
        _d *= ( n.abs(_d[0])/ (1- n.abs(_d[0])))  # these should be changed to the dc bin of the windowed data.

    return n.fft.fftshift(_dw), n.fft.fftshift(_d), n.fft.fftshift(tau)


colors = n.array([(31,119,180), (255,127,14), (44,160,44), (214,39,40), (148,103,189)])/255.

#file_base = sys.argv[1]
file_base = '../alldata/NC41_12'
amp = '_DB.csv'
phs = '_P.csv'
dfile_time = 'Time/set1/TXT102.csv'
dfile_freq= 'Freq/set1/TXT101.csv'
fq, amps = fromcsv(file_base + amp)
fq, phs= fromcsv(file_base + phs)
dns, ddb = fromcsvdaisy_time(dfile_time)
dfreq, dfreqdb, dphs = fromcsvdaisy_freq(dfile_freq)

valids = {
          '50 - 250 MHz'  : n.where(n.logical_and(fq>.05 ,fq<.25)), 
          '100 - 200 MHz' : n.where(n.logical_and(fq>.1 ,fq<.2)), 
          '140 - 160 MHz' : n.where(n.logical_and(fq>.140 ,fq<.160)),
          'no cage: 50 - 1000 MHz' : None,
          'no cage: 100 - 200 MHz' : n.where(n.logical_and(dfreq>.1, dfreq<.2))
          #'100 - 200 MHz (no cage)' : n.ones(len(dfreq), dtype=n.bool)#n.where(n.logical_and(dfreq>.1, dfreq<.2))
#          'second' : n.where(n.logical_and(fq>.250 ,fq<.500))
         }

#fig, axes = p.subplots(3,1,figsize=(5.4,9))
##fig2, axes2 = p.subplots(3,1,figsize=(10,8))
#for i,v in enumerate(valids.keys()):
#    dw, d, tau = take_delay(amps[valids[v]], phs[valids[v]], fq[valids[v]])
#    print v
#    axes[i].plot(tau, 10*n.log10(n.abs(dw)**2), linewidth=3, label='%s'%v)
#    axes[i].set_xlim(-30,350) 
#    axes[i].set_ylim(-100, 1)
#    axes[i].vlines(60, -100,2, linestyle='--', linewidth=3)
#    axes[i].hlines(-60,-100 ,500, linestyle='--', linewidth=3)
#    axes[i].set_xlabel('delay (ns)')
#    axes[i].set_ylabel('return loss (dB)')
#    axes[i].grid(1)
#    axes[i].legend() 
plots = []
names = []
for i,v in enumerate(valids.keys()):
    
    if v == 'no cage: 50 - 1000 MHz':
        plots.append(p.plot(dns, ddb, linewidth=2, label='%s'%v, color=colors[i]))
        names.append(v)
        #print dns,ddb
    elif v == 'no cage: 100 - 200 MHz':
        dw, d, tau = take_delay(dfreqdb[valids[v]], dphs[valids[v]], dfreq[valids[v]], window='hamming')
        plots.append(p.plot(tau, 10*n.log10(n.abs(dw)**2), linewidth=2, label='%s'%v, color=colors[i]))
        names.append(v)
    else:
        dw, d, tau = take_delay(amps[valids[v]], phs[valids[v]], fq[valids[v]], window='hamming')
        plots.append(p.plot(tau, 10*n.log10(n.abs(dw)**2), linewidth=2, label='%s'%v, color=colors[i]))
        names.append(v)
plots = n.array(plots)
names = n.array(names)
p.xlim(-30,350) 
p.ylim(-100, 1)
p.vlines(60, -100,100, linestyle='--', linewidth=2)
p.hlines(-60,-100 ,500, linestyle='--', linewidth=2)
p.xlabel('delay (ns)')
p.ylabel('return loss (dB)')
p.grid(1)
plots = plots[[2,1,3,0,4]]
print names
names = names[[2,1,3,0,4]]
print names
names = [i for i in names]
plots = [i[0] for i in plots]
print names, plots
#p.legend()
p.legend(plots, names) 

p.show()
