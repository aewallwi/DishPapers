#! /usr/bin/env python
import numpy as n, pylab as p, capo as C, aipy as a
import sys, scipy, csv

# 'bl' -- baseline vectors of size 30x3 for 30 baselines and three columns for x,y,z
# 'lst' -- LST during the drift scan. There are 40 LST samples at a cadence of 36 minutes thereby covering 24 hours
# 'fg_lags' -- delays appropriate for foreground data of size n(fg_lags)
# 'fg_power' -- Foreground power in power spectrum units and is an array of size n(bl) x n(fg_lags) x n(lst)
# 'eor_lags' -- delays appropriate for EoR data of size n(eor_lags)
# 'eor_power' -- EoR power in power spectrum units and is an array of size n(bl) x n(eor_lags) x n(lst)

npz = n.load(sys.argv[-1])
print npz.files

LST = 0

dk_deta = C.pspec.dk_deta(9.) 
bls = npz['bl']
tau1 = npz['fg_lags'] * 1e9 # ns
tau2 = npz['eor_lags'] * 1e9 # ns
fg = npz['fg_power'] * 1e6 # mK^2
eor = npz['eor_power'] * 1e6 # mK^2

#C.arp.waterfall(n.average(fg[:3], axis=0)); p.show()
fg = n.average(fg[...,:25], axis=-1)
eor = n.average(eor, axis=-1)

fg = n.average(fg[:3], axis=0)
eor = n.average(eor[:3], axis=0)
fg_mdl = n.where(fg > 1e7, fg, 0)

for kpl_cut in [0.08, 0.1, 0.12, 0.15][::-1]: # h Mpc^-1
    tau_cut = kpl_cut / dk_deta
    tau2_cut = n.argmin(n.abs(tau2 - tau_cut))
    SZ = eor.shape[0]
    eor_fold = 0.5 * (eor[SZ/2:] + eor[-SZ/2::-1])
    tau_fold = tau2[SZ/2:]
    eori = scipy.interpolate.interp1d(tau_fold, eor_fold, bounds_error=False)
    print tau_cut
    tau = n.linspace(0,500,1e3)
    resp_bound = n.ones_like(tau)
    for j,t in enumerate(tau1):
        if fg_mdl[j] <= 0: continue
        resp = n.sqrt(eori(tau+t) / fg_mdl[j])
        resp = n.where(tau+t < tau_cut, 1, resp)
        resp_bound = n.where(resp < resp_bound, resp, resp_bound)
    #p.plot(tau,10*n.log10(resp_bound))
    resp_bound = 10*n.log10(resp_bound)
    zeros = n.zeros_like(resp_bound)
    p.plot(tau,resp_bound,'--',label='$k_\\parallel=%4.2f$ spec'%kpl_cut)
    p.fill_between(tau, zeros, resp_bound, where=zeros>=resp_bound, facecolor='black', alpha=.5)
#p.fill_between(tau,0,10*n.log10(resp_bound), 'k', alpha=.5)

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
file_base = '/Users/aparsons/Documents/2015-08-27_GB_reflectometry/alldata/NC41_12'
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

p.xlim(-30,500) 
p.ylim(-100, 0)
#p.vlines(60, -100,100, linestyle='--', linewidth=2)
#p.hlines(-60,-100 ,500, linestyle='--', linewidth=2)
p.xlabel('delay (ns)')
p.ylabel('return loss (dB)')
p.grid()
p.legend() 

p.show()
