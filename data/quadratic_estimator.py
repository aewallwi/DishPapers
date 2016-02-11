#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a
import sys, optparse, capo

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

#window function to apply to data (in quad est) if necessary
WINDOW='none'
PLOT=True
HALF=True


def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

def get_Q(mode, n_k):
    '''Gets the Q matrix in the quadratic estimator formalism. 
       See Ali et. al. eq.12'''
    _m = n.zeros((n_k,), dtype=n.complex)
    _m[mode] = 1.
    m = n.fft.fft(n.fft.ifftshift(_m)) * a.dsp.gen_window(nchan, WINDOW)
    Q = n.einsum('i,j', m, m.conj())
    return Q

def read_npz(f, half=HALF):
    '''Return dictionary with keys [bl][cent] = vis'''
    data = {} # full visibility data (all chans)
    subdata = {} #a subset of the channels corresponding to the windowed channels but no window applied.
    datawind = {} #windowed data with blackman harris. Subset of chans
    wfreqs = {}
    hfreqs = {}
    npz = n.load(f)
    bls = npz['bl'] #30X3, 20 bls by x,y,z
    freq = npz['freq'] #256, in hz
    lst = npz['lst']*12/180. #len = 80, in hours
    subcent = npz['subbands'] #3 centers of the subbands.
    fwgts = npz['freq_wts'] #3X256 window functions for reach subband.
    vis = npz['skyvis_freq'] #30bls x 256 chs x 80 times
    bandvis = npz['subband_skyvis_freq'] #30bls x 3 bands x256 chs x 80 times
    for i,(x,y,z) in enumerate(bls):
        data[i] = {}
        subdata[i] = {}
        datawind[i] = {}
        data[i]=  vis[i,:,:]
        for k,midch in enumerate(subcent):
            mask = n.where(fwgts[k]>0.0)
            wfreqs[midch] = mask
            datawind[i][midch] = bandvis[i,k,mask,:]
            if half:
                #half the bandwidth for subdata
                nchan = len(mask[0])
                mask=mask[0][nchan-nchan/2-1:nchan+nchan/2]
                hfreqs[midch] = mask
                subdata[i][midch] = vis[i,mask,:]
            else:
                subdata[i][midch] = vis[i,mask,:]

    return data, datawind, subdata, lst, freq, bls, fwgts, wfreqs, hfreqs

def etas(freqs):
    return n.fft.fftshift(capo.pspec.f2eta(freqs))

#!!!!need to figure out normalization!!!!

files = args
#'bl' (size 30x3) baseline vectors (in m) in the simulation
#'freq' (size 256) frequencies (in Hz) covering the full bandwidth 
#'lst' (size 80) LST values (in degrees) in the simulation
#'subbands' (size 3) subband frequency centers (in Hz)
#'freq_wts' (size 3x256) frequency weights (Blackman-Harris) weights that will be applied on fullband visibilities to get subband visibilities
#'skyvis_freq' (size 30x256x80) Full band visibilities (in Jy) on all baselines, frequencies and LST
#'subband_skyvis_freq' (size 30x3x256x80) subband visibilities (in Jy) on all baselines, subbands, frequencies, and LST

#scalars to convert to K. (l^2/2kb)^2 * X^2Y/(omega*B) for each frequency window
scalardict = { 'achr': n.array([6.42416224e-09, 4.62634957e-09, 3.39232286e-09]),
               'chrm': n.array([6.40467344e-09, 4.57932240e-09, 3.87532971e-09]),
               'func': n.array([9.84182143e-09, 8.11682276e-09, 6.73087936e-09]) }

tmp = n.load(files[0])
nchan = len(tmp['freq'][tmp['freq_wts'][0]>0]) #take windowed frequencies.
if WINDOW=='none': nchan = nchan/2 + 1 #make it an odd number
nbls = len(tmp['bl']) #number of baselines
#get Q matrix
Q = [get_Q(i, nchan) for i in xrange(nchan)]

#q^hat = x^t*C^-1*Q_alpha*C^-1*x
#p = Mq where W is a normalization matrix related to the M matrix

if HALF:
    banddataC = n.zeros(shape=(3, nbls, 25)) #array that holds baseline pspec for three band
    banddataC_norm = n.zeros(shape=(3, nbls, 25)) #array that holds baseline pspec for three band
else:
    banddataC = n.zeros(shape=(3, nbls, 49)) #array that holds baseline pspec for three band
    banddataC_norm = n.zeros(shape=(3, nbls, 49)) #array that holds baseline pspec for three band
banddataI = n.zeros(shape=(3, nbls, 49)) #array that holds baseline pspec for three bands
banddataI_norm = n.zeros(shape=(3, nbls, 49)) #array that holds baseline pspec for three bands
#get power spectra for each baseline and frequency. 
for f in files:
    x, xw, xs, lst, freq, bls, fwgts, wfreqs, hfreqs = read_npz(f)
    #need to take half the channels for xs since not being windowed
    fin = {}# dictionary for final pspecs. bls only
    C,_C,_Cx,_CQ = {},{},{},{}
    for bl in xs.keys():
        fin[bl] = {}         
        C[bl] = {}
        _C[bl] = {}
        _Cx[bl] = {}
        _CQ[bl] = {}
        for band in xs[bl].keys():
            C[bl][band] = cov(xs[bl][band].squeeze())
            U,S,V = n.linalg.svd(C[bl][band].conj())
            _C[bl][band] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
            _Cx[bl][band] = n.dot(_C[bl][band], xs[bl][band].squeeze())
    
            if PLOT and 0:
                #p.subplot(311); capo.arp.waterfall(xs[bl][band].squeeze(), mode='phs', extent=(lst[0],lst[-1],0,49))
                p.subplot(311); capo.arp.waterfall(x[bl].squeeze(), mode='phs')
                p.colorbar(shrink=.5)
                p.subplot(323); capo.arp.waterfall(C[bl][band], mode='log')
                p.subplot(324); p.plot(n.einsum('ij,jk', n.diag(S), V).T.real)
                p.subplot(313); capo.arp.waterfall(_Cx[bl][band], mode='log')
                p.colorbar(shrink=.5)
                p.suptitle('[ '+', '.join(map(str,bls[bl])) + ' ] at ' + str(band))
                p.show()

            _CQ[bl][band]={}
            for ch in xrange(nchan):
                _CQ[bl][band][ch] = n.dot(_C[bl][band], Q[ch])

    FC = n.zeros((nchan,nchan), dtype=n.complex)
    qC = n.zeros((nchan, _Cx.values()[0].values()[0].shape[1]), dtype=n.complex)
    Q_Cx = {}
    for bl in xs.keys():
        if not Q_Cx.has_key(bl): Q_Cx[bl] = {}
        for bb,band in  enumerate(xs[bl].keys()):
            if not Q_Cx[bl].has_key(band): Q_Cx[bl][band] = [n.dot(Q[i], _Cx[bl][band]) for i in xrange(nchan)]
                        
            _qC = n.array([_Cx[bl][band].conj() * Q_Cx[bl][band][i] for i in xrange(nchan)])
            qC =+ n.sum(_qC, axis=1)
            if PLOT and 0:
                p.subplot(111); capo.arp.waterfall(qC); p.colorbar(shrink=.5)
                p.suptitle('[ '+', '.join(map(str,bls[bl])) + ' ] at ' + str(band))
                p.show()
            
            for i in xrange(nchan):
                for k in xrange(nchan):
                    FC[i,k]+= .5*n.einsum('ij,ji', _CQ[bl][band][i], _CQ[bl][band][k])
            if PLOT and 0: 
                p.subplot(111); capo.arp.waterfall(FC); p.colorbar(shrink=.5)
                p.suptitle('[ '+', '.join(map(str,bls[bl])) + ' ] at ' + str(band))

                p.show()



            # We are at the point where we have Fisher matrix and q-estimator.
            #NEed p = Mq. 
            #Take M = F^-1
            U,S,V = n.linalg.svd(FC.conj())
            #MC = n.dot(n.transpose(V), n.dot(n.diag(1./S), n.transpose(U)))
            MC = n.linalg.inv(FC)
            
            WC = n.dot( MC, FC)
            norm = WC.sum(axis=-1); norm.shape += (1,)
            MC /= norm; WC = n.dot(MC, FC)

            pC = n.dot(MC, qC)
            xs[bl][band]= xs[bl][band].squeeze()
            win = a.dsp.gen_window(nchan, window='blackman-harris')
            win.shape = (-1,1)
            #dp = n.fft.ifft(win*xs[bl][band],axis=0)*n.conj(n.fft.ifft(win*xs[bl][band],axis=0)) 
            dp = n.fft.ifft(xw[bl][band],axis=1)*n.conj(n.fft.ifft(xw[bl][band],axis=1)) 
            dp = dp.squeeze()
            banddataC[bb][bl] = n.average(pC, axis=-1)
            banddataI[bb][bl] = n.fft.fftshift(n.average(dp, axis=-1))

            scalar = scalardict[f[:4]][bb]
            banddataC_norm[bb][bl] = banddataC[bb][bl]*scalar*(nchan*(freq[2]-freq[1]))**2 #need extra factor of bw
            banddataI_norm[bb][bl] = banddataI[bb][bl]*scalar*(nchan*(freq[2]-freq[1]))**2
            
            #dp = dp.squeeze()
            #import IPython; IPython.embed()
            if PLOT and 0:
                p.subplot(311); capo.arp.waterfall(qC) ; p.colorbar(shrink=.5)
                p.subplot(323); capo.arp.waterfall(FC) ; p.colorbar(shrink=.5)
                p.subplot(324); capo.arp.waterfall(WC) ; p.colorbar(shrink=.5)
                p.subplot(313); capo.arp.waterfall(pC) ; p.colorbar(shrink=.5)
                p.suptitle('[ '+', '.join(map(str,bls[bl])) + ' ] at ' + str(band))
                p.show()

            if PLOT and 0 :
                p.subplot(111); p.semilogy(etas(freq[hfreqs[band]]),n.abs(banddataC[bb][bl]))
                p.subplot(111); p.semilogy(etas(freq[hfreqs[band]]),n.abs(banddataC_norm[bb][bl]))
                p.subplot(111); p.semilogy(etas(freq[wfreqs[band]]),n.abs(banddataI[bb][bl]), color='g')
                p.subplot(111); p.semilogy(etas(freq[wfreqs[band]]),n.abs(banddataI_norm[bb][bl]), color='m')
                p.suptitle('[ '+', '.join(map(str,bls[bl])) + ' ] at ' + str(band))
                p.show()

            

    n.savez('pspecs_'+f.split('/')[-1], bls=bls, pC=banddataC, pI=banddataI, pCnorm=banddataC_norm, pInorm=banddataI_norm)
