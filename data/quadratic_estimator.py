#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a
import sys, optparse, capo
import ipdb as db

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

#window function to apply to data (in quad est) if necessary
WINDOW='none'
PLOT=True
HALF=True

def Formatter(object):
    def __init__(self, im):
        self.im = im
    def __call__(self, x, y):
        z = self.im.get_array()[int(y), int(x)]
        return 'x={:1e-10f},y={:1e-10f},z={:1e-10f}'.format(x,y,z)

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex128)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

def get_Q(mode, n_k):
    '''Gets the Q matrix in the quadratic estimator formalism. 
       See Ali et. al. eq.12'''
    _m = n.zeros((n_k,), dtype=n.complex128)
    _m[mode] = 1.
    m = n.fft.fft(n.fft.ifftshift(_m)) * a.dsp.gen_window(nchan, WINDOW)
    Q = n.einsum('i,j', m.conj(), m)
    print n.diag(Q)
#    capo.arp.waterfall(Q, mode='phs')
#    p.colorbar(shrink=.5)
#    p.show()
    return Q

def read_npz(f, half=HALF,ntimes=None):
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
        data[i]=  vis[i,:,:ntimes]
        for k,midch in enumerate(subcent):
            mask = n.where(fwgts[k]>0.0)
            #mask = [ mask[0]-1, mask, mask[-1]+1 ]
            wfreqs[midch] = mask
            datawind[i][midch] = bandvis[i,k,mask,:ntimes]
            if half:
                #half the bandwidth for subdata
                nchan = len(mask[0])
                mask=mask[0][nchan-nchan/2-1:nchan+nchan/2]
                hfreqs[midch] = mask
                subdata[i][midch] = vis[i,mask,:ntimes]
            else:
                subdata[i][midch] = vis[i,mask,:ntimes]

    return data, datawind, subdata, lst, freq, bls, fwgts, wfreqs, hfreqs

def etas(freqs):
    return n.fft.fftshift(capo.pspec.f2eta(freqs))


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
nsubs = len(tmp['subbands'])
#ntimes = len(tmp['lst'])
ntimes = 50
#get Q matrix
Q = [get_Q(i, nchan) for i in xrange(nchan)]

#q^hat = x^t*C^-1*Q_alpha*C^-1*x
#p = Mq where W is a normalization matrix related to the M matrix

if HALF:
    banddataC = n.zeros(shape=(nsubs, nbls, nchan, ntimes), dtype=n.complex128) #unnorm data
    banddataC_norm = n.zeros(shape=(nsubs, nbls, nchan, ntimes), dtype=n.complex128) #norm data
#else:
#    banddataC = n.zeros(shape=(nsubs, nbls, nchan, ntimes)) 
#    banddataC_norm = n.zeros(shape=(nsubs, nbls, nchan, ntimes)) 
banddataI = n.zeros(shape=(nsubs, nbls, nchan*2 -1, ntimes), dtype=n.complex128) 
banddataI_norm = n.zeros(shape=(nsubs, nbls, nchan*2 -1, ntimes), dtype=n.complex128) 

#get power spectra for each baseline and frequency. 
for f in files:
    x, xw, xs, lst, freq, bls, fwgts, wfreqs, hfreqs = read_npz(f, ntimes=ntimes)
#    import IPython; IPython.embed()
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
            #C[bl][band] = n.identity(len(xs[bl][band]), dtype=n.complex64)
            U,S,V = n.linalg.svd(C[bl][band].conj())
            _C[bl][band] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
            _Cx[bl][band] = n.dot(_C[bl][band], xs[bl][band].squeeze())
    
            if PLOT and 0:
                #p.subplot(311); capo.arp.waterfall(xs[bl][band].squeeze(), mode='phs', extent=(lst[0],lst[-1],0,49))
                p.subplot(311); capo.arp.waterfall(x[bl].squeeze(), mode='log')
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

    FC = n.zeros((nchan,nchan), dtype=n.complex128)
    qC = n.zeros((nchan, _Cx.values()[0].values()[0].shape[1]), dtype=n.complex128)
    Q_Cx = {}
    hetas = n.zeros(shape=(3,25))
    wetas = n.zeros(shape=(3,49))
    for bl in xs.keys():
        if not Q_Cx.has_key(bl): Q_Cx[bl] = {}
        for bb,band in  enumerate(xs[bl].keys()):
            if not Q_Cx[bl].has_key(band): Q_Cx[bl][band] = [n.dot(Q[i], _Cx[bl][band]) for i in xrange(nchan)]
                        
            _qC = n.array([_Cx[bl][band].conj() * Q_Cx[bl][band][i] for i in xrange(nchan)])
            qC =+ n.sum(_qC, axis=1)
            #db.set_trace()
            if n.any(qC.imag/n.max(n.abs(qC)) > 1e-14): 
                raise ValueError('Significant imaginary component found!!!, bl=%s, band=%d'%(bl,bb))
            else:
                qC.imag = 0.0

            if PLOT and 0:
                print qC.imag / qC.real
                p.subplot(111); capo.arp.waterfall(qC.imag/qC.real, mode='lin'); p.colorbar(shrink=.5)
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
            #U,S,V = n.linalg.svd(FC.conj())
            #MC = n.dot(n.transpose(V), n.dot(n.diag(1./S), n.transpose(U)))
            
            MC = n.linalg.inv(FC) # normalization matrix

            if PLOT and 0:
                p.subplot(311); p.plot(qC.imag, 'o')
                p.subplot(312); p.plot(qC.real, 'o')
                p.subplot(313); p.plot(qC.imag/qC.real, 'o')
                #p.subplot(232); capo.arp.waterfall(MC, mode='imag'); p.colorbar(shrink=.5)
                #p.subplot(233); capo.arp.waterfall(FC, mode='real'); p.colorbar(shrink=.5)
                #p.subplot(234); capo.arp.waterfall(FC, mode='imag'); p.colorbar(shrink=.5)
                #p.subplot(235); capo.arp.waterfall(MC); p.colorbar(shrink=.5)
                #p.subplot(236); capo.arp.waterfall(n.dot(MC,FC), mode='lin'); p.colorbar(shrink=.5)
                p.show()
            
            WC = n.dot( MC, FC)
            norm = WC.sum(axis=-1); norm.shape += (1,)
            MC /= norm; WC = n.dot(MC, FC)

            print n.sum(WC, axis=-1)



            pC = n.dot(MC, qC)
            xs[bl][band]= xs[bl][band].squeeze()
            win = a.dsp.gen_window(nchan, window='blackman-harris')
            win.shape = (-1,1)
            #dp = n.fft.ifft(win*xs[bl][band],axis=0)*n.conj(n.fft.ifft(win*xs[bl][band],axis=0)) 
            sdf = freq[2] - freq[1]
            dp = n.fft.ifft(sdf*xw[bl][band].shape[1]*xw[bl][band],axis=1)*n.conj(n.fft.ifft(sdf*xw[bl][band].shape[1]*xw[bl][band],axis=1)) 
#            if bl==2:
            #import IPython; IPython.embed()
            dp = dp.squeeze()
            banddataC[bb][bl] = pC
            banddataI[bb][bl] = n.fft.fftshift(dp, axes=0)

            scalar = scalardict[f[:4]][bb]
            banddataC_norm[bb][bl] = banddataC[bb][bl]*scalar*(nchan*(freq[2]-freq[1]))**2 #need extra factor of bw
            banddataI_norm[bb][bl] = banddataI[bb][bl]*scalar#*(nchan*(freq[2]-freq[1]))**2
            
            
            #dp = dp.squeeze()
            #import IPython; IPython.embed()
            if PLOT and 0:
                #p.subplot(321); capo.arp.waterfall(qC.imag/qC.real, mode='lin') ; p.colorbar(shrink=.5)
                p.subplot(321); capo.arp.waterfall(FC.real) ; p.colorbar(shrink=.5)
                p.subplot(322); capo.arp.waterfall(FC.imag) ; p.colorbar(shrink=.5)
                p.subplot(323); capo.arp.waterfall(MC.real) ; p.colorbar(shrink=.5)
                p.subplot(324); capo.arp.waterfall(MC.imag) ; p.colorbar(shrink=.5)
                p.subplot(313); capo.arp.waterfall(pC.real) ; p.colorbar(shrink=.5)
                p.suptitle('[ '+', '.join(map(str,bls[bl])) + ' ] at ' + str(band))
                p.show()

            hetas[bb] = etas(freq[hfreqs[band]])
            wetas[bb] = etas(freq[wfreqs[band]])
            if PLOT and 1 :
                p.subplot(111); p.semilogy(etas(freq[hfreqs[band]]),n.abs(banddataC[bb][bl][:,12]),color='c')
                p.subplot(111); p.semilogy(etas(freq[hfreqs[band]]),n.abs(banddataC_norm[bb][bl][:,12]), color='k')
                p.subplot(111); p.semilogy(etas(freq[wfreqs[band]]),n.abs(banddataI[bb][bl][:,12]), color='g')
                p.subplot(111); p.semilogy(etas(freq[wfreqs[band]]),n.abs(banddataI_norm[bb][bl][:,12]), color='m')
                p.suptitle('[ '+', '.join(map(str,bls[bl])) + ' ] at ' + str(band))
                p.show()

            

    n.savez('pspecs_'+f.split('/')[-1], bls=bls, pC=banddataC, pI=banddataI, pCnorm=banddataC_norm, pInorm=banddataI_norm, hetas=hetas, wetas=wetas)
