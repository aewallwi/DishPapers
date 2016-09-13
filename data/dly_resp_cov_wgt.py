#! /usr/bin/env python
import numpy as np, pylab as plt
import capo, aipy
import capo.oqe as oqe
import sys

def gen_constraint(tau, achr, eor, tcut=200.):
    mn = {}
    for i,ti in enumerate(tau):
        for j,tj in enumerate(tau):
            if abs(tj) < tcut: continue
            dt = np.abs(np.around(5*(tj-ti), 0)/5)
            #mn[dt] = min(abs(eor[j] / achr[i]), mn.get(dt,1)) # if in amplitude
            mn[dt] = min(eor[j] - achr[i], mn.get(dt,0)) # if in dB
    dts = mn.keys(); dts.sort()
    mns = np.array([mn[dt] for dt in dts])
    return np.array(dts), mns

#def dB(sig): return 10*np.log10(np.average(sig.real, axis=1)) / 2. # from mK^2 to mK
def dB(sig): return 10*np.log10(np.abs(np.average(sig.real, axis=1))) / 2. # from mK^2 to mK
#def dB(sig): return 10*np.log10(np.abs(np.average(np.abs(sig), axis=1))) / 2. # from mK^2 to mK
#def dB(sig): return 10*np.log10(np.abs(np.median(sig, axis=1))) / 2. # from mK^2 to mK

#def wideband_clean(fg):
#    window = aipy.dsp.gen_window(fg.shape[-1], 'blackman-harris')
#    window.shape = (1,-1)
#    _fg = np.fft.ifft(window*fg)
#    _fg[:,:8] = 0; _fg[:,-7:] = 0
#    #_fg[:,:5] = 0; _fg[:,-4:] = 0
#    return np.fft.fft(_fg)/window

npz = np.load(sys.argv[-1])
dk_deta = capo.pspec.dk_deta(9.)
bls = npz['bl']
bl0 = 0 # do shortest baselines
#bl0 = 3 # do 2nd shortest baselines
print bls[bl0:bl0+3]
fqs = npz['freq'] / 1e9 # GHz
sky = npz['skyvis_freq'].astype(np.complex128) * 1e3 # mK (bl,fq,t)
sig = npz['vis_noise_freq'].astype(np.complex128) * 1e3 # mK (bl,fq,t)

npz = np.load('delayspectrum.npz')
tau_meas = npz['delay']
ref_meas = npz['amp']; ref_meas -= ref_meas[np.where(tau_meas == 0)]


CH0,NCHAN = 30, 61
#CH0,NCHAN = 0, 127
#fqs = np.linspace(.1,.2,NCHAN)
#fqs = np.linspace(.130,.160,NCHAN)
fqs = fqs[CH0:CH0+NCHAN]
dly = np.fft.fftfreq(NCHAN, fqs[1]-fqs[0])
k0 = ('even',(0,1),'I')
k1 = ('even',(1,2),'I')
k2 = ('even',(2,3),'I')
k3 = ('even',(3,4),'I')
fg,eor = {},{}
if True:
    fg[k0]  = sky[bl0+0,CH0:CH0+NCHAN].T
    fg[k1]  = sky[bl0+1,CH0:CH0+NCHAN].T
    fg[k2]  = sky[bl0+2,CH0:CH0+NCHAN].T
    eor[k0] = sig[bl0+0,CH0:CH0+NCHAN].T
    eor[k1] = sig[bl0+1,CH0:CH0+NCHAN].T
    eor[k2] = sig[bl0+2,CH0:CH0+NCHAN].T
else:
    NSAMP = 100
    ts = np.linspace(0,2*np.pi,NSAMP)
    NFG = 400
    scalar = 1

    def mk_dly_profile(dbs, slopes):
        prf = 0.
        for db,slope in zip(dbs,slopes):
            prf += 10**((db + slope * np.abs(dly))/10.)
        phs = np.random.uniform(0,2*np.pi,size=NCHAN)
        phs[0] = 0; phs[1:NCHAN/2+1] = phs[NCHAN/2+1:]
        return np.fft.fft(prf*np.exp(1j*phs))
    #bp = mk_dly_profile([0,-10],[-1.,-.1]))
    fg = 0
    for cnt in xrange(NFG):
        fg_ch = 1e3 * (fqs/.150)**np.random.uniform(-2,0); fg_ch.shape = (-1,1)
        bm = mk_dly_profile([np.random.uniform(-10,0),np.random.uniform(-40,0)],[np.random.uniform(-1,-2),np.random.uniform(-.1,-.2)]); bm.shape = (-1,1)
        fg_t = np.sin((cnt+1)*ts); fg_t.shape = (1,-1)
        fg += bm * fg_ch * fg_t
    eor = .01*oqe.noise(size=(NCHAN,NSAMP))

f = 0.2 # inject less eor so pspec goes below the eor signal used for computing ratios
dat = {}
for k in fg: dat[k] = (fg[k] + f*eor[k])
dat_cut = {}
for k in dat: dat_cut[k] = np.concatenate([dat[k][:54],dat[k][65:]], axis=0)
ds = oqe.DataSet(dsets=dat)

# gather C,iC matrices for baselines to try cross-application
Cs,iCs = {},{}
for k in dat:
    #Cs[k] = ds.C(k)
    Cs[k] = sum([ds.C(ki) for ki in dat if ki != k])
    #ds.set_C({k:Cs[k]+1e-1*np.identity(NCHAN)}) # regularize a bit with some diagonal
    ds.set_C({k:Cs[k]})
    iCs[k] = ds.iC(k)

ds.set_data(dat_cut)
tau = np.fft.fftshift(dly)
window = aipy.dsp.gen_window(NCHAN, 'blackman-harris'); window.shape = (-1,1)
for ki in iCs:
    qI = ds.q_hat(ki,ki,use_cov=False)
    FI = ds.get_F(ki,ki,use_cov=False)
    MI,WI = ds.get_MW(FI, mode='I')
    pI = ds.p_hat(MI,qI)
    ds.set_data(eor)
    qI_eor = ds.q_hat(ki,ki,use_cov=False)
    pI_eor = ds.p_hat(MI,qI_eor)
    ds.set_data(dat_cut)
    #pW = 1.6*2*np.abs(np.fft.fftshift(np.fft.ifft(window*fg[ki].T, axis=0), axes=0))**2
    pW = 1.6*2*np.abs(np.fft.fftshift(np.fft.ifft(window*dat_cut[ki].T, axis=0), axes=0))**2
    plt.figure(1)
    plt.plot(tau, dB(pI), 'b', label='I')
    plt.plot(tau, dB(pW), 'g', label='W')
    plt.plot(tau, dB(pI_eor), 'k', label='E')
    ds.set_iC({ki:iCs[ki]})
    qC = ds.q_hat(ki,ki)
    FC = ds.get_F(ki,ki)
    #MC,WC = ds.get_MW(FC, mode='F^-1/2')
    MC,WC = ds.get_MW(FC, mode='F^-1')
    pC = ds.p_hat(MC,qC)
    plt.figure(1)
    plt.plot(tau, dB(pC), 'r', label='C')
    plt.figure(2)
    for kcut in [.1,.15,.2]:
        tcut = kcut / dk_deta
        tau_c,prf_c = gen_constraint(tau, dB(pW), dB(pI_eor), tcut=tcut)
        plt.plot(tau_c, prf_c, 'g', label='W-%0.2f' % (kcut))
        tau_c, prf_c = gen_constraint(tau, dB(pC), dB(pI_eor), tcut=tcut)
        plt.plot(tau_c, prf_c, 'r', label='C-%0.2f' % (kcut))
#plt.figure(1); plt.legend(loc='best')
plt.figure(2)
plt.plot(tau_meas, ref_meas, 'k')
#plt.legend(loc='best')
plt.xlim(0,500)
plt.show()

