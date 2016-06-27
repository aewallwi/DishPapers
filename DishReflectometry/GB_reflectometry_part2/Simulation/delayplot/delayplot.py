import numpy as np
import matplotlib.pyplot as plt
import aipy as ap

######Read in data#####
ant_file = 'useHERAs11.csv'
feed_file= 'useHERA_feeds11.csv'
def readfile(fn):
    fp = open(fn,'r')
    i=0
    freq=[]
    s11_mag_dB=[]
    s11_phase_deg=[]
    s11 = []
    for line in fp:
        #print line,
        if i>0 and len(line)>2:
            data = line.split(',')
            freq.append(float(data[0]))
            m = np.power(10.0,float(data[1])/20.0)
            s11_mag_dB.append(float(data[1]))
            p = np.pi*float(data[2])/180.0
            s11_phase_deg.append(float(data[2]))
            s11.append(m*np.exp(1j*p))
        i+=1
    fp.close()
    freq = np.array(freq)
    s11_mag_dB = np.array(s11_mag_dB)
    s11_phase_deg = np.array(s11_phase_deg)
    s11=np.array(s11)
    return freq,s11_mag_dB,s11_phase_deg, s11
    
fa,ma,pa,S11=readfile(ant_file)
ff,mf,pf,Gf=readfile(feed_file)
Ef = 1.0 + Gf

plt.figure('S11')
plt.plot(fa,ma,label='antenna')
plt.plot(ff,mf,label='feed')
plt.legend()

def csq(c):
    cc = c*c.conjugate()
    return np.real(cc.real)

#####Compute Sigma#####
plt.figure('Sigma_f')
#####       S     #####
S = Gf*(S11-Gf)/Ef + Ef
plt.plot(fa,csq(S))

#####Analytical |S|^2#####
t1 = csq(Ef)
_t = Gf*(S11-Gf)/Ef
t2 = 2.0*np.real(_t.real)
t3 = csq(Gf)*csq(S11-Gf)/csq(Ef)

#S2 = t1 + t2 + t3
S2 = csq(S)
#S2 = S2*S2.conjugate()
plt.plot(fa,S2)

#####Get delay spectrum version#####
#####    Change type here     #####
S2use = 'square'  # 'linear' or 'square'
####################################
if S2use[0] == 'l':
    use_S = S
    Sscale = -20.0
    label = 'Eq 12'
    tloc = [40,25]
else:
    use_S = S2
    Sscale = -10.0
    label = 'Eq 13'
    tloc = [100,17]
WINDOW = 'blackman-harris'
plt.figure('delay_spec')
tau = np.fft.fftfreq(fa.size, (fa[1]-fa[0])*1e6)
tau = np.fft.fftshift(tau)
if WINDOW == 'uniform':
    window = np.ones(fa.size)
else:
    window = ap.dsp.gen_window(fa.size,WINDOW)
function = use_S*window
DS = np.fft.ifft(function)/window.mean()
DS = np.fft.fftshift(DS)
tau_hera = tau*1E9
dhera = Sscale*np.log10(np.abs(DS))
plt.plot(tau_hera,dhera)
#plt.plot(tau*1e9,Sscale*np.log10(np.abs(DS)))
#plt.text(tloc[0],tloc[1],label)
plt.show()
np.savez('delay_spectrum.npz',tau_hera=tau_hera,dhera=dhera)