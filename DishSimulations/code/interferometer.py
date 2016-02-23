#!/usr/bin/env python
#************************************************************
#interferometer class
#ignore polarization for now
#obs: an ephem observer with lon,lat,and date
#freqs: a list of frequencies for the interferometer (Hz)
#times: a list of ephem dates
#assume homogenous array for now
#************************************************************
import matplotlib.pyplot as plt
import astropy.wcs as wcs
import numpy as np
import ephem,os
from beam import Beam
from beam import gaussbeam,cosbeam,tophatbeam,airybeam,mwabeam
from bandpass import Bandpass
from bandpass import flatband
from arraytools import unique_rows
import ctypes as ct
from joblib import Parallel,delayed
import interferometer_utils as utils
import multiprocessing,cosmology
pi=np.pi
c=3e8
"""
#test diffcomb
bl=diffcomb(range(128))
print bl
#test compute visibilities
#with an interferometer with antennae spacing of .1 lambda and a source located at m=0.5
nuv=100
u=np.arange(-nuv/2,nuv/2)*.1
v=np.arange(-nuv/2,nuv/2)*.1
ug,vg=np.meshgrid(u,v)
u=ug.flatten()
v=vg.flatten()
w=np.zeros(len(u))
uvw=np.vstack([u,v,w]).T
ns=10
l=np.random.rand(ns)
lm=np.vstack([l,1-l**2.]).T
s=np.random.rand(ns)
print lm.shape
#lm=np.array([[0.,.5],[0.,0.]])
vis=computevisibilities(uvw,lm,s).reshape((nuv,nuv))
print vis.shape
theory=np.zeros((nuv,nuv))
for mm in range(lm.shape[0]):
    theory+=s[mm]/np.sqrt(1.-np.sum(lm[mm,:]**2.))*np.cos(-2*np.pi*(ug*lm[mm,0]+vg*lm[mm,1]))
plt.imshow((np.real(vis)-theory)/theory,interpolation='nearest')
plt.colorbar()
plt.show()
plt.imshow(np.real(vis),interpolation='nearest')
plt.colorbar()
plt.show()
"""



#************************************************************
#change coordinates of antennae in north-south,east-west
#at an observatory with latituded lat and local sidereal
#time lst to coordinates with X at 0h, Y at -6h and Z 
#through zenith (see Thompson Swenson Moran)
#************************************************************
def nsew2xyz(ewns,lst,lat):
    lat=np.radians(lat)
    lst=np.radians(lst)
    xyz=np.zeros(ewns.shape)
    xyz[:,2]=ewns[:,2]*np.sin(lat)+ewns[:,1]*np.cos(lat)
    xyz[:,1]=ewns[:,0]
    xyz[:,0]=ewns[:,2]*np.cos(lat)-ewns[:,1]*np.sin(lat)
    xyzt0=xyz[:,0]*np.cos(lst)+xyz[:,1]*np.sin(lst)
    xyzt1=xyz[:,1]*np.cos(lst)-xyz[:,0]*np.sin(lst)
    xyz[:,0]=xyzt0
    xyz[:,1]=xyzt1
    return xyz


def refModel(f,r,phi,tau):
    return 1./(1-np.exp(1j*(2*np.pi*f*tau+phi)))

    

class Interferometer():
    def __init__(self,lon,lat,antpos,freqs,df=40e3,date=ephem.now(),driftMode=True,beam=Beam(mwabeam,[None,[0.,0.]])):
        self.driftMode=driftMode
        self.obs=ephem.Observer()
        self.lat=lat
        self.lon=lon
        self.obs.lat=str(lat)
        self.obs.lon=str(lon)
        self.obs.date=date
        self.antpos=antpos
        #convert from n-s,e-w coords to XYZ coords where Y points to -6 hours, X to 0 hours and Z along the earth's rotation axis. 
        self.lst=np.degrees(float(repr(self.obs.sidereal_time())))
        self.xyz=nsew2xyz(antpos,self.lst,lat)
        self.nant=self.xyz.shape[0]
        self.freqs=freqs
        df=df
        nf=len(freqs)
        self.delays=np.arange(-nf/2,nf/2)/(nf*df)
        self.nvis=(self.nant-1)*self.nant/2
        self.model=np.zeros((len(freqs),self.nvis)).astype(complex)
        self.pcentre=np.array([self.lst,lat])
        if(self.driftMode):
            self.nunique,self.unique_map,self.unique_uvw,self.uvw=self.getUnique(self.pcentre)
#            print self.uvw
#            print len(freqs),self.nunique
            self.model_true=np.zeros((len(freqs),self.nunique)).astype(complex)
#        print self.unique_map
        self.modelspace='freq'
        self.data=np.zeros((len(freqs),self.nvis)).astype(complex)
        self.datastate='freq'
        #fifteen degrees at lowest frequency
        self.beam=beam
        self.bandpass=[Bandpass(flatband,[1]) for mm in range(self.xyz.shape[0])]#instantiate all antennae to same bandpass
               
        #return
    #************************************************************
    #compute gridded psf in uv space
    #************************************************************
    def computepsf_uv(self,nu,nv,du,dv,weight='natural'):
        nf=len(self.freqs)
        gridded_psf=np.zeros((nf,nv,nu))
        for mm in range(nf):
            sf=c/self.freqs[mm]
            gridded_psf[mm,:,:]=np.real(utils.griddata_nn_2D((self.uvw[:,:2]/sf).squeeze(),np.ones(self.uvw.shape[0]),du,dv,nu,nv))
            gridded_psf[mm,:,:]=gridded_psf[mm,:,:]+np.rot90(gridded_psf[mm,:,:],k=2)
        if(weight=='uniform'):
            gridded_psf[gridded_psf>0.]=gridded_psf[gridded_psf>0.]/gridded_psf[gridded_psf>0.]
        return gridded_psf

    #compute rotation synthesis psf quickly
    def computepsf_uv_fast():
        return

    def setBandpass(self,bandpass):
        self.bandpass=bandpass

    def setBeam(self,beam):
        self.beam=beam
    #************************************************************
    #give a different beam model
    #************************************************************
    def newbeam(self,beam):
        self.beam=beam
        return
    #************************************************************
    #update beam paramters
    #************************************************************
    def chgbeam(self,params):
        self.beam.set_params(params)
        return
    #new bandpass
    def newband(self,bandpass):
        self.bandpass=bandpass
        return
    def chgband(self,params):
        self.bandpass.params=bandpass
        return
    #************************************************************
    #return gain towards loc at frequency f
    #************************************************************
    def gain(self,loc,f):
        return self.beam.gain(loc,f)#*self.bandpass.gain(f)
    
    #************************************************************
    #phase to zenith
    #************************************************************
    def phasezenith(self):
        zenith=np.array([self.lst,self.lat])
        self.chgphasecentre(zenith)
        return
    #************************************************************
    #point primary beam
    #************************************************************
    def point(self,pcentre):
        if(self.beam.bfunc.__name__=='gaussbeam' or self.beam.bfunc.__name__=='tophatbeam' or self.beam.bfunc.__name__=='cosbeam' or self.beam.bfunc.__name__=='airybeam'):
            print 'pointing at '+str(pcentre)
            self.beam.chgparams([pcentre,self.beam.params[1]])
        elif(self.beam.bfunc.__name__=='mwabeam'):
            self.beam.chgparams([pcentre,self.obs])
        return
    #************************************************************
    #point beam to zenith
    #************************************************************
    def pointzenith(self):
        self.point(np.array([self.lst,self.lat]))

    #************************************************************
    #returns a list of uvw points for phase center pcentre with 
    #ra and dec
    #************************************************************
    def getuvw(self):
        return self.uvw
    #************************************************************
    #change date of interferometer
    #************************************************************
    def chgdate(self,date):
        self.obs.date=date
        self.lst=np.degrees(float(repr(self.obs.sidereal_time())))
        self.xyz=nsew2xyz(self.antpos,self.lst,self.lat)
        self.chgphasecentre(self.pcentre)
   #************************************************************
    #step iterate time by dt (seconds)
    #************************************************************
    def step(dt=2.):
        newdate=ephem.Date(self.obs.date+dt*ephem.second)
        self.chgdate(newdat)
    #************************************************************
    #radec is a 2-tuple for the phase centre of an interferometer
    #in degrees
    #xyz is a list of antenna or baseline coordinates in coord
    #system where Y is at -6h RA, Z is in the direction of the north pole
    #and X is in the 0h direction
    #returns uvw coordinates of antennae in coord system with 
    #phase centre at radec. 
    #uvw are still in whatever units xyz were in. 
    #************************************************************
    #************************************************************
    #phase array towards pcentre by updating uvw of physical antennae
    #does not recompute visibilities (must do in an additional step). 
    #************************************************************
    def chgphasecentre(self,pcentre):
        if(not(self.driftMode)):
            self.pcentre=pcentre
            pcentre=np.radians(pcentre)
            ra=pcentre[0]
            dec=pcentre[1]
            xyzrot=np.zeros(self.xyz.shape)
            xyzrot[:,0]=np.sin(ra)*self.xyz[:,0]+np.cos(ra)*self.xyz[:,1]
            xyzrot[:,1]=-np.sin(dec)*np.cos(ra)*self.xyz[:,0]+np.sin(dec)*np.sin(ra)*self.xyz[:,1]+np.cos(dec)*self.xyz[:,2]
            xyzrot[:,2]=np.cos(dec)*np.cos(ra)*self.xyz[:,0]-np.cos(dec)*np.sin(ra)*self.xyz[:,1]+np.sin(dec)*self.xyz[:,2]
            uvwpos=np.zeros((self.nvis,3))
            uvwpos[:,0]=utils.diffcomb(xyzrot[:,0])
            uvwpos[:,1]=utils.diffcomb(xyzrot[:,1])
            uvwpos[:,2]=utils.diffcomb(xyzrot[:,2])
            self.uvw=uvwpos
            self.unique_uvw=unique_rows(np.around(self.uvw,decimals=3))
            #        for mm in range(self.unique_uvw.shape[0]):
            #            self.redundancy_uvw[mm]=self.uvw[(self.uvw==self.unique_uvw[mm,:]).any(axis=1),:].shape[0]
            self.nunique=self.unique_uvw.shape[0]
#            print self.nunique
            self.model_true=np.zeros((len(self.freqs),self.nunique),dtype=complex)
            #        plt.scatter(uvwpos[:,0],uvwpos[:,1])
            #        plt.show()
            #        plt.close()
        else:
            self.pcentre=np.array([self.lst,self.lat])
            #print 'pcenter='+str(self.pcenter)
        return

    def getUnique(self,pcentre):
        pcentre=np.radians(pcentre)
        ra=pcentre[0]
        dec=pcentre[1]
        xyzrot=np.zeros(self.xyz.shape)
        xyzrot[:,0]=np.sin(ra)*self.xyz[:,0]+np.cos(ra)*self.xyz[:,1]
        xyzrot[:,1]=-np.sin(dec)*np.cos(ra)*self.xyz[:,0]+np.sin(dec)*np.sin(ra)*self.xyz[:,1]+np.cos(dec)*self.xyz[:,2]
        xyzrot[:,2]=np.cos(dec)*np.cos(ra)*self.xyz[:,0]-np.cos(dec)*np.sin(ra)*self.xyz[:,1]+np.sin(dec)*self.xyz[:,2]
        uvwpos=np.zeros((self.nvis,3))
        uvwpos[:,0]=utils.diffcomb(xyzrot[:,0])
        uvwpos[:,1]=utils.diffcomb(xyzrot[:,1])
        uvwpos[:,2]=utils.diffcomb(xyzrot[:,2])
        uvw=uvwpos
        unique_uvw=unique_rows(np.round(uvw,decimals=3))
        nunique=unique_uvw.shape[0]
        #generate a dictionary between each baseline and unique baseline
        unique_dictionary={}
        for mm in range(uvw.shape[0]):
            row=np.where((np.round(uvw[mm,:],decimals=3)==np.round(unique_uvw,decimals=3)).all(axis=1))[0][0]
            unique_dictionary[mm]=row
        return nunique,unique_dictionary,unique_uvw,uvw
    #compute model for single freq to be run in parallel
#    def setmodel_freq(self,

    #************************************************************
    #compute visiblity model from list of l,m,s values
    #only compute unique visibilities from sources
    #than apply bandpasses to unique visibilities to obtain full 
    #visibility set. 
    #also allow for gradients
    #************************************************************
    def setmodelcube(self,srclist,dl,gradients=None,residuals=None,smin=None):
#        print srclist[0].shape
        if(gradients is None):
            gradients=np.vstack([np.zeros(srclist[0].shape[0]),np.zeros(srclist[0].shape[0])]).T
        if(residuals is None):
            residuals=np.zeros(srclist[0].shape[0])
        if(smin is None):
            smin=0.
#        print gradients.shape
        self.modelspace='freq'
        #world coordinate system
        w=wcs.WCS(naxis=2)
        w.wcs.crpix=[1,1]
        w.wcs.cdelt=np.array([dl,dl]) #(set so pixcrd gives radians)
        w.wcs.ctype=["RA---SIN","DEC--SIN"]
        w.wcs.crval=[self.pcentre[0],self.pcentre[1]]
        num_cores=multiprocessing.cpu_count()/2
        cfactor=32.8*3.*1e9*3e5/(2*pi)
        lms=[srclist[mm][:,[0,1]] for mm in range(len(self.freqs))]
        slists=[srclist[mm][:,2]+residuals[mm] for mm in range(len(self.freqs))]
        gains=[np.zeros(srclist[mm].shape) for mm in range(len(self.freqs))]
        unique_uvw=[self.unique_uvw*self.freqs[mm]/c for mm in range(len(self.freqs))]
        #get max frequency mask
        sfilter=np.sqrt(np.sum(np.radians(lms[0])**2.,axis=1))<=1.
#        print lms[0]
#        print np.any(np.isnan(gradients.flatten()))
        gfilter=np.sqrt(np.sum((np.radians(lms[0])+gradients*cfactor/self.freqs.min()**2.)**2.,axis=1))>=1.
        radec0=w.wcs_pix2world(srclist[0][:,[0,1]]/dl+1,1)
        #plt.scatter(radec0[:,0],radec0[:,1])
        #plt.show()
        
        for mm,f in enumerate(self.freqs):
            #filter smin
            radec=w.wcs_pix2world(srclist[mm][:,[0,1]]/dl+1,1)
            lms[mm]=np.radians(lms[mm])
            grads=cfactor*gradients/(f*f)
            grads[gfilter,:]=0.
            lms[mm]=lms[mm]+grads
            gains[mm]=self.gain(radec,f)
            slists[mm]=slists[mm]*gains[mm]
            slists[mm]=slists[mm][sfilter]
            lms[mm]=lms[mm][sfilter]
            
        nx=np.sqrt(srclist[0].shape[0])
#        print nx
#        plt.imshow((gains[0]).reshape((nx,nx)))
#        plt.colorbar()
#        plt.show()
        model=np.array(Parallel(n_jobs=num_cores)(delayed(utils.computevisibilities)(unique_uvw[i],lms[i],slists[i]) for i in range(len(unique_uvw))))
#        self.model=np.array([processChan(uvw[i],lms[i],slists[i],i) for i in range(len(self.freqs))])
        self.model_true=model
        #now map to model
        modelInd=0
        for mm in range(self.nant):
            for nn in range(mm):
                self.model[:,modelInd]=self.bandpass[mm].gain(self.freqs)*np.conj(self.bandpass[nn].gain(self.freqs))*model[:,self.unique_map[modelInd]]
                modelInd+=1
        self.modelspace='freq'
        
        return

#    def refAmp(f,r,phi,tau):
#        return 1./np.sqrt(r*r+1.-2.*r*np.cos(2.*np.pi*f*tau+phi))
#    def refPha(f,r,phi,tau):
#        return np.arctan(r*np.sin(2.*np.pi*f*tau+phi)/(1.-r*np.cos(2*np.pi*f*tau+phi)))



 


    #************************************************************
    #compute visibility model from list ra,dec,s180,alpha
    #and pcentre(2-tuple)
    #for now assume spectrally flat
    #can easily do this fast in the future using shared memory 
    #arrays
    #gradients is the TEC gradient along the sight line to each source
    #in units of TECU/km
    #************************************************************
    def setmodel(self,srclist,gradients=None,additive=False):
        if(gradients is None):
            gradients=np.vstack([np.zeros(srclist.shape[0]),np.zeros(srclist.shape[0])]).T
       # print 'gradients'
       # print gradients
        self.modelspace='freq'
        #eliminate all sources that are below the horizon
        keepers=np.empty(srclist.shape[0],dtype=bool)
        keepers[:]=True
        for mm in range(srclist.shape[0]):
            pos=ephem.FixedBody()
            pos._ra=np.radians(srclist[mm,0])
            pos._dec=np.radians(srclist[mm,1])
            pos.compute(self.obs)
            if(float(repr(ephem.degrees(pos.alt))<0.)):
                keepers[mm]=False
            else:
                keepers[mm]=True
        srclist=srclist[keepers,:]
        #create world coordinate system
        w=wcs.WCS(naxis=2)
        w.wcs.crpix=[1,1]
        w.wcs.cdelt=np.array([1.,1.])*(180./pi) #(set so pixcrd gives radians)
        w.wcs.ctype=["RA---SIN","DEC--SIN"]
        w.wcs.crval=[self.pcentre[0],self.pcentre[1]]
#        print w.wcs.crval
        #convert all ra/dec to lm
 #       print srclist[:,[0,1]]
        lm=w.wcs_world2pix(srclist[:,[0,1]],0)
  #      print 'lm='+str(lm)
        #remove any sources that will go out of fov
        cfactor=2.8*3*1e9*3e5 #factor for converting TECU/km to phase in radians/km
        popout=np.sum((lm+gradients*cfactor/self.freqs.min()**2.)**2.,axis=1)>1.
        gradients[popout,:]=0.
        num_cores=multiprocessing.cpu_count()/2
        model=np.array(Parallel(n_jobs=num_cores)(delayed(utils.computevisibilities)(self.unique_uvw*self.freqs[i]/c,lm+gradients*cfactor/self.freqs[i]**2.,srclist[:,2]*(self.freqs[i]/180e6)**-srclist[:,3]*self.gain(srclist[:,[0,1]],self.freqs[i])) for i in range (len(self.freqs))))
        if(additive):
            self.model_true+=model
        else:
            self.model_true=model
        #now map to model
        for mm in range(self.nant):
            for nn in range(mm):
                modelInd=mm*self.nant+nn
                if(additive):
                    model[:,modelInd]+=self.bandpass[mm].gain(self.freqs)*np.conj(self.bandpass[nn].gain(self.freqs))*model[:,self.unique_map[modelInd]]
                else:
                    model[:,modelInd]=self.bandpass[mm].gain(self.freqs)*np.conj(self.bandpass[nn].gain(self.freqs))*model[:,self.unique_map[modelInd]]
        self.modelspace='freq'
        return
    



    #************************************************************
    #take delay transform of all visibilities. 
    #************************************************************
    def delaytransformmodel(self):
        if(self.modelspace=='freq'):
            df=self.freqs[1]-self.freqs[0]
            self.model=np.fft.ifftshift(np.fft.fft(np.fft.fftshift(self.model,axes=0),axis=0),axes=0)*df
            self.modelspace='delay'
        return
    #************************************************************
    #take inverse delay transform
    #************************************************************
    def idelaytransformmodel(self):
        if(self.modelspace=='delay'):
            df=self.freqs[1]-self.freqs[0]
            nf=self.freqs.size
            dd=1./(df*nf)
            self.model=np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(self.model,axis=0),axis=0),axis=0)*dd*nf#factor of nf added to compensate for numpy ifft normalization convention. 
            self.modelspace='freq'
        return

    #************************************************************
    #get noise power spectrum on delay modes in mK^2 Mpc^3 between
    #at frequency fc
    #tsys150 is the system temperature at 150 MHz
    #************************************************************
    def delayNoise(self,tsys150,tint,fc):
        redundancy=np.zeros(self.uvw.shape[0])
        #compute redundancy factor at fc
        for mm in range(self.uvw.shape[0]):
            redundancy[mm]=self.uvw[(np.around(self.uvw,decimals=3)==self.uvw[mm,:]).all(axis=1),:].shape[0]
        #print np.sum(redundancy)
        return cosmology.X(fc)**2.*cosmology.Y(fc)*tsys150**2.*(fc/150e6)**(-2.5*2.)*self.beam.beamP(fc)**2./self.beam.beamPP(fc)/(2*tint*redundancy)
        
#test interferometer class with simple E-W interferometer at the equator and a source at zenith
"""
xlist=np.arange(0,10,.1)
ylist=np.zeros(xlist.size)
zlist=np.zeros(xlist.size)
freqs=np.array([[180e6]])
xlist=xlist*c/freqs[0,0]
xyzlist=np.vstack([xlist,ylist,zlist]).T
print xyzlist.shape
lon=0.
lat=0.
myint=interferometer(lon,lat,xyzlist,freqs)
print myint.lst
#create a source at interferometer LST
lst=myint.lst
dl=10.
dm=0.
srclist=np.array([[lst+dl,0.+dm,1.,0.]])
print srclist[0,:2]
myint.chgcentre([lst,0.])
myint.setmodel(srclist)
print myint.uvw.shape
print myint.model.shape
#print myint.xyz
#print myint.uvw
#print myint.uvw
#print myint.model
dl=np.sin(dl*pi/180.)
dm=np.sin(dm*pi/180.)
uaxis=np.linspace(0,10.,1000)
model=np.cos(2*pi*(dl*uaxis))/np.sqrt(1.-dl**2.-dm**2.)
plt.plot(myint.uvw[0,:,0],np.real(myint.model[0,:]),'ko',uaxis,model)
plt.show()
plt.close()
"""
"""
#test time evolution with two element E-W interferometer phased to the north pole
xlist=np.arange(0,2)
ylist=np.zeros(xlist.size)
zlist=np.zeros(xlist.size)
"""
#************************************************************
#returns an interferometer with antpos and freqs at the MRO
#will be instantiated with date set to now
#************************************************************
def getMRO(antpos,freqs):
    obs=Interferometer(116.671,-26.701,antpos,freqs)
    return obs
#def getKahroo():
#    return kru


