#!/usr/bin/env python
import numpy as np
import geometry as geom
from mwapy import ephem_utils
import ephem
import scipy.special as special
import scipy.integrate as integrate
import math
from mwapy.pb import primary_beam as mwa_primary_beam
pi=np.pi
c=3e8
#load sweetspots for mwa beam
mwaSweetSpots=np.loadtxt('sweetspots.txt',skiprows=2)
mwaSweetSpots[:,1]=mwaSweetSpots[:,1]+90.
#************************************************************
#simple beam centred at pcentre. cpc and csrc are sky coords
#************************************************************
def tophatbeam(loc,f,params):
    cpc=params[0]
    afwhm=params[1]
    if(len(loc.shape)==1):
        loc=np.array([loc])
    output=np.ones(loc.shape[0])
    output[geom.arc_length_list(loc,cpc)>afwhm]=0.
    return output

#************************************************************
#simple gaussian beam centered at pcentre with fwhm at 180MHz
#************************************************************
def gaussbeam(loc,f,params):
    if(len(loc.shape)==1):
        loc=np.array([loc])
    cpc=params[0]
    afwhm=params[1]
    sep=geom.arc_length_list(loc,cpc)
    sigma=afwhm/np.sqrt(2.*np.log(2.))#*(180e6/f)
    return np.exp(-.5*(sep/sigma)**2.)
#************************************************************
#an airy disk beam
#make it achromatic with three params
#************************************************************
def airybeam(loc,f,params):
    if(len(loc.shape)==1):
        loc=np.array([loc])

    a=params[1]/2#params 1 is the aperture DIAMETER
    if(len(params)==3):
        fc=params[2]
        k=2*pi*fc/c
    else:
        k=2*pi*f/c
    cpc=params[0]
    sep=np.radians(geom.arc_length_list(loc,cpc))
    x=k*np.sin(sep)*a
    output=(2*special.jn(1,x)/x)**2.
    output[sep==0.]=1.
    return output
                  
    

#************************************************************
#a simpl sin beam. we keep it achromatic due to sharp cutoff
#************************************************************
def cosbeam(loc,f,params):
    if(len(loc.shape)==1):
        loc=np.array([loc])
    cpc=params[0]
    afwhm=params[1]#*180e6/f
    #sep=geom.arc_length_list(loc,cpc)
    sep1=np.abs(loc[:,0]-cpc[0])*np.cos(np.radians(cpc[1]))
    sep2=np.abs(loc[:,1]-cpc[1])
    output=np.cos(2.*pi*sep1/afwhm)*np.cos(2.*pi*sep2/afwhm)
    output[np.logical_or(sep1>afwhm/2.,sep2>afwhm/2.)]=0.
    return output
#************************************************************
#params include the array ephem observatory and the pointing center
#which is translated into an az and el and delay setting. 
#************************************************************

def mwabeam(loc,f,params):
    #az,alt=obs.az,el(ra,dec)
    obs=params[1]
    pcentre=params[0]
    #find az el of pointing centre
    pos=ephem.FixedBody()
    pos._ra=np.radians(pcentre[0])
    pos._dec=np.radians(pcentre[1])
    pos.compute(obs)
    az=float(repr(ephem.degrees(pos.az)))*180./np.pi
    alt=float(repr(ephem.degrees(pos.alt)))*180./np.pi
    #find closets sweet splot
    cpc=[az,alt]
    #    print 'center='+str(cpc)
    sep=geom.arc_length_list(mwaSweetSpots[:,:2],cpc)
    #print sep
    delays=mwaSweetSpots[np.where(sep==sep.min())[0][0],3:].astype(int)
    ra=loc[:,0]
    dec=loc[:,1]
    #azAlt=np.zeros(loc.shape)
    pos=ephem.FixedBody()
#    for ss in range(azAlt.shape[0]):
#        pos._ra=np.radians(loc[ss,0])
#        pos._dec=np.radians(loc[ss,1])
#        pos.compute(obs)
#        az=float(repr(ephem.degrees(pos.az)))
#        alt=float(repr(ephem.degrees(pos.alt)))
#        azAlt[ss,:]=az,alt
    lst=np.degrees(float(repr(obs.sidereal_time())))
    #print 'lat='+str(float(repr(obs.lat))*180./np.pi)
    lat=float(repr(obs.lat))*180./np.pi
    #print 'delays='+str(delays)
    az,alt=ephem_utils.eq2horz(lst-ra,dec,lat)
    theta=(90.-alt)*math.pi/180.
    phi=az*math.pi/180.
    dip_sep=mwa_primary_beam._DIPOLE_SEPARATION
    dipheight=mwa_primary_beam._DIPOLE_HEIGHT
    rX,rY=mwa_primary_beam.MWA_Tile_analytic(theta,phi,freq=f,delays=delays,dipheight=dipheight,dip_sep=dip_sep,zenithnorm=True,power=True)
    #for now let's return the stokes sum
    return (rX+rY)/2.

def heraBeam(loc,f,params):
    obs=params[0]
    df=params[1]
    fc=params[2]
    dPhi=np.radians(params[3])
    dTheta=np.radians(params[4])
    beamCube=params[5]
    pos=ephem.FixedBody()
    lst=np.degrees(float(repr(obs.sidereal_time())))
    lat=float(repr(obs.lat))*180./np.pi
    ra=loc[:,0]
    dec=loc[:,1]
    az,alt=ephem_utils.eq2horz(lst-ra,dec,lat)
    theta=(90.-alt)*math.pi/180.
    phi=az*math.pi/180.
    thInd=np.round(theta/dTheta).astype(int)
    phInd=np.round(phi/dPhi).astype(int)
    output=np.zeros(len(ra))
    fInd=np.round((f-fc)/df)+beamCube.shape[0]/2
    inFieldTh=np.logical_and(thInd>=0,thInd<beamCube.shape[2])
    inFieldPh=np.logical_and(phInd>=0,phInd<beamCube.shape[1])
    inField=np.logical_and(inFieldTh,inFieldPh)
    if(fInd>=0 and fInd<beamCube.shape[0]):
        output[inField]=beamCube[fInd,phInd[inField],thInd[inField]]
    return output



#************************************************************
#represents the beam of an antennae
#accepts a function bfunc with the form
#bfunc(location,f,params)
#params are any extra parameters that bfunc might need
#for example, the MWA beam would need delay settings while
#the gaussian beam would take a fwhm at a certain frequency
#and the pointing centre of the beam 
#************************************************************
class Beam():
    def __init__(self,bfunc,params):
        self.bfunc=bfunc
        self.params=params
    def chgparams(self,params):
        self.params=params
    def gain(self,loc,f):
        return self.bfunc(loc,f,self.params)
    #get integrated beam power
    def beamP(self,f):
    #get the integrated beam power
        if(self.bfunc.__name__=='tophatbeam'):
            afwhm=self.params[1]
            return pi*(np.radians(afwhm/2))**2.
        elif(self.bfunc.__name__=='gaussbeam'):
            afwhm=self.params[1]
            sigma=np.radians(afwhm/np.sqrt(2.*np.log(2.)))#*(180e6/f))
            return pi*sigma*sigma/2.
        elif(self.bfunc.__name__=='cosbeam'):
            afwhm=self.params[1]
            return np.radians(afwhm)**2./pi**2.
        elif(self.bfunc.__name__=='airybeam'):
            a=self.params[1]
            k=2*pi*f/c
            return 4*pi/(k*a)*4.*integrate.quad(lambda x: (special.jn(1,np.sqrt(1-x**2.))/np.sqrt(1-x**2.))**2.,0,1)[0]
        
        elif(self.bfunc.__name__=='mwabeam'):
            #integrate mwabeam over theta and phi
              obs=self.params[1]
              pcentre=self.params[0]
              pos=ephem.FixedBody()
              pos._ra=np.radians(pcentre[0])
              pos._dec=np.radians(pcentre[1])
              pos.compute(obs)
              az=float(repr(ephem.degrees(pos.az)))*180./np.pi
              alt=float(repr(ephem.degrees(pos.alt)))*180./np.pi
              cpc=[az,alt]
              sep=geom.arc_length_list(mwaSweetSpots[:,:2],cpc)
              delays=mwaSweetSpots[np.where(sep==sep.min())[0][0],3:].astype(int)
              dip_sep=mwa_primary_beam._DIPOLE_SEPARATION
              dipheight=mwa_primary_beam._DIPOLE_HEIGHT
              def mwaBeamI(theta,phi):
                  rX,rY=mwa_primary_beam.MWA_Tile_analytic(theta,phi,freq=f,delays=delays,dipheight=dipheight,dip_sep=dip_sep,zenithnorm=True,power=True)
                  return (rX+rY)/2.*np.sin(theta)
              return integrate.dblquad(mwaBeamI,0.,np.pi/2.,lambda x: 0., lambda x: 2.*np.pi)[0]
        elif(self.bfunc.__name__=='heraBeam'):
            fc=self.params[2]
            df=self.params[1]
            dPhi=np.radians(self.params[3])
            dTheta=np.radians(self.params[4])
            beamCube=self.params[5]
            fInd=np.round((f-fc)/df)+beamCube.shape[0]/2
            if(fInd>=0 and fInd<beamCube.shape[0]):
                beamCube=beamCube*dPhi*dTheta
                _,_,thetaCube=np.meshgrid(np.zeros(beamCube.shape[1]),np.zeros(beamCube.shape[0]),np.arange(beamCube.shape[2])*dTheta)
                thetaCube=np.sin(thetaCube)
                return np.sum((beamCube[fInd,:,:]*thetaCube[fInd,:,:]).flatten())
            

                

    def beamPP(self,f):
        if(self.bfunc.__name__=='tophatbeam'):
            return self.beamP(f)
        elif(self.bfunc.__name__=='gaussbeam'):
            return self.beamP(f)/2.
        elif(self.bfunc.__name__=='cosbeam'):
            afwhm=self.params[1]
            return np.radians(afwhm)**2./16.
        elif(self.bfunc.__name__=='airybeam'):
            a=self.params[1]
            k=2*pi*f/c
            return 4*pi/(k*a)*16.*integrate.quad(lambda x: (special.jn(1,np.sqrt(1-x**2.))/np.sqrt(1-x**2.))**4.,0,1)[0]

        elif(self.bfunc.__name__=='mwabeam'):
            #integrate mwabeam over theta and phi
            obs=self.params[1]
            pcentre=self.params[0]
            pos=ephem.FixedBody()
            pos._ra=np.radians(pcentre[0])
            pos._dec=np.radians(pcentre[1])
            pos.compute(obs)
            az=float(repr(ephem.degrees(pos.az)))*180./np.pi
            alt=float(repr(ephem.degrees(pos.alt)))*180./np.pi
            cpc=[az,alt]
            sep=geom.arc_length_list(mwaSweetSpots[:,:2],cpc)
            delays=mwaSweetSpots[np.where(sep==sep.min())[0][0],3:].astype(int)
            dip_sep=mwa_primary_beam._DIPOLE_SEPARATION
            dipheight=mwa_primary_beam._DIPOLE_HEIGHT
            def mwaBeamI(theta,phi):
                rX,rY=mwa_primary_beam.MWA_Tile_analytic(theta,phi,freq=f,delays=delays,dipheight=dipheight,dip_sep=dip_sep,zenithnorm=True,power=True)
                return ((rX+rY)/2.)**2.*np.sin(theta)
            return integrate.dblquad(mwaBeamI,0.,np.pi/2.,lambda x: 0., lambda x: 2*np.pi)[0]
        elif(self.bfunc.__name__=='heraBeam'):
            fc=self.params[2]
            df=self.params[1]
            dPhi=np.radians(self.params[3])
            dTheta=np.radians(self.params[4])
            beamCube=self.params[5]
            fInd=np.round((f-fc)/df)+beamCube.shape[0]/2
            if(fInd>=0 and fInd<beamCube.shape[0]):
                beamCube=beamCube**2.*dPhi*dTheta
                _,_,thetaCube=np.meshgrid(np.zeros(beamCube.shape[1]),np.zeros(beamCube.shape[0]),np.arange(beamCube.shape[2])*dTheta)
                thetaCube=np.sin(thetaCube)
                return np.sum((beamCube[fInd,:,:]*thetaCube[fInd,:,:]).flatten())
            
