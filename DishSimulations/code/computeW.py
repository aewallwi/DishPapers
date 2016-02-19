#!/usr/bin/env python

import numpy as np
import sys
import scipy.interpolate as interp
import re
import scipy.optimize as op
import matplotlib.pyplot as plt
import cosmology as cosmo

print 'cosmo=%s'%(cosmo.__file__)

littleh=0.68
config,foregrounds=sys.argv[1],sys.argv[2]

def pmLogPlot(y,ymin):
    posVals=y>0.
    negVals=y<0.
    
    



#numerical recipes 2nd edition page 661
def linFit(x,y,dy):
    s=np.sum(1./dy**2.)
    sx=np.sum(x/dy**2.)
    sy=np.sum(y/dy**2.)
    sxx=np.sum(x**2./dy**2.)
    sxy=np.sum(x*y/dy**2.)
    delta=s*sxx-(sx)**2.
    a=(sxx*sy-sx*sxy)/delta
    b=(s*sxy-sx*sy)/delta
    sa=sxx/delta
    sb=s/delta
    csab=-sx/delta
    covMat=np.array([[sb,csab],[csab,sa]])
    return np.array([b,a]),covMat


lowInd=4
highInd=2
def readConfig(fileName):
    """
    read parser configuration file

    """

    # read pipeline parameters
    print "=" * 50
    print "Read configuration file " + fileName
    config = {}
    with open(fileName, "r") as cfile:
        for line in cfile:
            if line[0] != "#":
                params = line.split("#")[0].split("=")  # split off comments
                print params
                config[params[0]] = eval(params[1])
    return config

if(foregrounds=='mod'):
    postfix='Mod'
elif(foregrounds=='opt'):
    postfix='Opt'
elif(foregrounds=='pess'):
    postfix='Pes'

config=readConfig(config)
steps=np.array(config['STEPFRACS'])
dStep=config['DSTEP']
nSteps=len(steps)
reStep=re.compile('step_[0-9]{1,2}')
reNum=re.compile('[0-9]{1,2}')
reGHz=re.compile('[0-9]{1,2}.[0-9]{5}')
dirname=config['STEPPARAM']+'_'+config['LABEL']+'_Output/data/'
plotdir=config['STEPPARAM']+'_'+config['LABEL']+'_Output/plots/'
trackMode=config['TRACKMODE']
if(trackMode=='TRACK'):
    track=trackMode+config['TRACK']+'hrs'
else:
    track=trackMode

senseList=dirname+'senseList'+postfix+'.txt'
print senseList
senseList=open(senseList).readlines()
dDeltaFile=config['STEPPARAM']+'_'+config['LABEL']+'_Output/data/dDelta_'+config['STEPPARAM']+'.npz'
derivMat=np.load(dDeltaFile)

zAxis=derivMat['zAxis']
kAxis=derivMat['kAxis']
nK=len(kAxis)
nZ=len(zAxis)


wMat=np.zeros((nZ,nK))
wMat1=np.zeros((nZ,nK))
wMat2=np.zeros((nZ,nK))
errMat=np.zeros((nZ,nK))

senseMat=[[ [] for nn in range(nZ)] for mm in range(nSteps)]


#load up files
for mm,fileName in enumerate(senseList):
    print 'filename=%s'%(fileName)
    step=reStep.findall(fileName)[0]
    stepNum=int(reNum.findall(step)[0])
    freq=float(reGHz.findall(fileName)[0])
    #get z
    z=np.round(cosmo.f2z(freq*1e9),decimals=2)
    print 'z=%f'%(z)
    print zAxis
    zInd=np.where(zAxis==np.round(z,decimals=2))[0][0]
    senseFile=np.load(dirname+fileName[:-1])
    kAx=derivMat['kAxis']/littleh
    print senseFile['ks']
    minK=senseFile['ks'].min()
    maxK=senseFile['ks'].max()
    dp21=interp.interp1d(senseFile['ks'],senseFile['errs'])
#    print kAx[lowInd:-highInd]
    #interpolate to ks of 21cmFAST sims
    inRange=np.logical_and(kAx>minK,kAx<maxK)
    senseTemp=np.ones(len(kAx))*np.inf
    senseTemp[inRange]=dp21(kAx[inRange])
    senseMat[stepNum][zInd]=senseTemp
senseMat=np.array(senseMat)
print senseMat.shape
ddPsList=np.zeros((nZ,nK))
print 'nz=%f'%(nZ)
print 'nk=%f'%(nK)
#now compute derivative of senseMat
print senseMat.shape
for mm in range(nZ):
    for nn in range(nK):
        y=senseMat[:,mm,nn]
        fitFlags=np.logical_and(np.abs(steps)<dStep,np.invert(np.isnan(y)))
        xFit=steps[fitFlags]
        yFit=y[fitFlags]
        if(len(xFit)>=2):
            params,dparams=linFit(xFit,yFit,np.ones(len(xFit)))
            #op.curve_fit(lambda x,*p:x*p[0]+p[1],xFit,yFit,p0=[0.,0.],maxfev=int(1e4))
#            print params
#            print dparams
            if(not(np.any(np.isinf(dparams)))):
                plt.plot(xFit,yFit,'bo',label='data')
                plt.plot(xFit,xFit*params[0]+params[1],'-r',label='fit')
                plt.savefig(config['STEPPARAM']+'_'+config['LABEL']+'_Output/plots/dSigma_k%d_z%d_%s_%s.png'%(nn,mm,config['ARRAY'],foregrounds))
                plt.close()
#                if(np.abs(params[0])>3.*np.sqrt(dparams[0,0])):
 #                   slope=params[0]
  #              else:
   #                 slope=0.
                slope=params[0]
            else:
                slope=0.
        else:
            slope=0.
        ddPsList[mm,nn]=slope

#save derivatives
np.savez(config['STEPPARAM']+'_'+config['LABEL']+'_Output/data/ddData_%s_%s.npz'%(config['ARRAY'],foregrounds),ddGrid=ddPsList,kAxis=kAxis,zAxis=zAxis)
#now go through senseList but only append if step is zero

for mm,fileName in enumerate(senseList):
    step=reStep.findall(fileName)[0]
    stepNum=int(reNum.findall(step)[0])
    freq=float(reGHz.findall(fileName)[0])
    #get z
    z=np.round(cosmo.f2z(freq*1e9),decimals=2)
    zInd=np.where(zAxis==np.round(z,decimals=2))[0][0]

    if(steps[stepNum]==0.):
        senseFile=np.load(dirname+fileName[:-1])
        kAx=derivMat['kAxis']/littleh#divide kaxis by little h
        wLine=derivMat['dGrid'][zInd,:]
#        print 'kaxis='+str(kAxis)
#        print 'ks='+str(senseFile['ks'])
        dp21=interp.interp1d(senseFile['ks'],senseFile['errs'])
        minK=senseFile['ks'].min()
        maxK=senseFile['ks'].max()
        inRange=np.logical_and(kAx>minK,kAx<maxK)
#        print senseFile['errs'].shape
        wMat1[zInd,inRange]=((wLine[inRange])/dp21(kAx[inRange]))
        wMat2[zInd,inRange]=np.sqrt(2)*((ddPsList[zInd,inRange])/dp21(kAx[inRange]))
        errMat[zInd,inRange]=dp21(kAx[inRange])
        errMat[zInd,np.invert(inRange)]=np.inf

#now run through z and k and make plots 

for nn in range(nK):
#    print wMat1[:,nn]
#    print wMat2[:,nn]
#    print np.any(np.abs(wMat1[:,nn])>0.)
#    print np.any(np.abs(wMat2[:,nn])>0.)
    if(np.any(np.abs(wMat1[:,nn])>0.) and np.any(np.abs(wMat2[:,nn])>0.)):
        plt.plot(zAxis,np.abs(wMat1[:,nn]),'-r',label='mean term')
        plt.plot(zAxis,np.abs(wMat2[:,nn]),'-k',label='variance term')
        plt.title('k=%f hMpc$^{-1}$'%(kAx[nn]))
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig(plotdir+'wTerms_%s_k%d.png'%(foregrounds,nn))
        plt.close()

wMat=np.array([wMat1,wMat2])
                 
                 

np.savez(config['STEPPARAM']+'_'+config['LABEL']+'_Output/data/wMatrix_%s_%s_%s_%s.npz'%(config['STEPPARAM'],config['ARRAY'],foregrounds,track),kAxis=derivMat['kAxis'],zAxis=derivMat['zAxis'],wMat=wMat,errMat=errMat)
np.savez(config['STEPPARAM']+'_'+config['LABEL']+'_Output/data/wMatrix1_%s_%s_%s_%s.npz'%(config['STEPPARAM'],config['ARRAY'],foregrounds,track),kAxis=derivMat['kAxis'],zAxis=derivMat['zAxis'],wMat=wMat1,errMat=errMat)
np.savez(config['STEPPARAM']+'_'+config['LABEL']+'_Output/data/wMatrix2_%s_%s_%s_%s.npz'%(config['STEPPARAM'],config['ARRAY'],foregrounds,track),kAxis=derivMat['kAxis'],zAxis=derivMat['zAxis'],wMat=wMat2,errMat=errMat)


