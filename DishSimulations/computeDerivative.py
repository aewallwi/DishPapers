 #May 14th 2015
#Aaron Ewall-Wice
#aaronew@mit.edu
#read config file
#get parameter, its fiducial value and it's fractional change. 
#for each change, select 21cmFAST run and compute power spectrum for all datacubes
import numpy as np
import os
import sys
import re
import glob
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import multiprocessing
import powerSpectrum as power 
import scipy.optimize as op
import scipy.interpolate as interp
num_cores=multiprocessing.cpu_count()

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

#get list of data cubes 
def getCubeNames(dirName):
    fileList=glob.glob(dirName+'/dTb/delta_T')
    #get redshifts form file list
    return fileList
        

 #read in data cube,return its side length as well. 
def readCube(cubeName):
    print 'cubename=%s'%(cubeName)
    #get side length from cube name
    sideFind=re.compile('_[0-9]{2,4}_')
    lenFind=re.compile('[0-9]{2,4}Mpc')
    zFind=re.compile('z[0-9]{3}.[0-9]{2}')
    nside=int(sideFind.findall(cubeName)[0][1:-1])
    dSide=np.round(float(lenFind.findall(cubeName)[0][:-3]))
    z=float(zFind.findall(cubeName)[0][1:])
    dCube=np.fromfile(cubeName,dtype='float32').reshape(nside,nside,nside)
    return z,dSide,dCube

#once I have the list, I will want to compute power spectra for each cube       
#round to nearest z
def computeSpectra(dName,kmin,kmax,nk,zListIn):
    print zListIn
    cubeList=glob.glob(dName+'/dTb/delta_T*')
    zList=[]
    psList=[]
    kaxList=[]
    #first generate list of z
    for fname in cubeList:
        zRe=re.compile('z[0-9]{3}.[0-9]{2}')
        z=float(zRe.findall(fname)[0][1:])
        zList.append(z)
    #find dz from zListIn
    zList=np.array(zList)
    dz=zListIn[1]-zListIn[0]
    zInds=np.ones(len(zList))*-99#np.round((zList-zListIn[0])/dz).astype(int)
    #discard duplicates
    #for each z in zListIn, find clost z and set zInd to that closest z
    cubeListUse=[]
#    for mm in range(len(zListIn)):
#        if(zListIn[mm]>=zList[0] and zListIn[mm]<zList[-1]):
#            useInd=np.where(np.abs(zListIn[mm]-zList)==np.abs(zListIn[mm]-zList).min())
#            #useInd=0
#            print useInd
#            useInd=useInd[0][0]
#            cubeListUse.append(cubeList[useInd])
            
        
        
    zList=[]
    #now only compute cubes where z has a valid zInd         
    for fname in cubeList:
        zRe=re.compile('z[0-9]{3}.[0-9]{2}')
        z=float(zRe.findall(fname)[0][1:])
        if(z>=zListIn[0]-2. and z<zListIn[-1]+2.):
            z,dSide,dCube=readCube(fname)
            kaxis,counts,ps=power.ps1D(dCube,dSide/dCube.shape[0],kmin,kmax,nk,logbins=True)
            ps=ps*kaxis**3./(2*np.pi)#convert to delta squared
            psList.append(ps)
            zList.append(z)
            kaxList.append(kaxis)
    return zList,kaxList,psList

#given a list of values compute power spectra
def computeSteps(configFile):
    config=readConfig(configFile)
    param=config['STEPPARAM']
    kmin=config['KMIN']
    kmax=config['KMAX']
    dlgk=config['DLKBIN']
    steps=config['STEPFRACS']
    zmin=config['ZMIN']
    zmax=config['ZMAX']
    label=config['LABEL']
    dstep=config['DSTEP']
    #compute dlogk
    kmin=np.log10(kmin)
    kmax=np.log10(kmax)
    nk=int(np.round((kmax-kmin)/dlgk))
    kmin=10.**kmin
    kmax=10.**kmax
    #load up file lists
    dirList=[]
    nstep=len(steps)
    for mm in range(nstep):
        dname='21cmFAST_%s_step%d'%(param,mm)
        dirList.append(dname)
    #now compute power spectra
    bandInfo=np.loadtxt(param+'_'+label+'_Output/bandInfo.txt')
    zList=bandInfo[:,0]
    psLists=Parallel(n_jobs=num_cores)(delayed(computeSpectra)(dirList[mm],kmin,kmax,nk,zList) for mm in range(len(dirList)))
    #make plots for each redshift. 
    #for each k bin, plot line for each step at each resdhift. 
    #create a standard z list that will be interpolated
    #zAxis=np.arange(zmin,zmax,0.2)
    #create zAxis from intersection of all z axes
    zAxis=psLists[0][0]
    for mm in range(1,len(psLists)):
        zAxis=np.intersect1d(zAxis,psLists[mm][0])
    print 'zAxis=%s'%zAxis
    #iterate through power spectra and interpolate 
    kMatrix=np.zeros((len(zAxis),nk))
    psMatrix=np.zeros((len(zAxis),nk,nstep))
    
    
    for mm in range(nk):
        for nn in range(nstep):
            print mm
            print nn
            print np.array(psLists[nn][2]).shape
            ps=np.array(psLists[nn][2])[:,mm]
            k=np.array(psLists[nn][1])[:,mm]
            zAxisTemp=np.array(psLists[nn][0])
            #only store z
            print 'k=%s'%(k)
            kMatrix[:,mm]=k[0]
            for pp,z in enumerate(zAxis):
                print np.where(zAxisTemp==z)
                print z
                print ps
                psMatrix[pp,mm,nn]=ps[np.where(zAxisTemp==z)[0][0]]

    print psLists[0]
    #plot power spectra in steps
    nz=len(zAxis)
    
    outDir=os.getcwd()+'/'+param+'_'+label+'_Output'
    plotDir=outDir+'/plots'
    dataDir=outDir+'/data'
    if(not(os.path.isdir(outDir))):
        os.mkdir(outDir)
    if(not(os.path.isdir(plotDir))):
        os.mkdir(plotDir)
    if(not(os.path.isdir(dataDir))):
        os.mkdir(dataDir)
    
    psmatOut=np.zeros((len(zList),psMatrix.shape[1],psMatrix.shape[2]))
    #save power spectra in 21cmFAST format for sensitivity calculation. 
    #interpolate power spectra to zListIn redshifts
    for mm in range(nstep):
        for nn in range(psMatrix.shape[1]):
            psvec=psMatrix[:,nn,mm]
            psFunc=interp.interp1d(zAxis,psvec)
            psmatOut[:,nn,mm]=psFunc(zList)
            
    for nn in range(len(zList)):
        for mm in range(nstep):
            kvec=kMatrix[nn,:]
            psvec=psmatOut[nn,:,mm]
            z=zList[nn]
            #for each k interpolate to zListIn redshifts
            psmat=np.vstack([kvec,psvec]).T
            np.savetxt(dataDir+'/ps_%s_step%d_z%.2f.txt'%(param,mm,z),psmat,fmt='%e')

            


    #compute derivatives and perform a fit on inner 10 %
    #create a derivative grid as a function of kbin and z
    #derivGrid=np.zeros((nz,nk))
    for nn in range(nz):
        plt.close()
        for mm in range(len(psLists)):
#            print steps[mm]
            if(not(np.all(np.abs(psMatrix[nn,:,mm])<1e-5)) and np.all(psMatrix[nn,:,mm]>0.)):
                plt.plot(kMatrix[nn,:],psMatrix[nn,:,mm],label='%s=%.1e'%(param,(1+steps[mm])))
#            plt.plot(psLists[mm][1][nn][:],psLists[mm][2][nn][:],label='%s=%.1e'%(param,(1+steps[mm])))
#            plt.text(psLists[mm][1][nn][0],psLists[mm][2][nn][0],'%s=%.1e'%(param,(1+steps[mm])),fontsize=8,ha='right')
            
                plt.text(kMatrix[nn,0],psMatrix[nn,0,mm],'%s=%.1e'%(param,(1+steps[mm])),fontsize=8,ha='right')
#            print kMatrix[nn,:]
#            print psMatrix[nn,:,mm]
        
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('k (Mpc)$^{-1}$')
        plt.ylabel('$\\Delta^2(k)$ (mK$^2$)')
        z=zAxis[nn]
        plt.title('z=%.2f'%(z))
#        print plotDir+'/pS_z%05.2f_%s.png'%(z,param)
        plt.savefig(plotDir+'/pS_z%.2f_%s.png'%(z,param))
        plt.close()
    derivGrid=np.zeros((nz,nk))
    for mm in range(nk):
        #dgrid is an nstep x nz matrix
        dGrid=np.zeros((nstep,nz))
        for nn in range(nz):
            for oo in range(nstep):
                dGrid[oo,nn]=psMatrix[nn,mm,oo]#psLists[oo][2][nn][mm]
        plt.close()
        if(nstep==1):
            indZero=0
        else:
            indZero=steps.index(0)
#        print 'zeroInd=%d'%(indZero)
        if(len(steps)>1):
            for nn in range(nz):
                z=zAxis[nn]
            #fit line to +/- ten percent
                fitInds=np.abs(steps)<=dstep
                fitSteps=np.array(steps)[fitInds]
                fitPS=dGrid[fitInds,nn]
                fitSteps=fitSteps[np.invert(np.abs(fitPS)<1e-5)]
                fitPS=fitPS[np.invert(np.abs(fitPS)<1e-5)]
#                print fitPS
#                print fitSteps
                if(len(fitSteps)>2):
                    params,dparams=linFit(fitSteps,fitPS,np.ones(len(fitPS)))
            
                    if(not(np.any(np.isinf(dparams)))):
                        slope=params[0]
            #derivative is params[1]
                    else:
                        slope=0.
                    derivGrid[nn,mm]=slope
                    plt.close()
                    plt.plot(fitSteps,fitPS,'bo',label='simulations')
                    plt.plot(fitSteps,params[0]*fitSteps+params[1],'-r',label='fit')
                    plt.xlabel('fractional change in %s'%(param))
                    plt.ylabel('$\\Delta^2$ (mK$^2$)')
                    plt.legend(loc='best')
                    plt.savefig(plotDir+'/fit_pS_%s_k%.2f_z%.2f.png'%(param,kMatrix[nn,mm],z))
                    plt.close()
                
            
            plt.plot(steps,(dGrid[:,nn]-dGrid[indZero,nn])/dGrid[indZero,nn],'o',label='z=%.2f'%(zAxis[nn]))
            plt.xlabel('Fractional Change in %s'%(param))
            plt.ylabel('Fractional Change in $\\Delta^2$')
            plt.legend()
            plt.ylim([-1.,1.])
            plt.savefig(plotDir+'/pSV%s_k%f_z%f.png'%(param,kMatrix[nn,mm],z))
            plt.close()
        
        #plot ps as a function of z
        #also plot dPs as a function of z
        kAxis=np.array(kMatrix[0,:])
        #interpolate derivGrid
        derivGridOut=np.zeros((len(zList),derivGrid.shape[1]))
        for nn in range(len(kAxis)):
            dpVec=derivGrid[:,nn]
            dpsFunc=interp.interp1d(zAxis,dpVec)
            derivGridOut[:,nn]=dpsFunc(zList)

        np.savez(dataDir+'/dDelta_%s.npz'%(param),dGrid=derivGridOut,zAxis=zList,kAxis=kAxis)
        for oo in range(nstep):
            plt.plot(zAxis,dGrid[oo,:])
            plt.text(zAxis[-1],dGrid[oo,-1],'%s=%.2f X fiducial'%(param,1.+steps[oo]),fontsize=8,ha='left')
        plt.xlabel('z')
        plt.ylabel('$\\Delta^2$ (mK)$^2$')
        plt.yscale('log')
        plt.ylim(1e-3,1e4)
        plt.savefig(plotDir+'/pSVZ_%s_k%f.png'%(param,kMatrix[0,mm]))
        plt.close()
        if(nstep>1):
            plt.plot(zAxis,np.abs(derivGrid[:,mm]))
            plt.xlabel('z')
            plt.ylabel('d$\\Delta^2$/d$\\theta$')
            plt.yscale('log')
            plt.savefig(plotDir+'/dPsVz_%s_k%.2f.png'%(param,kMatrix[0,mm]))
            plt.close()
            

    #for each k bin, plot evolution of different models
    



        
        
        
configfile=sys.argv[1]
computeSteps(configfile)
