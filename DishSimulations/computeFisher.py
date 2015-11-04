#!/usr/bin/env python
#compute fisher matrix for optimistic, pessimistic and fiducial foregrounds scenarios using both joint redshift constraints and for each redshift independently.Invert the fisher matrix to get the covariance matrix. Plot error ellipses and write out covariance matrices to text files for each
#redshift
#and joint over all redshifts. 

import scipy.misc as misc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys


labelIndex={'ZETA_X':'$\\frac{\\Delta \\zeta_X}{\\zeta_X}$','X_RAY_Tvir_MIN':'$\\frac{\\Delta T_{vir}}{T_{vir}}$','X_RAY_SPEC_INDEX':'$\\frac{\\Delta \\alpha}{\\alpha}$','NU_X_THRESH':'$\\frac{\\Delta \\nu_{min}}{\\nu_{min}}$'}
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
                config[params[0]] = eval(params[1])
    print config
    return config

#compute the fisher matrix for a list of wMatrices
def computeCov(wList):
    nParams=len(wList)
    fisherMat=np.zeros((nParams,nParams))
#    print 'wList=%s'%wList
    for ii in range(nParams):
        wList[ii][np.logical_or(np.isnan(wList[ii]),np.isinf(wList[ii]))]=0.#all nans set to zero. (where do they come from?)
        for jj in range(ii+1):
#            print 'product = %.2e'%(np.sum((wList[ii]*wList[jj]).flatten()))
            fisherMat[ii,jj]=np.sum((wList[ii]*wList[jj]).flatten())
            fisherMat[jj,ii]=fisherMat[ii,jj]
    #invert fisherMat to get covariance matrix
#    print fisherMat
    try:
        covMat=np.linalg.inv(fisherMat)
    except:
        print 'singular'
        covMat=None
    print covMat
    return covMat

def f2z(z):
    return f/1420.-1.


#compute the independent covariance matrices 
#scen is a string describing whether we are talking about the optimistic, moderat, or pessimistic foreground model. 

def computeCovariances(paramlist,scen,zmin,zmax,df,noFM=False):
    #for each parameter in list, load wMatrix, select redshifts that are df MHz apart between zmin and zmax, filter down wMatrix 
    wList=[]
    zListAll=[]
    for param in paramlist:
        wMat=[]
        zList=[]
        print 'param=%s'%param
        wFiles=np.load(param+'_Output/data/wMatrix_%s_%s.npz'%(param,scen))
        zAxis=wFiles['zAxis']
#        print zAxis
        zminL=zAxis[np.where(np.abs(zAxis-zmin)==np.abs(zAxis-zmin).min())]
        zmaxL=zAxis[np.where(np.abs(zAxis-zmax)==np.abs(zAxis-zmax).min())]
        z=zminL
        while(z<=np.min([zmaxL,zAxis.max()])):
            zInd=np.where(np.abs(zAxis-z)==np.abs(zAxis-z).min())[0]
            if(not(noFM) or (z>=1420./88.-1. and z<1420./108.-1.)):
                wMat.append(wFiles['wMat'][zInd,:])
                zList.append(zAxis[zInd])
            elif(z>=1420./80.-1. or z<1420./80.-1.):
                wMat.append(wFiles['wMat'][zInd,:])
                zList.append(zAxis[zInd])
            dz=df/(1420.)*(1+z)**2.
            z+=dz
        print 'wMat=%s'%wMat
        wMat=np.array(wMat)
        zList=np.array(zList)
#        print zList
        zListAll.append(zList)
        wList.append(wMat)
    #now compute the covariance matrices for each wMat
    zListAll=np.array(zListAll).squeeze()
    print zListAll
    print zListAll.shape
    zTest=True
    for mm in range(zListAll.shape[0]):
        if(not(np.all(zListAll[mm,:]==zListAll[0,:]))):
            print 'warning z does not match!'
    #create full covMat from all redshift
    covMatFull=computeCov(wList)
    
    #create list of covMatrices for each redshift
    covMatList=[]
 
    for zInd in range(len(zListAll[0])):
        wListT=[]
        for wMat in wList:
            wRow=wMat[zInd,:]
            wListT.append(wRow)

        covMatList.append(computeCov(wListT))

    #return z axis, a joint covariance matrix, and a list of covariance matrices. 
    return zListAll[0,:],covMatFull,covMatList



#basic ellipse plotting. Will likely make production plots in an ipynotebook
def drawEllipses(paramNames,covList,plotName):
    #create 2 sigma error ellipses from a list of covariance matrices
     nParams=len(paramNames)
     nPlots=int(np.round(misc.comb(nParams,2)))
     print 'nPlots=%d'%nPlots
     fig,axarr=plt.subplots(nParams-1,nParams-1)
     print len(axarr)
     for ii in range(nParams):
         for jj in range(ii+1,nParams):
             ax=axarr[ii,jj-1]
             for cov in covList:
                 print cov
                 a=np.sqrt((cov[ii,ii]+cov[jj,jj])/2.+np.sqrt((cov[ii,ii]-cov[jj,jj])**2./4.+cov[ii,jj]**2.))
                 b=np.sqrt((cov[ii,ii]+cov[jj,jj])/2.-np.sqrt((cov[ii,ii]-cov[jj,jj])**2./4.+cov[ii,jj]**2.))
                 theta=.5*np.arctan(2.*cov[ii,jj]/(cov[ii,ii]-cov[jj,jj]))
                 print 'a=%.2f,b=%.2f,theta=%.2f'%(a,b,theta)
                 ell=matplotlib.patches.Ellipse(xy=[0.,0.],width=a,height=b,angle=theta*180./np.pi)
                 ax.set_xlim(-.2,.2)
                 ax.set_ylim(-.2,.2)
                 ax.add_artist(ell)
                 ell.set_clip_box(ax.bbox)
                 ell.set_facecolor('none')
                 ell.set_edgecolor('k')
             ax.set_xlabel(labelIndex[paramNames[ii]])
             ax.set_ylabel(labelIndex[paramNames[jj]])
     fig.set_size_inches([(nParams-1)*5,(nParams-1)*5.])
     plt.savefig(plotName+'.png')

#first read config
config=sys.argv[1]
config=readConfig(config)
#get list of parameters
paramList=open(config['PARAMS']).readlines()
for mm in range(len(paramList)):
    paramList[mm]=paramList[mm][:-1]#remove new lines
ZMIN=config['ZMIN']
ZMAX=config['ZMAX']
DNU=config['DELTAF']
NAME=config['DIRNAME']
covListAll=[]
for scen in ['opt','mod','pess']:
    #compute the Fisher matrix
    zList,covJoint,covList=computeCovariances(paramList,scen,ZMIN,ZMAX,DNU)
    covListAll.append(covJoint)
    #plot it
    drawEllipses(paramList,covList,NAME+'/ConfidenceEllipses_'+scen+'_zmin%.2f_zmax%.2f'%(ZMIN,ZMAX))
    drawEllipses(paramList,[covJoint],NAME+'/ConfidenceEllipses_Joint_'+scen+'_zmin%.2f_zmax%.2f'%(ZMIN,ZMAX))
    #save covariances
    np.savez(NAME+'/CovarianceJoint_'+scen+'_zmin%.2f_zmax%.2f.npy'%(ZMIN,ZMAX),cov=np.array(covList),zList=zList)
    np.savez(NAME+'/Covariances_'+scen+'_zmin%.2f_zmax%.2f.npz'%(ZMIN,ZMAX),cov=covJoint,zList=zList)
#now draw joint ellipses for opt,mod,pess scenarios
drawEllipses(paramList,covListAll,NAME+'/ConfidenceEllipses_Joint_Comparison_zmin%.2f_zmax%.2f'%(ZMIN,ZMAX))
    
    

