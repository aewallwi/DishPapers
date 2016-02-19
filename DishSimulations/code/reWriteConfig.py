#!/usr/bin/env python
#Aaron Ewall-Wice
#aaronew@mit.edu
#May 13th 2015
#Program to parse a 21cmFAST config file and replace parameters
import numpy as np
import sys
import os.path
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

#replace a parameter in 21cmFAST config with val
def replaceParam(fileName,outName,paramName,val):
    sVal=str(val)
       
    fIn=open(fileName)
    lines=fIn.readlines()
    fIn.close()
    
    for mm,line in enumerate(lines):
        tokens=line.split(' ')
        if(len(tokens)>=4):
            if(tokens[1]==paramName):
                if(tokens[3][-1]=='\n'):
                    tokens[3]=tokens[3][:-1]
                line=line.replace(tokens[3],'('+sVal+'*'+tokens[3]+')')
                lines[mm]=line
    fOut=open(outName,'w')
    for line in lines:
        fOut.write(line)
    fOut.close()

def getParam(fileName,paramName):
    fIn=open(fileName)
    lines=fIn.readlines()
    fIn.close()
    
    output=0.
    for line in lines:
        tokens=line.split(' ')
        if(len(tokens)>4):
            if(tokens[1]==paramName):
                output=float(tokens[3])
    return output


def stepParam(configName):
    #make a fractional step
    config=readConfig(configName)
    inFile='21cmFAST/Parameter_files/%s'%(config['PARAMFILE'])
    paramName=config['STEPPARAM']
    steps=config['STEPFRACS']
    stepNums=config['STEPNUMS']
    #generate steps
    stepVals=[1.+steps[mm] for mm in range(len(steps))]
    for mm in range(len(steps)):
        dirName=os.getcwd()+'/'+paramName+'_step%d_Parameter_files'%(stepNums[mm])
        print dirName
        if(not(os.path.isdir(dirName))):
            os.mkdir(dirName)
        outName=dirName+'/'+inFile.split('/')[-1]
        replaceParam(inFile,outName,paramName,stepVals[mm])

    #if param name is Tvir, than need to also rewrite ANAL PARAM
    if(paramName=='X_RAY_Tvir_MIN'):
        paramNameExtra='ION_Tvir_MIN'
        inFileExtra='21cmFAST/Parameter_files/ANAL_PARAMS.H'
        for mm in range(len(steps)):
            dirName=os.getcwd()+'/'+paramName+'_step%d_Parameter_files'%(stepNums[mm])
            outName=dirName+'/'+inFileExtra.split('/')[-1]
            replaceParam(inFileExtra,outName,paramNameExtra,stepVals[mm])



stepParam(sys.argv[1])
