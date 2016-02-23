#************************************************************
#supplementary c functions for interferometer class
#aaronew@mit.edu
#makes things go fase
#************************************************************
#!/usr/bin/env python
#************************************************************
#load ctype functions
#************************************************************
import os,sys
import numpy as np

_path=os.path.dirname('__file__')
print _path
if(sys.platform=='darwin'):
    libstr='interferometer_darwin'
elif(sys.platform=='linux2'):
    libstr='interferometer_linux'
lib=np.ctypeslib.load_library(libstr,_path)
dcomb=getattr(lib,'diffcomb')
cvis=getattr(lib,'computevisibilities')
griddata_nn=getattr(lib,'griddata_nn')
count_nn=getattr(lib,'count_nn')
#set ctype function attributes
dcomb.restype=None
dcomb.argtypes=[np.ctypeslib.ndpointer(float,flags='aligned, contiguous'),np.ctypeslib.ndpointer(float,flags='aligned, contiguous'),np.ctypeslib.c_intp]
cvis.restype=None
cvis.argtypes=[np.ctypeslib.ndpointer(float,flags='aligned, contiguous'),
               np.ctypeslib.ndpointer(float,flags='aligned, contiguous'),
               np.ctypeslib.ndpointer(float,flags='aligned, contiguous'),
               np.ctypeslib.ndpointer(float,flags='aligned, contiguous'),
               np.ctypeslib.ndpointer(float,flags='aligned, contiguous'),
               np.ctypeslib.ndpointer(float,flags='aligned, contiguous'),
               np.ctypeslib.ndpointer(complex,flags='aligned, contiguous'),
               np.ctypeslib.c_intp,
               np.ctypeslib.c_intp]

griddata_nn.restype=None
griddata_nn.argtypes=[np.ctypeslib.ndpointer(np.ctypeslib.c_intp,flags='aligned,contiguous'),
                      np.ctypeslib.ndpointer(np.ctypeslib.c_intp,flags='aligned,contiguous'),
                      np.ctypeslib.ndpointer(complex,flags='aligned,contiguous'),
                      np.ctypeslib.ndpointer(complex,flags='aligned,contiguous'),
                      np.ctypeslib.c_intp,
                      np.ctypeslib.c_intp,
                      np.ctypeslib.c_intp]
count_nn.restype=None
count_nn.argtypes=[np.ctypeslib.ndpointer(np.ctypeslib.c_intp,flags='aligned,contiguous'),
                   np.ctypeslib.ndpointer(np.ctypeslib.c_intp,flags='aligned,contiguous'),
                   np.ctypeslib.ndpointer(int,flags='aligned,contiguous'),
                   np.ctypeslib.c_intp,
                   np.ctypeslib.c_intp,
                   np.ctypeslib.c_intp]
                   

#wrapper functions functions contained in interferometer.so
def diffcomb(x):
    #convert x to float
    requires=['CONTIGUOUS','ALIGNED']
    x=np.require(x,float,requires)
    y=np.empty(x.size*(x.size-1)/2,dtype=float)
    y=np.require(y,float,requires)
    lib.diffcomb(x,y,x.size)
    return y

#compute the visibilites for sources with 
#flux s (1xn array)
#located at locations lm (nx2 numpy array)
#at uvw given by uvw (nx3 numpy array)
def computevisibilities(uvw,lm,s):
    #only compute redundent visibilities
    #round visibilities to millions of a wavelength, otherwise bit errors will falsely miss redundant baselines
    #filter nans, zeros,and lm**2>1.
    nanzero=np.invert(np.logical_or(np.logical_or(np.isinf(s),np.logical_or(np.isnan(s),s==0.)),np.sum(lm**2.,axis=1)>=1.))
    s=s[nanzero]
    lm=lm[nanzero,:]
    
    
    uvw=np.around(uvw,6)
    requires=['CONTIGUOUS','ALIGNED']
    u=np.require(uvw[:,0],float,requires)
    v=np.require(uvw[:,1],float,requires)
    w=np.require(uvw[:,2],float,requires)
    l=np.require(lm[:,0],float,requires)
    m=np.require(lm[:,1],float,requires)
    s=np.require(s,float,requires)
    output=np.empty(uvw.shape[0],dtype=complex)
    output=np.require(output,complex,requires)
    lib.computevisibilities(u,v,w,l,m,s,output,uvw.shape[0],lm.shape[0])
    return output


def count_nn_2D(xy,dx,dy,nx,ny,cx=0.,cy=0.):
    requires=['CONTIGUOUS','ALIGNED']
    x=(np.round((xy[:,0]-cx)/dx)+nx/2).astype(int)
    y=(np.round((xy[:,1]-cy)/dy)+ny/2).astype(int)
    x=np.require(xy[:,0],int,requires)
    y=np.require(xy[:,1],int,requires)
    counts=np.require(np.zeros(nx*ny),int,requires)
    lib.count_nn(x,y,counts,x.shape[0],int(nx),int(ny))
    return counts


def griddata_nn_2D(xy,data,dx,dy,nx,ny,cx=0.,cy=0.):
    requires=['CONTIGUOUS','ALIGNED']
    output=np.require(np.zeros(nx*ny),complex,requires)
    x=(np.round((xy[:,0]-cx)/dx)+nx/2.)
    y=(np.round((xy[:,1]-cy)/dy)+ny/2.)
    x=np.require(x,int,requires)
    y=np.require(y,int,requires)
    data=np.require(data,complex,requires)
    nx=int(nx)
    ny=int(ny)
    lib.griddata_nn(x,y,data,output,x.shape[0],nx,ny)
    output=output.reshape((ny,nx))
    return output
