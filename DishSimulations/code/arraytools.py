#!/usr/bin/env python
import numpy as np
#************************************************************
#Aaron Ewall-Wice
#aaronew@mit.edu
#October 20th 2014
#methods to sort 2-D numpy arrays by column or row.
#************************************************************
#sort a 2d array of floats by a specified column
#************************************************************
def columnsort(a,column=0,order='desc'):
    print a.shape
    asorted=np.zeros(a.shape)
    dtype=[('x',float),('y',int)]
    sortkey=np.zeros(a.shape[0],dtype=dtype)
    sortkey['x']=a[:,column]
    indices=range(a.shape[0])
    sortkey['y']=indices
    sortkey.sort(order='x')
    for col in range(a.shape[1]):
        asorted[:,col]=a[sortkey['y'],col]
    if(order=='desc'):
        asorted=np.flipud(asorted)
    return asorted
        
        
#************************************************************
#sort a 2d array of floats by a specified column
#************************************************************
def rowsort(a,row=0,order='desc'):
    return columnsort(a.T,row,order).T


#************************************************************
#get unique rows of an array, courtesy bfroehle
#https://github.com/numpy/numpy/issues/2871
#************************************************************
def unique_rows(A, return_index=False, return_inverse=False):
    """
    Similar to MATLAB's unique(A, 'rows'), this returns B, I, J
    where B is the unique rows of A and I and J satisfy
    A = B[J,:] and B = A[I,:]

    Returns I if return_index is True
    Returns J if return_inverse is True
    """
    A = np.require(A, requirements='C')
    assert A.ndim == 2, "array must be 2-dim'l"

    B = np.unique(A.view([('', A.dtype)]*A.shape[1]),
               return_index=return_index,
               return_inverse=return_inverse)

    if return_index or return_inverse:
        return (B[0].view(A.dtype).reshape((-1, A.shape[1]), order='C'),) \
            + B[1:]
    else:
        return B.view(A.dtype).reshape((-1, A.shape[1]), order='C')

#************************************************************
#get unique columns
#************************************************************
#def unique_columns(A,return_index=False,return_inverse=False):
#    if(return_index or return_inverse):
#        
#    unique_rows(A.T,return_index=return_index,return_inverse=return_inverse)
