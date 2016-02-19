/*************************************************************
Aaron Ewall-Wice
aaronew@mit.edu
October 18th 2014
c routines for interferometer. written at performance 
sensitive steps to make it go faster
************************************************************/
#include <stdio.h>
#include <math.h>
#define PI (3.141592653589793)
typedef struct {double real; double imag;} cdouble;
typedef struct {float real; float imag;} cfloat;
/*************************************************************
grid data to a 2D grid stored as a vector 
using nearest neighbor gridding
dx,dy: grid spacing
cxout,cyout: center values 
************************************************************/
void griddata_nn(long * xi, long * yi, cdouble * input, cdouble * output,long ndata, long nx, long ny)
{
  int m;
  int x,y;
  for(m=0;m<ndata;m++)
    {
      x=(int)xi[m];
      y=(int)yi[m];
      if(y<ny && x<nx && x>=0 && y>=0)
	{
	  output[y*nx+x].real=output[y*nx+x].real+input[m].real;
	  output[y*nx+x].imag=output[y*nx+x].imag+input[m].imag;
	}
     }
}
/***********************************************************
get the number uniformly gridded 
************************************************************/
void count_nn(long * xi, long * yi, long * counts, long ndata, long nx, long ny)
{
  long m,x,y;
  for(m=0;m<ndata;m++)
    {
      x=xi[m];
      y=yi[m];
      if(y<ny && x<nx && y>=0 && x>=0)
	{
	  counts[y*ny+x]+=1;
	}
    }
}

/************************************************************
//take differences in values of 1d array
//n is the length of the input array
//output array is assumed to be n*(n-1)/2 as long 
//as output
//only write function for doubles. must cast 
//other data types in python
************************************************************/
void diffcomb(double * x, double * out, long n)
{
  int l,m;
  for (l=0;l<n;l++)
    {
      for (m=0;m<l;m++)
	{
	  *out=x[l]-x[m];
	  out++;
	}
    }
}
/************************************************************
compute the visibilities for a list of uvw coordinates
from a list of l,m coords with flux s
 ************************************************************/
void computevisibilities(double * u, double * v, double * w, double * l, double * m, double * s,cdouble * out,long nvis,long nsrc)
{    
  double lt,mt,st,ut,vt,wt,vistr,visti,exparg,n;
  int i,j;
  for(i=0;i<nvis;i++)
    {
	out->real=0.;
	out->imag=0.;
	ut=*u;
	vt=*v;
	wt=*w;
	for(j=0;j<nsrc;j++)
	  {
	    lt=l[j];
	    mt=m[j];
	    st=s[j];
	    n=sqrt(1.-lt*lt-mt*mt);
	    exparg=-2.*PI*(lt*ut+mt*vt+(n-1.)*wt);
	    vistr=st*cos(exparg);
	    visti=st*sin(exparg);
	    out->real=out->real+vistr;
	    out->imag=out->imag+visti;
	  }
	u++;
	v++;
	w++;
	out++;
    }
}  


