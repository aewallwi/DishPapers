import numpy as np, pylab as plt, aipy as a
import sys, csv

def fromcsv1(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = np.array(list(d)[18:-1], dtype=np.float)
    return x[:,0]/1e9, x[:,1]
    
    
if True:  
    
    fq_feed,db_feed = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_FEED_DB.csv') # Reading the magnitude and phase of the feed only datafile to calibrate for the zero point/ 
    fq_feed,ph_feed = fromcsv1('/Users/Punthipatra/WORKAREA1/DATA_Area/Project_HERA/DishPapers/DishReflectometry/GB_reflectometry_part2/HERA_FEED_P.csv')
    print fq_feed