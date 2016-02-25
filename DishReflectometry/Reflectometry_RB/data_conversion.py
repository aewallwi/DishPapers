import numpy as np, pylab as plt, aipy as a
import sys, csv

def fromcsv1(filename):
    print 'Reading', filename
    d = csv.reader(open(filename,'r'), delimiter=',')
    x = np.array(list(d)[18:-1], dtype=np.float)
    return x[:,0], x[:,1]
    
    fq_rb,db_rb = fromcsv1('measurement_dishAndFeed_DB.csv')
    fq_rb,ph_rb = fromcsv1('measurement_dishAndFeed_P.csv')
    
    fq_rb = fq_rb*1000000
    
    with open('measurement_dishAndFeed1_DB.csv', 'wb') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow(['Spam'] * 5 + ['Baked Beans'])
    spamwriter.writerow(['Spam', 'Lovely Spam', 'Wonderful Spam'])