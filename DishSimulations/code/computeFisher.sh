#!/bin/bash
#compute fisher matrix from w vectors of various parameters
#permit performing cuts in k and z and whether one wants to perform independent constraints 

#config file has a list of parameters, a delta in frequency, a minimum redshift and a maximum redshift. Derive both independent constraints with each redshift interval or derive constraints from sums over redshift intervals. 
config=$1

#make directory to store outputs
source ${config}
echo ${DIRNAME}
mkdir ./${DIRNAME}
python computeFisher.py ${config}
