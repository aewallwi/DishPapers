#!/bin/bash
echo "running on ${HOSTNAME}"
thisdir=$(pwd)
folder=$1
ssh ${HOSTNAME} "cd ${thisdir}; cd ${folder}/Programs; make clean; make; ./drive_logZscroll_Ts"
mkdir ${thisdir}/${folder}/dTb
mv ${thisdir}/${folder}/Boxes/delta_T* ${thisdir}/${folder}/dTb/
#remove all boxes except for brightness temperature output
rm -rf ${thisdir}/${folder}/Boxes/*
