#!/bin/bash
#Aaron Ewall-Wice
#May 12th 2015
#script for stepping a single parameter in 21cmFAST
#assumes init.c has already be run
#than steps given parameter by delta, (given fiducial heating PARAM file)
#copies directory/data cubes
#builds programs
#and submits grid engine jobs

config=$1
thisDir=$(pwd)
#initialize config
source ${config}
#change CPU and RAM settings
#cat 21cmFAST/Parameter_files/INIT_PARAMS.H | sed 's/NUMCORES (int) [0-9]\{1,2\}/NUMCORES (int) '${NCPU}'/' > 21cmFAST/Parameter_files/INIT_PARAMS.H
cat 21cmFAST/Parameter_files/INIT_PARAMS_BACKUP.H | sed 's/NUMCORES (int) [0-9]\{1,2\}/NUMCORES (int) '${NCPU}'/' > 21cmFAST/Parameter_files/INIT_PARAMS.H
#cat 21cmFAST/Parameter_files/INIT_PARAMS.H | sed 's/RAM (float) [0-9]\{2\}/NUMCORES (float) '${RAM}'/' > 21cmFAST/Parameter_files/INIT_PARAMS.H
cat 21cmFAST/Parameter_files/INIT_PARAMS.H | sed 's/RAM (float) [0-9]\{2\}/RAM (float) '${RAM}'/' > 21cmFAST/Parameter_files/INIT_PARAMS_NEW.H
mv 21cmFAST/Parameter_files/INIT_PARAMS_NEW.H 21cmFAST/Parameter_files/INIT_PARAMS.H
#set TS_VERBOSE to zero, we need to conserve disk space for this. 
cat 21cmFAST/Parameter_files/HEAT_PARAMS.H | sed 's/Ts_verbose (int) 1/Ts_verbose (int) 0/'
#initialize DM boxes and PS will change this to qsub later
cd 21cmFAST/Programs
#make programs
#make clean
#make
#initialize PS
#./initOnly
cd ${thisDir}
#initialize step param files
python reWriteConfig.py ${config}
#get number of steps from number of config files

#need to ajdust this to allow for multiple param files to be written. Simple Soluton: Make a directory with param files and copy everything in that
ls -d *${STEPPARAM}*_Parameter_files > paramList.txt
nsteps=$(wc -l paramList.txt)
#copy 21cmFAST directory to different step dirs
#make soft links to the DM boxes and logs
mkdir geOut
while read paramDir
do
    step=$(grep -o "[0-9]\{1,3\}"<<<${paramDir})
    echo ${step}
    newDir=21cmFAST_${STEPPARAM}_step${step}
    if ! [ -f "${newDir}" ] 
    then
	mkdir ${newDir}
    fi
    cp -r 21cmFAST/* ${newDir}/
    mv ${paramDir}/* ${newDir}/Parameter_files/
    rm -rf ${paramDir}
    #qsub, will run when grid isn't as full
    qsub -l h_vmem=3G,h_rt=120:00:00 -S /bin/bash -cwd -V -j y -N ${STEPPARAM}_step_${step} -o $(pwd)/geOut/ -pe chost ${NCPU} run_logZscroll.sh ${newDir}
    #./run_logZscroll.sh ${newDir}
done<paramList.txt
