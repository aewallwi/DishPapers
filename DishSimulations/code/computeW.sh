#!/bin/bash
#driver for computing the fisher matrix for a list of parameters. 

#start with 1 param for now
config=$1
source ${config}
thisdir=$(pwd)
if ! [ -d "${thisdir}/${STEPPARAM}_${LABEL}_Output" ]
then
    mkdir ${thisdir}/${STEPPARAM}_${LABEL}_Output 
    python computeDerivative.py ${config}
fi
#first determine the redshifts to compute sensitivity and power spectra at
#than for each model, match redshifts and calculate sensitivity. 
python redshiftFreqs.py ${ZMIN} ${ZMAX} ${DZ} ${thisdir}/${STEPPARAM}_${LABEL}_Output


#now run 21cmSense for each reshift
#get number of steps by listing all run directories
ls -d 21cmFAST_${STEPPARAM}_step*/>${STEPPARAM}_steplist.txt
NSTEP=$(wc -l ${STEPPARAM}_steplist.txt | grep -o "[0-9]\{1,2\}")
cd ${STEPPARAM}_${LABEL}_Output/data
#now need to compute sensitivity for all steps in ps
 #generate array file for array
#if ! [ -f "${ARRAY}.track_${TRACK}hr_arrayfile.npz" ]
if [ "${ARRAY}" = "mwa" ]  || [ "${ARRAY}" = "lofar" ] 
then
    if ! [ -f "${ARRAY}.track_6.0hr_arrayfile.npz" ]
    then
	if [ "${ARRAY}" = "mwa" ] 
	then
	    echo "using 200 max bl"
	    python ${thisdir}/21cmSense-master/mk_array_file.py -C ${ARRAY} --bl_max 200 --track 6
	else
	    python ${thisdir}/21cmSense-master/mk_array_file.py -C ${ARRAY} --track 6
	fi

    fi
else
    if ! [ -f "${ARRAY}.drift_arrayfile.npz" ]
    then

	python ${thisdir}/21cmSense-master/mk_array_file.py -C ${ARRAY} #--track ${TRACK} #--bl_max 200.
    fi
fi

for stepnum in `seq 0 $[${NSTEP}-1]`
do
echo "nstep=${NSTEP}"
echo "stepnum=${stepnum}"
ls ps_${STEPPARAM}_step${stepnum}_*>psList.txt
while read psName
do
    z=$(grep -o "[0-9]\{1,3\}.[0-9]\{1,2\}"<<<$(grep -o "z[0-9]\{1,3\}.[0-9]\{1,2\}"<<<${psName}))
    echo "z=${z}"
    freq=$(python ../../z2f.py ${z})
    echo "f=${freq}"
    #get nchan and bandwidth 
    bwidth=$(python ../../z2df.py ${z} ${DZ})
    nchan=$(python ../../nChan.py ${bwidth} ${DF})
    echo "bwidth=${bwidth}"
    echo "nchan=${nchan}"
    if [ "${SAMPLEV}" = "False" ]
    then
	eorFlag="--eor ps_zeros.txt"
    else
	eorFlag="--eor ${psName}"
    fi

    mink=$(python ../../mink.py ${z} ${MINDELAY})
    echo "mink=${mink}"
    #generate sensitivity for pess,mod, and opt scenarios and 1 year of observing. 
    if [ "${ARRAY}" = "mwa" ]  || [ "${ARRAY}" = "lofar" ] 
    then
	if ! [ -f "step_${stepnum}_${ARRAY}.track_6.0hr_mod_${freq}.npz" ] 
#    if ! [ -f "step_${stepnum}_${ARRAY}.track_${TRACK}hr_mod_${freq}.npz" ] 
	then
#	python ${thisdir}/21cmSense-master/calc_sense.py -m mod -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth} ${ARRAY}.track_${TRACK}hr_arrayfile.npz
	    python ${thisdir}/21cmSense-master/calc_sense.py -m mod -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth} ${ARRAY}.track_6.0hr_arrayfile.npz
	fi
#    if ! [ -f "step_${stepnum}_${ARRAY}.track_${TRACK}hr_pess_${freq}.npz" ]
	if ! [ -f "step_${stepnum}_${ARRAY}.track_6.0hr_pess_${freq}.npz" ]
	then
#	python ${thisdir}/21cmSense-master/calc_sense.py -m pess -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth}  ${ARRAY}.track_${TRACK}hr_arrayfile.npz
	    python ${thisdir}/21cmSense-master/calc_sense.py -m pess -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth}  ${ARRAY}.track_6.0hr_arrayfile.npz
	fi
#    if ! [ -f "step_${stepnum}_${ARRAY}.track_${TRACK}hr_opt_${freq}.npz" ]
	if ! [ -f "step_${stepnum}_${ARRAY}.track_6.0hr_opt_${freq}.npz" ]
	then
#	python ${thisdir}/21cmSense-master/calc_sense.py -m opt -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth} ${ARRAY}.track_${TRACK}hr_arrayfile.npz
	    python ${thisdir}/21cmSense-master/calc_sense.py -m opt -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth} ${ARRAY}.track_6.0hr_arrayfile.npz
	fi
	
    else
	if ! [ -f "step_${stepnum}_${ARRAY}.drift_mod_${freq}.npz" ] 
#    if ! [ -f "step_${stepnum}_${ARRAY}.track_${TRACK}hr_mod_${freq}.npz" ] 
	then
#	python ${thisdir}/21cmSense-master/calc_sense.py -m mod -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth} ${ARRAY}.track_${TRACK}hr_arrayfile.npz
	python ${thisdir}/21cmSense-master/calc_sense.py -m mod -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth} ${ARRAY}.drift_arrayfile.npz
	fi
#    if ! [ -f "step_${stepnum}_${ARRAY}.track_${TRACK}hr_pess_${freq}.npz" ]
	if ! [ -f "step_${stepnum}_${ARRAY}.drift_pess_${freq}.npz" ]
	then
#	python ${thisdir}/21cmSense-master/calc_sense.py -m pess -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth}  ${ARRAY}.track_${TRACK}hr_arrayfile.npz
	    python ${thisdir}/21cmSense-master/calc_sense.py -m pess -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth}  ${ARRAY}.drift_arrayfile.npz
	fi
#    if ! [ -f "step_${stepnum}_${ARRAY}.track_${TRACK}hr_opt_${freq}.npz" ]
	if ! [ -f "step_${stepnum}_${ARRAY}.drift_opt_${freq}.npz" ]
	then
#	python ${thisdir}/21cmSense-master/calc_sense.py -m opt -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth} ${ARRAY}.track_${TRACK}hr_arrayfile.npz
	    python ${thisdir}/21cmSense-master/calc_sense.py -m opt -b ${mink} -f ${freq} ${eorFlag} --ndays 180. --n_per_day 6. --nchan ${nchan} --bwidth ${bwidth} ${ARRAY}.drift_arrayfile.npz
	fi
    fi

    done<psList.txt
if [ "${ARRAY}" = "mwa" ]  || [ "${ARRAY}" = "lofar" ] 
then
#ok we've processed all redshifts for given step, now move them
    ls ${ARRAY}.track_6.0hr_mod_*.npz>senseListMod.txt
    ls ${ARRAY}.track_6.0hr_opt_*.npz>senseListOpt.txt
    ls ${ARRAY}.track_6.0hr_pess_*.npz>senseListPes.txt
else
    ls ${ARRAY}.drift_mod_*.npz>senseListMod.txt
    ls ${ARRAY}.drift_opt_*.npz>senseListOpt.txt
    ls ${ARRAY}.drift_pess_*.npz>senseListPes.txt
fi
while read senseName
do
    #get frequency
    freq=$(grep -o "[0-9]\{1,2\}.[0-9]\{5\}"<<<${senseName})
    if [ "${ARRAY}" = "mwa" ]  || [ "${ARRAY}" = "lofar" ] 
    then
	mv ${ARRAY}.track_6.0hr_mod_${freq}.npz step_${stepnum}_${ARRAY}.track_6.0hr_mod_${freq}.npz
    else
	mv ${ARRAY}.drift_mod_${freq}.npz step_${stepnum}_${ARRAY}.drift_mod_${freq}.npz
    fi
    done<senseListMod.txt
while read senseName
do
    #get frequency
    freq=$(grep -o "[0-9]\{1,2\}.[0-9]\{5\}"<<<${senseName})
    if [ "${ARRAY}" = "mwa" ]  || [ "${ARRAY}" = "lofar" ] 
    then
	mv ${ARRAY}.track_6.0hr_opt_${freq}.npz step_${stepnum}_${ARRAY}.track_6.0hr_opt_${freq}.npz
    else
	mv ${ARRAY}.drift_opt_${freq}.npz step_${stepnum}_${ARRAY}.drift_opt_${freq}.npz
    fi
    done<senseListOpt.txt

while read senseName
do
    #get frequency
    freq=$(grep -o "[0-9]\{1,2\}.[0-9]\{5\}"<<<${senseName})
    if [ "${ARRAY}" = "mwa" ]  || [ "${ARRAY}" = "lofar" ] 
    then
	mv ${ARRAY}.track_6.0hr_pess_${freq}.npz step_${stepnum}_${ARRAY}.track_6.0hr_pess_${freq}.npz
    else
	mv ${ARRAY}.drift_pess_${freq}.npz step_${stepnum}_${ARRAY}.drift_pess_${freq}.npz
    fi
    done<senseListPes.txt
done

#compute wi as a matrix on which cuts in k and z can be made later. 
if [ "${ARRAY}" = "mwa" ]  || [ "${ARRAY}" = "lofar" ] 
then
    ls step*${ARRAY}.track_6.0hr_mod_*.npz>senseListMod.txt
    ls step*${ARRAY}.track_6.0hr_opt_*.npz>senseListOpt.txt
    ls step*${ARRAY}.track_6.0hr_pess_*.npz>senseListPes.txt
else
    ls step*${ARRAY}.drift_mod_*.npz>senseListMod.txt
    ls step*${ARRAY}.drift_opt_*.npz>senseListOpt.txt
    ls step*${ARRAY}.drift_pess_*.npz>senseListPes.txt
fi
cd ../../

python ${thisdir}/computeW.py ${config} mod
python ${thisdir}/computeW.py ${config} opt
python ${thisdir}/computeW.py ${config} pess