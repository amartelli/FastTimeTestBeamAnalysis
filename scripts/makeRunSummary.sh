#!/bin/bash

outDir="~/www/HGCal/FastTimeTB/MIPCalibrationRuns"
siChargeEst=("wave_fit_smallw_ampl" "charge_integ_smallw_mcp")
siTriggerThr=(2 3 4)
fitMode=(0 1 2)
for fit in ${fitMode[@]}; do 
    for est in ${siChargeEst[@]}; do
	for thr in ${siTriggerThr[@]}; do
	    echo "${fit} ${thr} ${outDir} ${est}"
	    root -l -b -q ../macros/makeRunSummary.C+\(${fit},${thr},\"${outDir}\",\"${est}\"\)
	done
    done
done

a=(`ls ~/www/HGCal/FastTimeTB`)
for i in ${a[@]}; do
    echo $i
    python drawRunSummary.py ~/www/HGCal/FastTimeTB/${i}/summary.root;
done


