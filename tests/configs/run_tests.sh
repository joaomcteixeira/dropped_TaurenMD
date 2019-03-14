#! /bin/bash

rm _*
rm *.csv
rm plot*
rm taurenmd.log
rm traj_OUTPUT.dcd

source activate taurenmd

CONFIGS=(
    # basic test with RMSDS using MDTraj
    test_1.json
    ## basic test with RMSDs using MDAnalysis
    test_2.json
    ## basic test with atom_selection selecting chain 1 - MDTraj
    test_3.json
    ## basic test with atom_selection selecting chain 1 - MDTraj
    test_4.json
    ## savetraj with atom and frame selection - MDTraj
    test_5.json
    ## savetraj with atom and frame selection - MDAnalysis
    test_6.json
    ## combined RMSDs with atom selection (chain) - MDTraj
    test_7.json
    ## same as above with empty selection - MDTraj
    test_8.json
    ## combined RMSDs with atom selection (chain) - MDAnalysis
    test_9.json
    ## same as above with empty selection - MDAnalysis
    test_10.json
    # combined RMSDs with atom selection (resid) - MDAnalysis
    test_11.json
    )


for i in ${CONFIGS[@]}
do
    echo "************************************${i}"
    python $1 -c ${i}
done
