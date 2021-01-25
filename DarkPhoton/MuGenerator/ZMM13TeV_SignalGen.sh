#!/bin/bash

STARTDIR=$(pwd)
OUTDIR=$STARTDIR/ZMM_Signal/
job=$1
scriptname=$2

mkdir -p $OUTDIR
scratch=$STARTDIR/scratch.$job
mkdir -p $scratch
cd $scratch
source /local/grid/cmssoft/cms/cmsset_default.sh
cmsenv
cmsRun $STARTDIR/$scriptname $job 
cmsRun $STARTDIR/ZMM13TeV_pythia8_SIMtoL1_cfg.py $job
#cmsRun $STARTDIR/ZMM13TeV_pythia8_L1toRECO_cfg.py $job
cmsRun $STARTDIR/ZMM13TeV_pythia8_L1toAOD_cfg.py $job
rm $scratch/ZMM13TeV+signal_"$job"_gensim.root
rm $scratch/ZMM13TeV+signal_"$job"_step1.root
mv $scratch/ZMM13TeV+signal_"$job"_AOD.root $OUTDIR
cd $STARTDIR
rm -rf $scratch
# -lh /local/cms/user/revering/madgraph/scratch  | grep revering
