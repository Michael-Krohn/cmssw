#!/bin/bash

STARTDIR=$(pwd)
OUTDIR=$STARTDIR/ZTT_13TeV/
job=$1
scriptname=$2
job=$((job+400))
mkdir -p $OUTDIR
scratch=$STARTDIR/scratch.$job
mkdir -p $scratch
cd $scratch
source /local/grid/cmssoft/cms/cmsset_default.sh
cmsenv
cmsRun $STARTDIR/$scriptname $job 
cmsRun $STARTDIR/ZTT13TeV_pythia8_SIMtoL1_cfg.py $job
cmsRun $STARTDIR/ZTT13TeV_pythia8_L1toAOD_cfg.py $job
rm $scratch/ZTT13TeV_pythia8"$job"_GENSIM.root
rm $scratch/ZTT13TeV_pythia8"$job"_step1.root
mv $scratch/ZTT13TeV_pythia8_"$job"_AOD.root $OUTDIR
cd $STARTDIR
rm -rf $scratch
# -lh /local/cms/user/revering/madgraph/scratch  | grep revering
