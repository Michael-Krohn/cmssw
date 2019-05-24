#!/bin/bash

STARTDIR=$(pwd)
OUTDIR=/hdfs/cms/user/revering/dphoton/ZMM_14TeV_sim/
job=$1
scriptname=$2
secondscriptname=$3

mkdir -p $OUTDIR
scratch=$STARTDIR/scratch.$job
mkdir -p $scratch
cd $scratch
cmsRun $STARTDIR/$scriptname $job 
cmsRun $STARTDIR/$secondscriptname $job
rm ZMM14TeV_RAW_"$job".root
mv ZMM14TeV_RECO_"$job".root $OUTDIR
cd $STARTDIR
rmdir $scratch
# -lh /local/cms/user/revering/madgraph/scratch  | grep revering
