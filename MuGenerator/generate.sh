#!/bin/bash

STARTDIR=$(pwd)
OUTDIR=/hdfs/cms/user/revering/dphoton/ZMM_14TeV_sim
job=$1
scriptname=$2

mkdir -p $OUTDIR
scratch=$STARTDIR/scratch.$job
mkdir -p $scratch
cd $scratch
cmsRun $STARTDIR/$scriptname $job 
mv ZMM_14TeV_file"$job".root $OUTDIR
rm -rf $scratch
# -lh /local/cms/user/revering/madgraph/scratch  | grep revering
