#!/bin/bash
declare -a SignalPoints=("1p0" "0p4" "0p2")
#declare -a SignalPoints=("0p2")

source setup.sh

echo "Preparing histograms"
python prepareHistograms.py -i dataFiles/darkPhotonCombineInput_23_12_20_InvariantMass.root
#python prepareHistograms.py -i dataFiles/darkPhotonCombineInput_23_10_9_InvariantMass.root

echo "Extract histograms"
python extractHistograms.py -i Modified/dataFiles/darkPhotonCombineInput_23_12_20_InvariantMass.root
#python extractHistograms_BDT_invariantMass_nProbes.py -i Modified/dataFiles/darkPhotonCombineInput_23_10_9_InvariantMass.root

echo "Make datacards"
python makeDatacards.py

echo "Run limits - fitDiagnostic"
for signal in ${SignalPoints[@]};
  do
    echo $signal
    combine -M FitDiagnostics CombineInput/Combined_DBrem_$signal".txt" --saveNormalizations  --saveShapes --saveWithUncertainties --saveOverallShapes --keepFailures -n _DBrem_$signal -t -1
    mv fitDiagnostics_DBrem_$signal".root" rootfiles/Run2Legacy_v4__Default/FitDiagnostics/
  done

echo "Extracting pre- and post-fit histograms"
python python/MakePrePostFitHistogram.py

echo "Draw pre-fit histograms"
python python/Draw_Prefit_1Canvas.py --blind

echo "Draw post-fit histograms"
python python/Draw_Postfit_1Canvas.py --blind
