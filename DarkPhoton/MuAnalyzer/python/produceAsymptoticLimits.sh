#!/bin/bash
declare -a SignalPoints=("0p1" "0p05" "1p0" "0p8" "0p2" "0p01" "0p6" "0p4" "0p025")

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
    combine -M AsymptoticLimits CombineInput/Combined_DBrem_$signal".txt" --run blind > LogsAsymptoticLimits/Log_Combined_DBrem_$signal".txt"

    combine -M AsymptoticLimits CombineInput/completeDisappearance_DBrem_$signal".txt" --run blind > LogsAsymptoticLimits/Log_CompleteDisappearance_DBrem_$signal".txt"

    combine -M AsymptoticLimits CombineInput/partialDisappearanceBDT_DBrem_$signal".txt" --run blind > LogsAsymptoticLimits/Log_PartialDisappearanceBDT_DBrem_$signal".txt" 
  done

echo "Extracting pre- and post-fit histograms"
python python/read_Asymptotic.py -i LogsAsymptoticLimits

echo "Plotting limits"
python python/Draw_Limit.py -i LogsAsymptoticLimits.txt
