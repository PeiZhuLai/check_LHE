#!/bin/bash

default=HZaTo2l2g
BASE_DIR="/afs/cern.ch/work/p/pelai/HZa/gridpacks/check_LHE/run2_zebing"
eos_DIR="/eos/home-p/pelai/HZa/ALP/gridpacks/check_LHE/run2/rootfile"

# Handle 0.1 to 0.9,  1 to 30
for m in {1..10} {15,20,25,30}
# for m in 2
do
    # format m file name (Replace . with p)
    m_formatted=$(echo "$m" | tr '.' 'p')

    # echo Create directory M$m
    # outputdir="$eos_DIR"/"$default"_M"$m_formatted"
    # mkdir -p $outputdir

    # Define source and target file names
    input_lhe="$BASE_DIR"/"LHEfile/${default}_M${m_formatted}/cmsgrid_final.lhe"
    output_rootfile="$eos_DIR"/"ALP_M${m_formatted}.root"

    echo "$m_formatted"
    python ./LHEReader/LHEReader.py --input "$input_lhe" --output "$output_rootfile"

done