#!/bin/bash

default=HZaTo2l2g
BASE_DIR="/afs/cern.ch/work/p/pelai/HZa/gridpacks/check_LHE/run2_zebing"
# Handle 0.1 to 0.9,  1 to 30
for m in {1..10} {15,20,25,30}
# for m in $(seq 0.1 0.1 0.9)
# for m in {0.4,2}
do
    # format m file name (Replace . with p)
    m_formatted=$(echo "$m" | tr '.' 'p')

    # echo Create directory M$m
    newdir=./LHEfile/"$default"_M"$m_formatted"
    mkdir -p $newdir
    
    # In case to delete
    # rm -fr ./LHEfile/"$default"_M"$m_formatted"

    # Define source and target file names
    # /afs/cern.ch/work/p/pelai/HZa/gridpacks/Zebing_run2_gridpacks/2017/HZaTo2l2g_M1/v1/HZaTo2l2g_M1_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.tar.xz
    source_file="/afs/cern.ch/work/p/pelai/HZa/gridpacks/Zebing_run2_gridpacks/2017/${default}_M${m_formatted}/v1/${default}_M${m_formatted}_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.tar.xz"
    target_file="${default}_M${m_formatted}_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.tar.xz"

    # # Copy the file
    cp "$source_file" "$newdir/"

    cd "$newdir"

    echo "Untar M$m"
    # Untar the specific file
    tar -xf "$target_file"

    echo "Current directory: $(pwd)"
    echo "generate LHE file for M$m"
    sh runcmsgrid.sh 1000 12345 10

    cd "$BASE_DIR"
done