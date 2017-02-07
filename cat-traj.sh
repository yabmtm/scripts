#!/bin/bash

# This script doesn't work with clones > 100. I'll get there eventually.
# Output can be found in ~/server2/data/SVR166219/TRAJECTORIES.

EXTRACT="FALSE" # extract xtc's, run pbc correction, and sort them
COMBINE="TRUE" # combine each clone and each run depending on whether there are time-step mismatches
NOCHECK="FALSE" # this will combine entire runs without checking for time-step mismatches
group=24 # this references the group in the index file (for my sims this is the custom ligand/protein group)
index_file="ndx.ndx" # change these if incorrect
structure_file="xtc.gro"


while test $# -gt 2 # parse input for arguments
do
    case "$1" in
        PROJ*) proj=$(echo "$1" | tr -dc '0-9') # get just project ID
            ;;
        --extract) EXTRACT=TRUE
            ;;
        --combine) COMBINE=TRUE
            ;;
        --no-check) NOCHECK=TRUE
            ;;
        *) echo "First argument should be project directory, e.g.: PROJ1337"
                echo 'Other arguments include: --extract --combine --no-check'; exit
            ;;
    esac
    shift
done


if [[ $1 == "PROJ"* ]]; then # make sure argument is valid
    proj=$(echo "$1" | tr -dc '0-9')
    work_dir="/home/server/server2/data/SVR166219"
    folders=$(ls ${work_dir}/$1 | grep RUN | wc -l); runs=$(($folders-1)) # set run number
    folders=$(ls ${work_dir}/$1/RUN0 | grep CLONE | wc -l); clones=$(($folders-1)) # set clone number
    echo $runs $clones
    if [ "$EXTRACT" == "TRUE" ]; then

        # Extract all frames correct their periodic boundary conditions, and sort them together in a separate directory

        echo -e "\033[0;32mSET GROUP NUMBER TO BE EXTRACTED FROM XTC BEFORE RUNNING SCRIPT.\033[0m"; sleep 2
        mkdir TRAJECTORIES/${1}_TRAJ
        cd ${work_dir}/$1; echo -e "\033[0;33mPROCESSING ${1}...\033[0m"; sleep 1
        for i in $(seq 0 $runs); do
            cd RUN$i; echo -e "\t\033[0;33mGOING INTO RUN${i}...\033[0m"; sleep 1
            for j in $(seq 0 $clones); do
                if [ -e CLONE$j/results0/traj_comp.xtc ]; then # make sure an xtc exists to be extracted
                    echo -e "\t\t\033[0;33mGOING INTO CLONE${j}...\033[0m"; sleep 0.5
                    frames=$(ls CLONE$j | grep results | wc -l); results=$(($frames-1))
                    for k in $(seq 0 $results); do
                        if [[ $k -lt 10 ]]; then # do a pbc correction on each trajectory and put output in trajectory folder
                            echo $group | gmx trjconv -f CLONE${j}/results${k}/traj_comp.xtc -n ~/server2/projects/Gromacs/p${proj}/ndx.ndx\
                            -s CLONE${j}/frame${k}.tpr -pbc whole -o ../../TRAJECTORIES/${1}_TRAJ/RUN${i}_CLONE${j}_FRAME00${k}.xtc; fi
                        if [[ $k -ge 10 ]] && [[ $k -lt 100 ]]; then
                            echo $group | gmx trjconv -f CLONE${j}/results${k}/traj_comp.xtc -n ~/server2/projects/Gromacs/p${proj}/ndx.ndx\
                            -s CLONE${j}/frame${k}.tpr -pbc whole -o ../../TRAJECTORIES/${1}_TRAJ/RUN${i}_CLONE${j}_FRAME0${k}.xtc; fi
                        if [[ $k -ge 100 ]]; then
                            echo $group | gmx trjconv -f CLONE${j}/results${k}/traj_comp.xtc -n ~/server2/projects/Gromacs/p${proj}/ndx.ndx\
                            -s CLONE${j}/frame${k}.tpr -pbc whole -o ../../TRAJECTORIES/${1}_TRAJ/RUN${i}_CLONE${j}_FRAME${k}.xtc; fi
                    echo -e "\t\t\t\033[0;33mPROCESSED $k FRAMES.\033[0m\n"; done
                else echo -e "\t\t\033[0;33mSKIPPING EMPTY CLONE$j\033[0m"; fi; done; done; fi

    if [ "$COMBINE" == "TRUE" ]; then

        # Combine all frames from each clone and create gmx check log files.

        cd $work_dir
        echo -e "\t\t\033[0;33mCONCATENATING TRAJECTORIES...$j\033[0m"; sleep 2
        cd TRAJECTORIES/${1}_TRAJ; mkdir -p CLONES/{bad_traj,RUNS,tmp} # create all subdirectories at once
        for i in $(seq 0 $runs); do
            for j in $(seq 0 $clones); do
                yes c | gmx trjcat -f RUN${i}_CLONE${j}_FRAME*.xtc -settime -o CLONES/RUN${i}_CLONE${j}.xtc # combine all frames for each clone
                echo "CHECKING RUN${i}, CLONE${j} FOR ERRORS."
                if [[ $j -lt 10 ]]; then # make a gmx check log file for each clone
                    gmx check -f CLONES/RUN${i}_CLONE${j}.xtc &> CLONES/RUN${i}_CLONE00${j}.log; fi
                if [[ $j -ge 10 ]] && [[ $j -lt 100 ]]; then
                    gmx check -f CLONES/RUN${i}_CLONE${j}.xtc &> CLONES/RUN${i}_CLONE0${j}.log; fi; done
                if [[ $j -ge 100 ]]; then
                    gmx check -f CLONES/RUN${i}_CLONE${j}.xtc &> CLONES/RUN${i}_CLONE${j}.log; fi; done
            for j in {0..9}; do mv CLONES/RUN${i}_CLONE${j}.xtc CLONES/RUN${i}_CLONE00${j}.xtc; done
            for j in {10..99}; do mv CLONES/RUN${i}_CLONE${j}.xtc CLONES/RUN${i}_CLONE0${j}.xtc done; done
            echo "FINISHED COMBINING CLONES FROM RUN${i}."; sleep 2
        for i in $(seq 0 $runs); do
            yes c | gmx trjcat -f CLONES/RUN${i}_CLONE*.xtc -settime -o CLONES/RUNS/RUN${i}.xtc
            gmx check -f CLONES/RUNS/RUN${i}.xtc &> CLONES/RUNS/RUN${i}.log; done; fi

        # Check all combined clones for mismatch timestep errors

        cd ${work_dir}/TRAJECTORIES/${1}_TRAJ/CLONES

        for i in *log; do # find all bad clones and move them into separate directory
            count=$( cat $i | grep match | wc -l)
            if [[ $count -ge 1 ]]; then
                file=`echo "$i" | cut -d'.' -f1`
                mv ${file}* bad_traj; fi; done

        for i in $(seq 0 $runs); do
            begin=0
            for j in $(seq 0 $clones); do
                if [[ $j -lt 10 ]]; then file="RUN${i}_CLONE00${j}.xtc"; fi
                if [[ $j -ge 10 ]] && [[ $j -lt 100 ]]; then file="RUN${i}_CLONE0${j}.xtc"; fi
                if [[ $j -gt 100 ]]; then file="RUN${i}_CLONE${j}.xtc"; fi
                if [ ! -f $file ]; then # copy each set of good xtc's over one at a time
                    for k in $(seq $begin $j); do
                        if [[ $begin -lt 10 ]]; then
                            cp RUN${i}_CLONE00${k}.xtc tmp; fi
                        if [[ $begin -ge 10 ]] && [[ $begin -lt 100 ]]; then
                            cp RUN${i}_CLONE0${k}.xtc tmp; fi
                        if [[ $begin -ge 100 ]]; then
                            cp RUN${i}_CLONE0${k}.xtc tmp; fi; done
                    gmx trjcat -f tmp/* -o RUNS/RUN${i}_CLONES_${begin}-$(($j-1)).xtc # concatenate each good set
                    rm -f tmp/*

                    begin=$(($j+1)); fi; done

            last_clone=$(ls | grep CLONE | sed "s/^.*CLONE//" | sed "s/\..*$//" | tail -n 1) # maybe last valid clone is not 99

            if [[ $begin -ne $last_clone ]]; then # concatenate the last good set
                for k in $(seq $begin $last_clone); do
                    if [[ $begin -lt 10 ]]; then
                        cp RUN${i}_CLONE00${k}.xtc tmp; fi
                    if [[ $begin -ge 10 ]] && [[ $begin -lt 100 ]]; then
                        cp RUN${i}_CLONE0${k}.xtc tmp; fi
                    if [[ $begin -ge 100 ]]; then
                        cp RUN${i}_CLONE${k}.xtc tmp; fi; done
                gmx trjcat -f tmp/* -o RUNS/RUN${i}_CLONES_${begin}-${clones}.xtc; fi; done

#        rm -rf tmp # clean-up
        cp ~/server2/projects/Gromacs/p${proj}/${structure_file} RUNS; fi # for convenience

if [ -z "$1" ] || [[ $1 != "PROJ"* ]]; then # if bad argument is passed
    echo "Invalid argument. Run as: $0 PROJ1337"
    echo "Alternatively, you can run $0 -h"; fi

