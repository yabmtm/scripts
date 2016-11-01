#!/bin/bash

# this script will run multiple iterations of preparing an energy minimization, checking where the bad water molecule is each time

# reset

rm -f BADATOMS
for i in 6.8 6.85 6.9 6.95 7 7.05 7.1 7.15 7.2; do # try some different
    for j in 6.8 6.85 6.9 6.95 7 7.05 7.1 7.15 7.2; do # box indices to get
        for k in 6.8 6.85 6.9 6.95 7 7.05 7.1 7.15 7.2; do # variation in where water molecules go
            rm -f box.gro em* ions.tpr mdout.mdp solv.gro step* YIN_minim* topol.top \#* && cp .topol.top topol.top # clean-up from previous runs
            gmx editconf -f conf.gro -o box.gro -box $i $j $k; sleep 1 # prepare and run the minimization
            gmx solvate -cp box.gro -cs spc216.gro -o solv.gro  -p topol.top; sleep 1
            gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr; sleep 1
            echo SOL | gmx genion -s ions.tpr -o solv.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1; sleep 1
            gmx grompp -f minim.mdp -c solv.gro -p topol.top -o em.tpr; sleep 1
            qsub gmx_minim.sh # submit the minimization

            while [ ! -f YIN_minim.out ] # wait until it's done
                do
                echo "minimizing..."
                sleep 5; done
            if [ -f em.gro ]; then echo "DONE" >> BADATOMS; exit; fi # break if it successfully minimizes
            BADATOM=$(cat YIN* | grep settled | sed "s/^..............//" | sed "s/[a-z][A-Z]*//g" | sed "s/ //g" | sed "s/\.//") # otherwise find which atom is preventing it

            MATCH=$(cat solv.gro | grep "$BADATOM" | grep OW) # show which atom it is, as well as its coordinates
            echo $MATCH >> BADATOMS; cat BADATOMS; sleep 2; done; done; done # write the information to a file and start the next trial
