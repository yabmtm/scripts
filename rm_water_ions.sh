#!/bin/bash

for i in SOL NA CL PO4 SO4; do
    for i in f_methyl f_non g_methyl g_non; do cd $i; sed -i '/$i/d' ./structure.gro; done; done
