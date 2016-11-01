#!/bin/bash

for file in lambda*;do
    cd $file
    echo $file
    tail md.log | grep Performance
    cd ..
done
