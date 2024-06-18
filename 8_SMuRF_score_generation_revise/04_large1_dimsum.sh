#!/bin/bash

# Loop through blocks 1 to 10
for i in {1..10}; do
    projectName="LARGE1_$i"
    countPath="./large1_variantCounts_block$i.txt"
    wildtypeSequence=$(awk 'NR==2 {print $1}' "large1_variantCounts_block$i.txt")

    ./DiMSum --projectName "$projectName" --experimentDesignPath ./large1_experimentDesign.txt --startStage 4 --numCores 12 --countPath "$countPath" --wildtypeSequence "$wildtypeSequence" --sequenceType noncoding
done
