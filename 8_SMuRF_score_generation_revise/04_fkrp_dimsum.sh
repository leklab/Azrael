#!/bin/bash

# Loop through blocks 1 to 5
for i in {1..5}; do
    projectName="FKRP_$i"
    countPath="./fkrp_variantCounts_block$i.txt"
    wildtypeSequence=$(awk 'NR==2 {print $1}' "fkrp_variantCounts_block$i.txt")

    ./DiMSum --projectName "$projectName" --experimentDesignPath ./fkrp_experimentDesign.txt --startStage 4 --numCores 12 --countPath "$countPath" --wildtypeSequence "$wildtypeSequence" --sequenceType noncoding
done

# Process separately for FKRP block 6 (which exludes FKRP-0 sample)
./DiMSum --projectName FKRP_6 --experimentDesignPath ./fkrp_experimentDesign_block6.txt --startStage 4 --numCores 12 --countPath ./fkrp_variantCounts_block6.txt --sequenceType noncoding --wildtypeSequence "$(awk 'NR==2 {print $1}' fkrp_variantCounts_block6.txt)"
