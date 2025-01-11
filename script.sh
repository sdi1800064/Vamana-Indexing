#!/bin/bash

commands=()
for t in {5..1}; do
    commands+=("./stitchedVamana -b testSets/dummy-data.bin -R 60 -a 1.2 -L 120 -t $t")
    commands+=("./recallSticthedVamana -b testSets/dummy-data.bin -q testSets/dummy-queries.bin -g testSets/groundtruth.ivecs -R 60 -a 1.2 -L 120 -k 100 -t $t --graph graphs/stitchedGraph")
    commands+=("./stitchedVamana -b testSets/dummy-data.bin -R 20 -a 1.2 -L 200 -t $t")
    commands+=("./recallSticthedVamana -b testSets/dummy-data.bin -q testSets/dummy-queries.bin -g testSets/groundtruth.ivecs -R 20 -a 1.2 -L 200 -k 100 -t $t --graph graphs/stitchedGraph")
    
    commands+=("./stitchedVamana -b testSets/contest-data-release-1m.bin -R 60 -a 1.2 -L 120 -t $t")
    commands+=("./recallSticthedVamana -b testSets/contest-data-release-1m.bin -q testSets/contest-queries-release-1m.bin -g testSets/neighbors1m.ivecs -R 60 -a 1.2 -L 120 -k 100 -t $t --graph graphs/stitchedGraph")
    commands+=("./stitchedVamana -b testSets/contest-data-release-1m.bin -R 20 -a 1.2 -L 200 -t $t")
    commands+=("./recallSticthedVamana -b testSets/contest-data-release-1m.bin -q testSets/contest-queries-release-1m.bin -g testSets/neighbors1m.ivecs -R 20 -a 1.2 -L 200 -k 100 -t $t --graph graphs/stitchedGraph")
done

for i in "${commands[@]}"; do
    echo "Running command: $i"
    $i
    if [ $? -ne 0 ]; then
        echo "Command failed."
        exit 1
    fi
done

echo "All commands executed successfully."

