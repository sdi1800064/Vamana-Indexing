#!/bin/bash

# Run the first command
echo "Running command 1..."
command1
if [ $? -ne 0 ]; then
    echo "Command 1 failed."
    exit 1
fi

# Run the second command
echo "Running command 2..."
command2
if [ $? -ne 0 ]; then
    echo "Command 2 failed."
    exit 1
fi

# Run the third command
echo "Running command 3..."
command3
if [ $? -ne 0 ]; then
    echo "Command 3 failed."
    exit 1
fi

echo "All commands executed successfully."


./stitchedVamana -b testSets/contest-data-release-1m.bin -R 20 -a 1.2 -L 200 -t 3

