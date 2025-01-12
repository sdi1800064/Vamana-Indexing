## Καραγιάννη Θωμάς sdi1800064@di.uoa.gr - 1115201800064
## Μαυραΐδης Κωνσταντίνος sdi0700101@di.uoa.gr - 1115200700101

# Execute

## Installation of test Framework
To install the check library run
```sudo apt-get install check```


## Executables
To create the test executable run
```make test```

To create the recall programm for Filter Vamana  run
``` make recallf ```

To create the Filtered Vamana executable that creates the Filtered Vamana Graphs run
```make createFilteredVamana```

To create the Recall program for the Stitched Vamana run
``` make recallst ```

To create the Stitched Vamana executable that creates the Stitched Vamana Graphs run
``` make stitched ```

## Run

### create filtered Vamana
In order to run the create for  Filtered Vamana sample command
```./createFilteredVamana -b contest-data-release-1m.bin --graph filteredVamana1mParallel4ThreadsL120R60Testdsa -L 120 -a 1.2 -R 60 -t 3```

### Recall
In order to run the Recall program of Filter Vamana run(sample command)
#### ex:
```./recallFilteredVamana -b contest-data-release-1m.bin -q contest-queries-release-1m.bin --graph filteredVamana1mParallel4Threads -g neighbors1m.ivecs -R 60 -a 1.2 -L 120 -k 100 -t 4``


### create stiched  Vamana
In order to run the create for  Stiched Vamana sample command
```./stitchedVamana -b contest-data-release-1m.bin -a 1.2 -R 60 -L 120 -t 4``


Similarly, to run the Recall program of Stitched Vamana run sample command
#### ex:
```./recallStitchedVamana -b testSets/dummy-data.bin -q testSets/dummy-queries.bin -g testSets/groundtruth.ivecs --graph stitchedGraph -R 30 -a 1.3 -k 100 -L 200 -t 4```

You need to use the following flags:

```-b <base_file_location>```

```-q <query_file_location>```

```-g <groundtruth_file_location>```

```--graph <graph_file_location>```

```-a alpha``` float

```-L Lamda``` integer

```-R R``` integer

```-k k``` integer

The arguments a and R are mandatory in case the graph file doesnt exist. In this case the program will create the graph(s) needed and save the graph(s) on this path ```<graph_file_location_R<r>.bin```

### Graph Creation
#### Stitched Graphs
To create the Sticthed Graphs you can either use the ```recallStitchedVamana``` executable or the ```stitchedVamana``` one.

To run the ```stitchedVamana``` you need to use the following arguments:

```-b <base_file_location>```

```-a alpha``` float

```-L Lamda``` integer

```-R R``` integer

In this case the graph will be stored in a file named __stitchedGraph_R<r>.bin__ 

#### Filtered Graph
The Filtered graph can be created by running the ```recallFilteredVamana``` executable.


# Cleaning

To clean the object files run
```make clean```
