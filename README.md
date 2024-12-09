# Execute

## Installation of test Framework
To install the check library run
```sudo apt-get install check```


## Executables
To create the test executable run
```make test```

To create the Filter Vamana executable run
``` make recallf ```

To create the Recall program for the Stitched Vamana run
``` make recalls ```

To create the Stitched Vamana executable that creates the Stitched Vamana Graphs run
``` make stitched ```

## Run
### Recall
In order to run the Recall program of Filter Vamana run
#### ex:
```./recallFilteredVamana -b testSets/dummy-data.bin -q testSets/dummy-queries.bin -g testSets/groundtruth.ivecs --graph filterGraph -R 30 -a```

Similarly, to run the Recall program of Stitched Vamana run
#### ex:
```./recallStitchedVamana -b testSets/dummy-data.bin -q testSets/dummy-queries.bin -g testSets/groundtruth.ivecs --graph stitchedGraph -R 30 -a 1.3 -k 100 -L 200```

You need to use the following flags:

```-b <base_file_location>```

```-q <query_file_location>```

```-g <groundtruth_file_location>```

```--graph <graph_file_location>```

```-a alpha``` float

```-L Lamda``` integer

```-R R``` integer

```-k k``` integer

The arguments a and R are mandatory in order the graph file doesnt exist. In this case the program will create the graph(s) needed and save the graph(s) on a file with this name ```<graph_R<r>.bin```

### Graph Creation
To create the Sticthed Graphs you can either use the ```recall``` executable or the ```stitchedVamana``` one.

To run the ```stitchedVamana``` you need to use the following arguments:

```-b <base_file_location>```

```-a alpha``` float

```-L Lamda``` integer

```-R R``` integer

In this case the graph will be stored in a file named __stitchedGraph_R<r>.bin__ 



# Cleaning

To clean the object files run
```make clean```
