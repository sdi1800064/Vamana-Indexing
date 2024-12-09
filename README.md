# Execute

To create the main executable run 
```make```

To create the vamanaFilteringExe executable run~~~~
``` make recallFilteredVamana ```

To install the check library run
```sudo apt-get install check```

To create the test executable run
```make test```

In order to run the main executable 
```./my_program```

In order to run the test executable 
```./test```

In order to run the main executable for the vamanafiltered dataset run
```./recallFilteredVamana -b dummy-data.bin -q dummy-queries.bin --graph graphVamanaFiltered.txt -g groundtruth.ivecs -R 30 -a 1.2 -L 60 -k 100 command line```



you can use the following flags for the main:

```-b <base_file_location>```

```-q <query_file_location>```

```-g <groundtruth_file_location>```

```-a alpha``` float

```-L Lamda``` integer

```-R R``` integer

```-k k``` integer


# Cleaning

To clean the object files run
```make clean```
