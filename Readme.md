# Effective Indexing for Dynamic Structural Graph Clustering

This repository implements the Bottom-$k$ and bucket indexing (BOTBIN) data structure in paper *Effective Indexing for Dynamic Structural Graph Clustering*
. It can support the construction of similarity index and bucket index from a graph, so as to quickly obtain approximate results of graph clustering. And for inserting or deleting edges, it can support fast modification, and return good clustering results with theoretical guarantees.

## Compile the code
Please note that the path of the input dataset needs to be configured by modifying the code. After you have configured the file path, you can execute the following commands to generate the executable file *botbin*.

```sh
$ cd src
$ make clean
$ make
```
After compiling the code, an executable file called *botbin* is generated.

## Run the code

Clean the dataset {file_id}
```sh
$ ./botbin -action 0 -file {file_id}
```

Build the BOTBIN of dataset {file_id} with the paramet $\rho$ and $\delta$.
```sh
$ ./botbin -action 1 -file {file_id} -rho {rho} -delta {delta} -build
```
Evenly sample some edges for update test.
```sh
$ ./botbin -action 4 
```
Update test under BOTBIN with specific parameters.
```sh
$ ./botbin -action 1 -file {file_id} -rho {rho} -delta {delta} 
```

After configuring the query parameters for *query_config* manually, you can generate the query results by invoking the following command:
```sh
$ ./botbin -action 1 -file {file_id} -rho {rho} -delta {delta} -query 
```

*Read the code for more details.*


## Data format
[vertex_id1] [vertex_id2]

[vertex_id1] [vertex_id2]

[vertex_id1] [vertex_id2]

...

## Acknowledgements

Thanks to the author of GS-index for some code.