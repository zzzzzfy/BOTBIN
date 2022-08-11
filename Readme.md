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

For example, you can download the dataset of twitter from http://konect.cc/networks/twitter/. Then you need to change the user_name and create the appropriate folder to ensure that the dataset is in the same path as in *Config.h* (you can also specify your own file path). Then change the downloaded twitter file to out.twitter and make sure its path is the same as *src_* in *Config.h* under the condition that the *file_id* is 0. The *file_id* is 0 corresponds to the twitter dataset, and other datasets can be configured in the same way. Then just run the following command to clean the dataset corresponding to *file_id*.

Clean the dataset {file_id}
```sh
$ ./botbin -action 0 -file {file_id}
```
After cleaning the dataset, you can run the following command directly. Build the BOTBIN of dataset {file_id} with the paramet $\rho$ and $\delta$.
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


## Acknowledgements

Thanks to the author of GS-index for some code.