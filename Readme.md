# Effective Indexing for Dynamic Structural Graph Clustering

This repository implements the Bottom-$k$ and bucket indexing (BOTBIN) data structure in paper *Effective Indexing for Dynamic Structural Graph Clustering*
. It can support the construction of similarity index and bucket index from a graph, so as to quickly obtain approximate results of graph clustering. And for inserting or deleting edges, it can support fast modification, and return good clustering results with theoretical guarantees.

## Compile the code

You can compile the code any way you like.

## Run the code
1. Please modify the path and file name for the program to read the data correctly.
2. First use **format** to clean the original image
3. Use **BOTBIN** to construct index
4. You can directly generate query parameters, and then make the program **query** to get the clustering results, and the clustering results will be output to each file separately.
5. You can use **sample_edges** to generate modified edges, and then set appropriate parameters to dynamically modify the index.
6. *Read the code for more details.*


## Data format
[vertex_id1] [vertex_id2]

[vertex_id1] [vertex_id2]

[vertex_id1] [vertex_id2]

...

## Acknowledgements

Thanks to the author of GS-index for some code.