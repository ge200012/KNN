The above is the source code for the paper “GTRSS: Graph-based Top-k Representative Similar
Subtrajectory Query”. Some of the related datasets are under Build, but due to the file size being too large, they have been uploaded to Released. 
The dataset and build files have also been uploaded, and can be downloaded from Tag. 
The above code can be used to test the performance of querying on the Chengdu dataset. Currently, 
the paper is under submission. Please feel free to communicate and exchange ideas with me.

When executing the above code, there are some important parameters that need to be adjusted in order to test different comparative experiments. The following introduces some of the key parameters
1. In the gnns.h
(1) Set the number of neighboring nodes for each node on the graph.
(2) Set True/Flash to indicate whether to record information about nodes retrieved during the retrieval process.
(3) In the execute() function, set the relevant algorithm used to "effieicentAlgorithm" or "exactS", "PSS", "POS“
2. In the knn.h
(1) Set the number of grids $nums_{grid}$ * $nums_{grid}$ for the top-level graph.
(2) Select the method of selecting nodes, $within/near$.
(3) The proportion of different types of neighboring nodes for each node in the top-level graph，$portion-from-front、portion-from-middle、portion-from-back$.
(4) In the naive_struct() function, set the construction method of the r_tree and some of the parameters mentioned above
4. Set the number of neighboring edges $E$ to be searched during the search process in the params.h
5. Set metrics such as "dtw", "FC", "WED" in SimilarTrajectory.cpp
