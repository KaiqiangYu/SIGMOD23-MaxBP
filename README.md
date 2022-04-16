# SIGMOD23-MaxBP
Our code, data and additional materials are avaliable here.
## Index  
```shell
.
|- README.md
|- Code
|  | - ...
|
|- SIGMOD_MaxBP_report.pdf (technical report)
```


# FastBB
An efficient algorithm for finding the Maximum k-BiPlex (MaxBP).

## Setup
```shell
g++ main.cpp -O3 -o FastBB
```

## Usage
  ./FastBB {OPTIONS}

    FastBB, an efficient algorithm for finding MaxBP

  OPTIONS:

      -h, --help                          Display this help menu
      -f[dataset], --file=[dataset]       Path to dataset
      -k[para k], --k=[para k]            The parameter k
      -u[para theta_l] --u=[para theta_l] The threshold of number of vertices at left side
      -v[para theta_r] --v=[para theta_r] The threshold of number of vertices at right side
      -K[para K], --K=[para K]            Number of MaxBPs to be found


## Input Graph Format
The input graph  should follow the following format.

 Example.graph

    3 1 2
    0 1 2
    1 0
    2 0
    (File ends with an empty line)

(1) Given an input graph G=(L,R), vertices are represented by non-negtive integers from 0 to |V(G)|. By default, {0,1,...,|L|-1} denotes the left side and {|L|,|L|+1,...,|L|+|R|-1} denotes the right side. 

(2) The first line includes 3 non-neigtive integers, e.g., 3 1 2, that denote the number of vertices, the id of the first vertex in the right side and the number of edges, respectively. To illustrate, consider the first line 3 2 1. The input graph has three vertices {0, 1, 2}, and two edges (0, 1) and (0, 2). {0} denotes the left side and {1,2} denotes the right sie.

(3) The following lines represent an adjacent list of the input graph. To illustrate, consider the second line 0 1 2. The vertex with id 0 is adjacent with two vertices 1 and 2.

## Output Format
By default, FastBB returns the MaxBP and the corresponding running time.

    Running Time: 1ms
    |E|: 11  |L|: 3  |R|: 4
    L: 0 2 1
    R: 4 5 6 7



## Running Example

```shell
> g++ main.cpp -O3 -o FastBB
>
> ./FastBB -f "./Example.graph" -k 1 -u 3 -v 3 -K 2
Running Time: 0ms
|E|: 10  |L|: 4  |R|: 3
L: 3 0 2 1
R: 4 6 7

|E|: 11  |L|: 3  |R|: 4
L: 0 2 1
R: 4 5 6 7
>
> cat Example.graph
8 4 13
0 4 5 6 7
1 4 5 6
2 4 5 6 7
3 6 7
4 0 1 2
5 0 1 2
6 0 1 2 3
7 0 2 3

>
>
```
