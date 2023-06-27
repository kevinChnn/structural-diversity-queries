# structural-diversity-queries
## Dataset
### Dataset Sources
[SNAP](https://snap.stanford.edu/) and the [KONECT](http://konect.cc/) Project.

### Dataset Format
First Line:         (n, m, t)  
                    [where n is (maximum vertex ID + 1), m is the number of edges, and t is the maximum edge timestamp]  
Following Lines:    (u, v, t)  
                    [where u and v are vertices and u < v]  

We have included a sample dataset, CollegeMsg-r.txt, which is procssed from the original CollegeMsg dataset from [SNAP](https://snap.stanford.edu/).

## Compile & Run
### Compile
```
g++ -O 3 main.cpp -o main
```

### Run
```
./main
usage: ./main [setting] [method] [interval] [graph] [tau] [theta](optional)
[setting] 1 -> Arbitrary Window Query
[setting] 2 -> Sliding Window Query
[method] 1 -> Online
[method] 2 -> Baseline
[method] 3 -> Our
```
The [setting] refers to the Arbitrary/Sliding Window Query.  
The [method] refers to the Online/Baseline/Our solution.  
The [interval] refers to the time label unit (3600s->1h).  
The [graph] refers to the dataset.  
The [tau] refers to the size threshold.  
The [theta] refers to the ratio of the sliding window size and the ratio of the query window size.  

### Examples:
Arbitrary Window Query (Online)
```
./main 1 1 3600 CollegeMsg-r.txt 2
```

Arbitrary Window Query (Our)
```
./main 1 2 3600 CollegeMsg-r.txt 2
```

Arbitrary Window Query (Online)
```
./main 1 3 3600 CollegeMsg-r.txt 2
```

Sliding Window Query (Online)
```
./main 2 1 3600 CollegeMsg-r.txt 2
```

Sliding Window Query (Our)
```
./main 2 2 3600 CollegeMsg-r.txt 2
```

Sliding Window Query (Online)
```
./main 2 3 3600 CollegeMsg-r.txt 2
```
