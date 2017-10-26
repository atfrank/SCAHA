Custom Hungarian Algorithm

Original Source:
https://www.topcoder.com/community/data-science/data-science-tutorials/assignment-problem-and-hungarian-algorithm/

All commands assume you are running in $HOME/R

How to Run:
```
R CMD SHLIB hungarian.c hungarian.h
Rscript hungarian_c.R    
```

C++ version (requires text input):
```
time echo n "$(cat test.txt)" | cpp_version/./hungarian
```
Where n = the dimension of the array.

For now all values are integers and cost matrix must be nxn.
