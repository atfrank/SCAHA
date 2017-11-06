Custom Hungarian Algorithm

Original Source:
https://www.topcoder.com/community/data-science/data-science-tutorials/assignment-problem-and-hungarian-algorithm/

All commands assume you are running in $HOME/R.

Supports doubles (up to two decimal places) and cost matrices with number of columns >= number of rows.

Testing the Hungarian algorithm:

How to Run:
```
# command line

# compile c code, create shared object
R CMD SHLIB hungarian.c hungarian.h
# test hungarian algorithm
Rscript hungarian_c.R   
# Prints the total cost, and time execution
```
Using the Custom Hungarian R code:

```
# command line

# compile c code, create shared object
R CMD SHLIB hungarian.c hungarian.h

# In R file

source('hungarian.R')
assignments <- hungarian(cost_matrix, maximum=FALSE)
# hungarian function is similar to solve_LSAP
```
