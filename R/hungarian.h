#include <stdio.h>
#include <string.h>
#define PRECISION 0.01
#define N 10000 //max number of vertices in one part
#define INF 100000000 //just infinity
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })


double * cost; //cost matrix
int n, max_match; //n workers and n jobs
double lx[N], ly[N]; //labels of X and Y parts
int xy[N]; //xy[x] - vertex that is matched with x,
int yx[N]; //yx[y] - vertex that is matched with y
int S[N], T[N]; //sets S and T in algorithm
double slack[N]; //as in the algorithm description
int slackx[N]; //slackx[y] such a vertex, that
 // l(slackx[y]) + l(y) - w(slackx[y],y) = slack[y]
int prev[N]; //array for memorizing alternating paths

double get_cost(int x, int y);
void init_labels();

void update_labels();

void add_to_tree(int x, int prevx);
void augment();
void run();
