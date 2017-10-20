#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "hungarian.h"
using namespace std;

int main(){
  int n;
  cin >> n;
  vector<vector<int>> a(n, vector<int>(n));
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      cin >> a[i][j];
    }
  }
  Hungarian h(std::move(a), n);
  cout << h.hungarian() << endl;
  h.print_assignment();
}
