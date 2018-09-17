#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  int N; // size of matrix

  if (argc <= 1) {
    cout << "Missing arguments:" << argv[0] << " specify size of matrix N"
         << endl;
    exit(1);
  } else {
    // reading value of N
    N = atoi(argv[1]) // N from acsi
  }

  double h = 1.0 / (N + 1);
}
