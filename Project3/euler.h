#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class euler {

public:

double h;   // step length
double vel, pos, accel;
euler(double accel);

}
