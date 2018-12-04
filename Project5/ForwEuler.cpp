#include <ForwEuler.hpp>

void FWSolver(int n, double alpha, int tmax){
        vec u = zeros<vec>(n);
        vec unew = zeros<vec>(tmax);

        // Boundary condition (u(0) set by zeros)
        u(n) = 1.0;
        unew_t(n) = 1.0;

        for(int j = 1; j < tmax; j++) {
                for(int i = 1; i < n; i++) {
                        u(i) = (1.0 - 2.0*alpha)*u(i) + alpha*u(i+1) + alpha*u(i-1);
                }
        }
}
