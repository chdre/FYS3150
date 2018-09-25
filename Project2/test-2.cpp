#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "eigvals.h"

TEST_CASE("Max value A(i,j)"){
        /* testing that the function for finding the max off-diagonal element
           works. Creating a toeplitz matrix with random elements. Comparing the
           result from function max_offdiag with the max value found using
           Armadillo's "max()". */
        int k, l;   // declaring k, l for max_offdiag.
        int n = 8;
        vec a = randu<vec>(n);
        mat A_test = toeplitz(a);
        A_test.diag().fill(0);    // setting the diagonal to 0

        double actualmax = A_test.max();
        double maxelm = max_offdiag(A_test, n, &l, &k);

        REQUIRE(maxelm==Approx(actualmax));
}

TEST_CASE("Eigenvalues"){
        /* testing that the function for finding eigenvals using Jacobi's method
           and comparing the solution with Armadillo's eig_sym. Creating a random
           Toeplitz matrix. Finally checking that the values are correct by
           looping over the eigenvalues.*/
        int n = 8;
        double h = 1.0/double(n);

        vec a = randu<vec>(n);
        mat A_test = toeplitz(a);
        A_test.diag().fill(0);    // setting the diagonal to 0

        mat R_test; // matrix for eigenvectors needed for Jacobi's method

        vec eigvals_jac = jacobi(A_test, R_test, n, h);
        eigvals_jac = sort(eigvals_jac);  // sorting the values by size

        // finding eigenvalues with Armadillo
        vec eigval;
        mat eigvec;
        eig_sym(eigval, eigvec, A_test);

        for(int i = 0; i < n; i++) {
                // looping over eigenvalues and comparing
                REQUIRE(eigval[i]==Approx(eigvals_jac[i]));
        }


}
