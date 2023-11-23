/**
 * The code in this file was obtained from this repository:
 * https://github.com/hmatuschek/eigen3-nnls/
 * and is a port of Fortran code made freely available on netlib.
 * It is now part of Eigen3, which is released under the Mozilla
 * Public License 2 (https://www.mozilla.org/en-US/MPL/2.0/)
 */
#ifndef NNLS_H
#define NNLS_H
#ifdef __cplusplus
extern "C"{
#endif
int nnls_c(double* a, const int* mda, const int* m, const int* n, double* b,
    double* x, double* rnorm, double* w, double* zz, int* index,
    int* mode);
#ifdef __cplusplus
}
#endif
#endif
