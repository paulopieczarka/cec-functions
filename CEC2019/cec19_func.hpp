#pragma once

#include <iostream>
#include <cstdlib>
#include <math.h>

#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

class Cec19 {
  public:
    Cec19();
    ~Cec19();

    void Lennard_Jones_func(double *, int, double *); /* Lennard Jones */
    void Hilbert_func(double *, int, double *); /* Hilbert */
    void Chebyshev_func(double *, int, double *); /* Chebyshev */
    void Schaffer_F7_func(double *, double *, int, double *, double *, int, int); /* Schwefel's F7 */
    void Ackley_func(double *, double *, int, double *, double *, int, int); /* Ackley's */
    void Rastrigin_func(double *, double *, int, double *, double *, int, int); /* Rastrigin's */
    void Step_rastrigin_func (double *, double *, int, double *, double *, int, int); /* Noncontinuous Rastrigin's */
    void Weierstrass_func(double *, double *, int, double *, double *, int, int); /* Weierstrass's */
    void Schwefel_func(double *, double *, int, double *, double *, int, int); /* Schwefel's */
    void EScaffer6_func(double *, double *, int, double *, double *, int, int); /* Expanded Scaffer's F6 */
    void Happycat_func(double *, double *, int, double *, double *, int, int); /* HappyCat */
    void Griewank_func(double *, double *, int, double *, double *, int, int); /* Griewank's */

    void test_func(double *, double *, int, int, int);
    void free_func();

  private:
    double *OShift;
    double *M;
    double *y;
    double *z;
    double *x_bound;

    int ini_flag;
    int n_flag;
    int func_flag;
    int *SS;

    void shift_func(double*, double*, int, double*);
    void rotate_func(double*, double*, int, double*);
    void sr_func(double *, double *, int, double*, double*, double, int, int); /* shift and rotate */
    void asy_func(double *, double *x, int, double);
    void osz_func(double *, double *, int);
};
