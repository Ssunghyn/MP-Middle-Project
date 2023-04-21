#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 3 // change this value to match the size of your matrix

void matrix_vector_mult(double A[N][N], double x[N], double Ax[N]) {
    int i, j;
    for (i = 0; i < N; i++) {
        Ax[i] = 0.0;
        for (j = 0; j < N; j++) {
            Ax[i] += A[i][j] * x[j];
        }
    }
}

double vector_dot_product(double x[N], double y[N]) {
    int i;
    double result = 0.0;
    for (i = 0; i < N; i++) {
        result += x[i] * y[i];
    }
    return result;
}

void vector_scalar_mult(double x[N], double a) {
    int i;
    for (i = 0; i < N; i++) {
        x[i] *= a;
    }
}

void vector_subtract(double x[N], double y[N], double result[N]) {
    int i;
    for (i = 0; i < N; i++) {
        result[i] = x[i] - y[i];
    }
}

void normalize_vector(double x[N]) {
    double norm = sqrt(vector_dot_product(x, x));
    vector_scalar_mult(x, 1.0 / norm);
}

void rayleigh_quotient_iteration(double A[N][N], double x[N], double* lambda) {
    double Ax[N], x_old[N], Ax_old[N], lambda_old, error;
    int i, iter;
    *lambda = vector_dot_product(x, Ax);
    normalize_vector(x);
    iter = 0;
    do {
        iter++;
        lambda_old = *lambda;
        for (i = 0; i < N; i++) {
            x_old[i] = x[i];
            Ax_old[i] = Ax[i];
        }
        matrix_vector_mult(A, x, Ax);
        *lambda = vector_dot_product(x, Ax);
        vector_scalar_mult(x, *lambda);
        vector_subtract(Ax, x, x);
        normalize_vector(x);
        error = fabs((*lambda - lambda_old) / lambda_old);
    } while (error > 1e-6 && iter < 100);
    printf("Rayleigh quotient iteration converged after %d iterations.\n", iter);
}

int main() {
    double A[N][N] = {{1, 2, 3},
                      {1, 2, 1},
                      {3, 2, 1}};
    double x[N] = {1.0, 1.0, 1.0};
    double lambda;
    rayleigh_quotient_iteration(A, x, &lambda);
    printf("The largest eigenvalue is approximately %f.\n", lambda);
    return 0;
}