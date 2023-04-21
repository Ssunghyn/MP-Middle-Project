#include <cmath>
#include <stdio.h>
#include <stdlib.h>

void MultiTransform(double** A, double** AT, int n);
void Multicopy(double** input, double** output, int n);
double Multinorm(double* v, int n);
void MultimatrixMul(double** A, double** B, double** output, int n);
double Multidots(double* A, double* B, int n);
void MultiQR_Factorization(double** A, double** Q, double** R, double** H, int n);
void MultifindEigenQR(double** A, double* eigenValue, double** eigenVector, int maxIter, int n, double tol = 1e-2);