#include <cmath>
#include <stdio.h>
#include <stdlib.h>

void Transform(double** A, double** AT, int n);
void copy(double** input, double** output, int n);
void swap(double& x1, double& x2);
void sort(double* A, int n);
double norm(double* v, int n);
void matrixMul(double** A, double** B, double** output, int n);
void matrixMul(double** A, double** output, double lambda, int n);
double dots(double* A, double* B, int n);
void QR_Decomposition(double** A, double** Q, double** R, int n);
double getDifference(double** A, double** B, int n);
void findEigenQR(double** A, double* eigenValue, double** eigenVector, int maxIter, int n, double tol=1e-3);