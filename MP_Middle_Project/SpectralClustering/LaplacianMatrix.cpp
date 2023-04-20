#include <vector>
#include <cmath>
#include "SpectralClustring.h"
#include "DS_definitions.h"
#include "DS_timer.h"
#include <omp.h>

#define THREADS (8)

double** generateLaplacianMatrix(double* affinityMatrix[], int n) {
	double** result = new double* [n];
	for (int i = 0; i < n; i++)
	{
		result[i] = new double[n];
		double sum = 0;
		for (int j = 0; j < n; j++)
		{
			if (i != j)
			{
				sum += (affinityMatrix)[i][j];
				result[i][j] = -(affinityMatrix)[i][j];
			}
		}
		result[i][i] = sum;
	}
	for (int i = 0; i < n; i++)
	{
		delete[] affinityMatrix[i];
	}
	delete[] affinityMatrix;
	return result;
}

double** generateLaplacianMatrixParallel(double* affinityMatrix[], int n) {
	double** result = new double* [n];
#pragma omp parallel num_threads(THREADS) 
	{
#pragma omp for
		for (int i = 0; i < n; i++)
		{
			result[i] = new double[n];
			double sum = 0;
			for (int j = 0; j < n; j++)
			{
				if (i != j)
				{
					sum += (affinityMatrix)[i][j];
					result[i][j] = -(affinityMatrix)[i][j];
				}
			}
			result[i][i] = sum;
		}
	}
	for (int i = 0; i < n; i++)
	{
		delete[] affinityMatrix[i];
	}
	delete[] affinityMatrix;
	return result;
}