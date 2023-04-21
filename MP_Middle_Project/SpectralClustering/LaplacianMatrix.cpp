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

double** generateLaplacianMatrixParallelCase1(double* affinityMatrix[], int n) {
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

double** generateLaplacianMatrixParallelCase2(double* affinityMatrix[], int n) {
	double** result = new double* [n];
	for (int i = 0; i < n; i++)
	{
		result[i] = new double[n];
		double sum = 0;
		#pragma omp parallel num_threads(THREADS)
		{
			#pragma omp for reduction(+:sum)
			for (int j = 0; j < n; j++)
			{
				if (i != j)
				{
					sum += (affinityMatrix)[i][j];
					result[i][j] = -(affinityMatrix)[i][j];
				}
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

double** generateLaplacianMatrixParallelCase3(double* affinityMatrix[], int n) {
	double** result = new double* [n];
	for (int i = 0; i < n; i++)
	{
		result[i] = new double[n];
		double sum = 0;
		#pragma omp parallel num_threads(THREADS)
		{
			#pragma omp for schedule(static, (n / THREADS)) nowait
			for (int j = 0; j < n; j++)
			{
				if (i != j) result[i][j] = -(affinityMatrix)[i][j];
			}
			#pragma omp for schedule(static, n / (THREADS)) reduction(+:sum)
			for (int j = 0; j < n; j++)
			{
				sum += affinityMatrix[i][j];
			}
		}
		result[i][i] = sum - (affinityMatrix)[i][i];
	}
	for (int i = 0; i < n; i++)
	{
		delete[] affinityMatrix[i];
	}
	delete[] affinityMatrix;
	return result;
}

double** generateLaplacianMatrixParallelCase4_1(double* affinityMatrix[], int n) {
	int threadCount1 = THREADS / 2;
	int	threadCount2 = THREADS - threadCount1;
	double** result = new double* [n];
	omp_set_nested(1);
	for (int i = 0; i < n; i++)
	{
		result[i] = new double[n];
		#pragma omp parallel num_threads(2)
		{
			#pragma omp single nowait
			{
				#pragma omp parallel num_threads(threadCount1)
				{
					#pragma omp for
					for (int j = 0; j < n; j++)
					{
						if (i != j) result[i][j] = -(affinityMatrix)[i][j];
					}
				}
			}
			#pragma omp single
			{
				double sum = 0;
				#pragma omp parallel num_threads(threadCount2)
				{
					int id = omp_get_thread_num();
					#pragma omp for reduction(+:sum)
					for (int j = 0; j < n; j++)
					{
						if (i != j) sum += (affinityMatrix)[i][j];
					}
				}
				result[i][i] = sum;
			}
		}
	}
	omp_set_nested(0);
	for (int i = 0; i < n; i++)
	{
		delete[] affinityMatrix[i];
	}
	delete[] affinityMatrix;
	return result;
}

double** generateLaplacianMatrixParallelCase4_2(double* affinityMatrix[], int n) {
	double** result = new double* [n];
	int threadCount1 = THREADS / 2;
	int threadCount2 = THREADS - threadCount1;
	for (int i = 0; i < n; i++)
	{
		result[i] = new double[n];
		double sum = 0;
		#pragma omp parallel num_threads(THREADS) reduction(+:sum)
		{
			int tid = omp_get_thread_num();
			if (tid < threadCount1)
			{
				int start = tid * n / threadCount1;
				int end = (tid + 1) * n / threadCount1;
				if (tid == threadCount1 - 1) end = n;
				for (int j = start; j < end; j++)
				{
					result[i][j] = -(affinityMatrix)[i][j];
				}
			}
			else {
				tid -= threadCount1;
				int start = tid * n / threadCount2;
				int end = (tid + 1) * n / threadCount2;
				if (tid == threadCount2 - 1) end = n;
				for (int j = start; j < end; j++)
				{
					if (i != j) sum += (affinityMatrix)[i][j];
				}
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