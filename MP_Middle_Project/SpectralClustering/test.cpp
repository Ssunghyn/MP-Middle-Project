#include <vector>
#include "DS_definitions.h"
#include "DS_timer.h"
#include <omp.h>
#include "DS_timer.h"
#include "DS_definitions.h"

double** generateLaplacianMatrix(double *affinityMatrix[], int n);
double** generateLaplacianMatrixParallel(double* affinityMatrix[], int n);
double** generateLaplacianMatrixParallel2(double* affinityMatrix[], int n);
double** generateLaplacianMatrixParallel3(double* affinityMatrix[], int n);
double** generateLaplacianMatrixParallel4(double* affinityMatrix[], int n);
double** generateLaplacianMatrixParallel5(double* affinityMatrix[], int n);

#define GenRandom ((double)(rand() % 100) / 10)
#define SIZE (10000)
#define THREADS (8)

std::vector<std::vector<double> > sample(SIZE);

int main(int argc, char** argv) {
	DS_timer timer(6);
	timer.setTimerName(0, (char*)"Serial");
	timer.setTimerName(1, (char*)"Parallel Custom nested");
	timer.setTimerName(2, (char*)"Parallel Custom thread indexing");
	timer.setTimerName(3, (char*)"Parallel Custom two for");
	timer.setTimerName(4, (char*)"Parallel Simple for");
	timer.setTimerName(5, (char*)"Parallel Simple inner for");


	for (int i = 0; i < SIZE; i++) {
		std::vector<double> temp(SIZE);
		for (int j = 0; j < SIZE; j++)
		{
			temp[j] = GenRandom;
		}
		sample[i] = temp;
	}
	double** data = new double* [SIZE];
	for (int i = 0; i < SIZE; i++)
	{
		data[i] = new double[SIZE];
		for (int j = 0; j < SIZE; j++)
		{
			data[i][j] = sample[i][j];
		}
	}
	timer.onTimer(0);
	double** result = generateLaplacianMatrix(data, SIZE);

	timer.offTimer(0);
	for (int i = 0; i < SIZE; i++)
	{
		delete[] result[i];
	}
	delete[] result;

	double** data2 = new double* [SIZE];
	for (int i = 0; i < SIZE; i++)
	{
		data2[i] = new double[SIZE];
		for (int j = 0; j < SIZE; j++)
		{
			data2[i][j] = sample[i][j];
		}
	}
	timer.onTimer(1);

	result = generateLaplacianMatrixParallel(data2, SIZE);
	timer.offTimer(1);
	for (int i = 0; i < SIZE; i++)
	{
		delete[] result[i];
	}
	delete[] result;

	double** data3 = new double* [SIZE];
	for (int i = 0; i < SIZE; i++)
	{
		data3[i] = new double[SIZE];
		for (int j = 0; j < SIZE; j++)
		{
			data3[i][j] = sample[i][j];
		}
	}
	timer.onTimer(2);

	result = generateLaplacianMatrixParallel2(data3, SIZE);
	timer.offTimer(2);
	for (int i = 0; i < SIZE; i++)
	{
		delete[] result[i];
	}
	delete[] result;


	double** data4 = new double* [SIZE];
	for (int i = 0; i < SIZE; i++)
	{
		data4[i] = new double[SIZE];
		for (int j = 0; j < SIZE; j++)
		{
			data4[i][j] = sample[i][j];
		}
	}
	timer.onTimer(3);

	result = generateLaplacianMatrixParallel3(data4, SIZE);
	timer.offTimer(3);
	for (int i = 0; i < SIZE; i++)
	{
		delete[] result[i];
	}
	delete[] result;

	double** data5 = new double* [SIZE];
	for (int i = 0; i < SIZE; i++)
	{
		data5[i] = new double[SIZE];
		for (int j = 0; j < SIZE; j++)
		{
			data5[i][j] = sample[i][j];
		}
	}
	timer.onTimer(4);

	result = generateLaplacianMatrixParallel4(data5, SIZE);
	timer.offTimer(4);
	for (int i = 0; i < SIZE; i++)
	{
		delete[] result[i];
	}
	delete[] result;

	double** data6 = new double* [SIZE];
	for (int i = 0; i < SIZE; i++)
	{
		data6[i] = new double[SIZE];
		for (int j = 0; j < SIZE; j++)
		{
			data6[i][j] = sample[i][j];
		}
	}
	timer.onTimer(5);

	result = generateLaplacianMatrixParallel5(data6, SIZE);
	timer.offTimer(5);
	for (int i = 0; i < SIZE; i++)
	{
		delete[] result[i];
	}
	delete[] result;

	timer.printTimer();


	EXIT_WIHT_KEYPRESS;
}

double** generateLaplacianMatrix(double *affinityMatrix[], int n) {
	double** result = new double* [n];
	for (int i = 0; i < n; i++)
	{
		result[i] = new double [n];
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
						sum += (affinityMatrix)[i][j];
					}
				}
				result[i][i] = sum - (affinityMatrix)[i][i];
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		delete[] affinityMatrix[i];
	}
	delete[] affinityMatrix;
	return result;
}

double** generateLaplacianMatrixParallel2(double* affinityMatrix[], int n) {
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
					sum += (affinityMatrix)[i][j];
				}
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

double** generateLaplacianMatrixParallel3(double* affinityMatrix[], int n) {
	double** result = new double* [n];
	for (int i = 0; i < n; i++)
	{
		result[i] = new double[n];
		double sum = 0;
		#pragma omp parallel num_threads(THREADS)
		{
			#pragma omp for schedule(static, (n / THREADS))
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

double** generateLaplacianMatrixParallel4(double* affinityMatrix[], int n) {
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

double** generateLaplacianMatrixParallel5(double* affinityMatrix[], int n) {
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