#include <vector>
#include "DS_definitions.h"
#include "DS_timer.h"
#include <omp.h>
#include "DS_timer.h"
#include "DS_definitions.h"

double** _generateLaplacianMatrix(double *affinityMatrix[], int n);
double** _generateLaplacianMatrixParallelCase1(double* affinityMatrix[], int n);
double** _generateLaplacianMatrixParallelCase2(double* affinityMatrix[], int n);
double** _generateLaplacianMatrixParallelCase3(double* affinityMatrix[], int n);
double** _generateLaplacianMatrixParallelCase4_1(double* affinityMatrix[], int n);
double** _generateLaplacianMatrixParallelCase4_2(double* affinityMatrix[], int n);

#define GenRandom ((double)(rand() % 100) / 10)
#define SIZE (1000)
#define THREADS (8)

std::vector<std::vector<double>> sample(SIZE);

int lapacianTestMain()
//int main(int argc, char** argv) 
{
	DS_timer timer(6);
	timer.setTimerName(0, (char*)"Serial");
	timer.setTimerName(1, (char*)"Parallel Case 1");
	timer.setTimerName(2, (char*)"Parallel Case 2");
	timer.setTimerName(3, (char*)"Parallel Case 3");
	timer.setTimerName(4, (char*)"Parallel Case 4_1");
	timer.setTimerName(5, (char*)"Parallel Case 4_2");


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
	double** result = _generateLaplacianMatrix(data, SIZE);

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

	result = _generateLaplacianMatrixParallelCase1(data2, SIZE);
	timer.offTimer(1);
	/*for (int i = 0; i < SIZE; i++)
	{
		delete[] result[i];
	}
	delete[] result; */

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

	result = _generateLaplacianMatrixParallelCase2(data3, SIZE);
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

	result = _generateLaplacianMatrixParallelCase3(data4, SIZE);
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

	result = _generateLaplacianMatrixParallelCase4_1(data5, SIZE);
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

	result = _generateLaplacianMatrixParallelCase4_2(data6, SIZE);
	timer.offTimer(5);
	for (int i = 0; i < SIZE; i++)
	{
		delete[] result[i];
	}
	delete[] result;

	printf("Data Size = %d\n", SIZE);
	timer.printTimer();


	EXIT_WIHT_KEYPRESS;
}

double** _generateLaplacianMatrix(double *affinityMatrix[], int n) {
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

double** _generateLaplacianMatrixParallelCase4_1(double* affinityMatrix[], int n) {
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
						if(i != j) sum += (affinityMatrix)[i][j];
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

double** _generateLaplacianMatrixParallelCase4_2(double* affinityMatrix[], int n) {
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
					if(i != j) sum += (affinityMatrix)[i][j];
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

double** _generateLaplacianMatrixParallelCase3(double* affinityMatrix[], int n) {
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

double** _generateLaplacianMatrixParallelCase1(double* affinityMatrix[], int n) {
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

double** _generateLaplacianMatrixParallelCase2(double* affinityMatrix[], int n) {
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