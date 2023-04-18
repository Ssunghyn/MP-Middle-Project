#include <vector>
#include <cmath>
#include "SpectralClustring.h"
#include "DS_definitions.h"
#include "DS_timer.h"
#include <omp.h>

#define THREADS (8)

void generateLaplacianMatrix(double* affinityMatrix[], int n) {
	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int j = 0; j < n; j++)
		{
			if (i != j)
			{
				sum += (affinityMatrix)[i][j];
				(affinityMatrix)[i][j] = -(affinityMatrix)[i][j];
			}
		}
		(affinityMatrix)[i][i] = sum;
	}
}

void generateLaplacianMatrixParallel2(double* affinityMatrix[], int n) {
	#pragma omp parallel num_threads(THREADS) 
	{
		#pragma omp for
		for (int i = 0; i < n; i++)
		{
			double sum = 0;
			for (int j = 0; j < n; j++)
			{
				if (i != j)
				{
					sum += (affinityMatrix)[i][j];
					(affinityMatrix)[i][j] = -(affinityMatrix)[i][j];
				}
			}
			(affinityMatrix)[i][i] = sum;
		}
	}
}