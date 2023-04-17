#include <vector>
#include <cmath>
#include "SpectralClustring.h"
#include "DS_definitions.h"
#include "DS_timer.h"
#include <omp.h>

#define THREADS (8)

void generateLaplacianMatrix(std::vector<std::vector<double>> *affinityMatrix) {
	std::vector<std::vector<double>> d(affinityMatrix->size());
	for (int i = 0; i < (affinityMatrix->size()); i++)
	{
		double sum = 0;
		for (int j = 0; j < (affinityMatrix->size()); j++)
		{
			sum += (*affinityMatrix)[i][j];
			if (i != j) (*affinityMatrix)[i][j] = -(*affinityMatrix)[i][j];
		}
		(*affinityMatrix)[i][i] = sum - (*affinityMatrix)[i][i];
	}
}

void generateLaplacianMatrixParallel(std::vector<std::vector<double>>* affinityMatrix) {
	std::vector<std::vector<double>> d(affinityMatrix->size());
	for (int i = 0; i < (affinityMatrix->size()); i++)
	{
		double sum = 0;
		#pragma omp parallel num_threads(THREADS)
		{
			#pragma omp for nowait
			for (int j = 0; j < (affinityMatrix->size()); j++)
			{
				if (i != j) (*affinityMatrix)[i][j] = -(*affinityMatrix)[i][j];
			}

			#pragma omp for reduction(+:sum)
			for (int j = 0; j < (affinityMatrix->size()); j++)
			{
				sum += (*affinityMatrix)[i][j];
			}
			(*affinityMatrix)[i][i] = sum - (*affinityMatrix)[i][i];
		}
	}
}