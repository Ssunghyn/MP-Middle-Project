#include <vector>
#include <cmath>
#include "SpectralClustring.h"
#include "DS_definitions.h"
#include "DS_timer.h"
#include <omp.h>

#define THREADS (8)

void generateLaplacianMatrix(std::vector<std::vector<double>> *affinityMatrix) {
	std::vector<std::vector<double>> d(affinityMatrix -> size());
	for (int i = 0; i < (affinityMatrix -> size()); i++)
	{
		double sum = 0;
		for (int j = 0; j < (affinityMatrix->size()); j++)
		{
			if (i == j)
			{
				double sum = 0;
				for (int k = 0; k < (affinityMatrix->size()); k++)
				{
					sum += (*affinityMatrix)[i][k];
				}
				(*affinityMatrix)[i][j] = sum - (*affinityMatrix)[i][j];
			}
			else (*affinityMatrix)[i][j] = - (*affinityMatrix)[i][j];
			
		}
	}
}

void generateLaplacianMatrixParallel(std::vector<std::vector<double>>* affinityMatrix) {
	std::vector<std::vector<double>> d(affinityMatrix->size());
	#pragma omp parallel num_threads(THREADS) 
	{
		#pragma omp for schedule(static, affinityMatrix->size() / THREADS)
		for (int i = 0; i < (affinityMatrix->size()); i++)
		{
			double sum = 0;
			for (int j = 0; j < (affinityMatrix->size()); j++)
			{
				if (i == j)
				{
					double sum = 0;
					for (int k = 0; k < (affinityMatrix->size()); k++)
					{
						sum += (*affinityMatrix)[i][k];
					}
					(*affinityMatrix)[i][j] = sum - (*affinityMatrix)[i][j];
				}
				else (*affinityMatrix)[i][j] = -(*affinityMatrix)[i][j];

			}
		}
	}
}