#include <stdio.h>
#include <omp.h>
#include "eigen-3.4.0/Eigen/Dense"
#include "AffinityMatrix.h"
#include "LaplacianMatrix.h"
#include "DataProcessing.h"
#include "DS_definitions.h"
#include "DS_timer.h"
#include <utility>
#include <vector>

using namespace Eigen;

#define TOL 1e-6

double** SingleMethod(Point* points,int n) {
	double** AffinityMatrix = generateAffinityMatrix(points, n);
	return generateLaplacianMatrix(AffinityMatrix, n);
}
double** MultiMethod1(Point* points, int n) {
	double** AffinityMatrix = generateAffinityMatrix_parallel(points, n);
	return generateLaplacianMatrixParallelCase1(AffinityMatrix, n);
}
/*double** MultiMethod2(Point* points, int n) {
	double** AffinityMatrix = generateAffinityMatrix_parallel(points, n);
	return generateLaplacianMatrixParallelCase2(AffinityMatrix, n);
}
double** MultiMethod3(Point* points, int n) {
	double** AffinityMatrix = generateAffinityMatrix_parallel(points, n);
	return generateLaplacianMatrixParallelCase3(AffinityMatrix, n);
}
double** MultiMethod4_1(Point* points, int n) {
	double** AffinityMatrix = generateAffinityMatrix_parallel(points, n);
	return generateLaplacianMatrixParallelCase4_1(AffinityMatrix, n);
}
double** MultiMethod4_2(Point* points, int n) {
	double** AffinityMatrix = generateAffinityMatrix_parallel(points, n);
	return generateLaplacianMatrixParallelCase4_2(AffinityMatrix, n);
}*/
void swap(double& x1, double& x2) {
	double temp = x1;
	x1 = x2;
	x2 = temp;
}

void labelDecomposition(MatrixXd& eigenVectors, int* results, int n) {
	std::vector<std::pair<double, int>> dictionary;
	int size = 0;
	for (int i = 0; i < n; i++) {
		if (dictionary.empty()) {
			std::pair<double, int> item(eigenVectors(i, 1), 1);
			dictionary.push_back(item);
			size++;
		}
		else {
			for (int j = 0; j < size; j++) {
				if (abs(dictionary[j].first - eigenVectors(i, 1)) < 1e-3) {
					dictionary[j].second += 1;
					break;
				}
				if (j == size - 1) {
					std::pair<double, int> item(eigenVectors(i, 1), 1);
					dictionary.push_back(item);
					size++;
				}
			}
		}
	}

	for (int i = 0; i < 2; i++) {
		int max = i;
		for (int j = i + 1; j < size; j++) {
			if (dictionary[max].second < dictionary[j].second) {
				max = j;
			}
		}
		std::swap(dictionary[max], dictionary[i]);
	}

	double labels[2];
	labels[0] = dictionary[0].first;
	labels[1] = dictionary[1].first;

	//printf("Size = %d\n", size);
	//printf("label 1 %lf\nlabel 2 %lf\n", labels[0], labels[1]);

	for (int i = 0; i < n; i++) {
		if (abs(eigenVectors(i, 1) - labels[0]) > abs(eigenVectors(i, 1) - labels[1])) {
			results[i] = 1;
		}
		else {
			results[i] = -1;
		}
	}
}

void find_2nd_Min(VectorXd& eigenValue, MatrixXd& eigenVector) {
	int n = eigenVector.rows();

	for (int i = 0; i < 2; i++) {
		int min = i;
		for (int j = 1; j < n; j++) {
			if (eigenValue[min] > eigenValue[j]) {
				min = j;
			}
		}
		swap(eigenValue[min], eigenValue[i]);
		for (int j = 0; j < n; j++) {
			swap(eigenVector(j,min), eigenVector(j,i));
		}
	}
}

//int realMain(int argc, char** argv)
int main(int argc, char** argv)
{
	std::string fname = argv[1];
	std::string file_name = ".\\data\\"+ fname;
	int n = atoi(argv[2]);
	FILE* fp;
    #ifdef _WIN64
	fopen_s(&fp, file_name.c_str(), "r");
    #else
	fp = fopen(file_name.c_str(), "r");
    #endif
	Point* points = new Point[n];
	bool isOpen = getData(fp, points, n);
	if (!isOpen) {
		return -1;
	}

	std::string name[6] = { "Single Method", "Multi Method1","Multi Method2" ,"Multi Method3" ,"Multi Method4_1" ,"Multi Method4_2" };
	DS_timer timer(6);
	for (int i = 0; i < 6; i++)
		timer.setTimerName(i, name[i]);

	timer.onTimer(0);
	double** dataSingle = SingleMethod(points, n);
	timer.offTimer(0);

	timer.onTimer(1);
	double** dataMuliti1 = MultiMethod1(points, n);
	timer.offTimer(1);

	/*
	timer.onTimer(2);
	double** dataMuliti2 = MultiMethod2(points, n);
	timer.offTimer(2);

	timer.onTimer(3);
	double** dataMuliti3 = MultiMethod3(points, n);
	timer.offTimer(3);

	timer.onTimer(4);
	double** dataMuliti4_1 = MultiMethod4_1(points, n);
	timer.offTimer(4);

	timer.onTimer(5);
	double** dataMuliti4_2 = MultiMethod4_2(points, n);
	timer.offTimer(5);
	*/

	bool isCorrect = true;
	int idx[2] = { 0 };
	/*
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double single = dataSingle[i][j];
			if (abs(single - dataMuliti1[i][j]) > TOL || abs(single - dataMuliti2[i][j]) > TOL || abs(single - dataMuliti3[i][j]) > TOL
				|| abs(single - dataMuliti4_1[i][j]) > TOL || abs(single - dataMuliti4_2[i][j]) > TOL) {
				isCorrect = false;
				idx[0] = i; idx[1] = j;
				break;
			}
		}
	}
	*/

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double single = dataSingle[i][j];
			if (abs(single - dataMuliti1[i][j]) > TOL) {
				isCorrect = false;
				idx[0] = i; idx[1] = j;
				break;
			}
		}
	}

	if (isCorrect) {
		printf("Data is correct.\n");
	}
	else {
		printf("Data is not correct. ");
		printf("Single[%d][%d] : %lf\n", idx[0], idx[1],dataSingle[idx[0]][idx[1]]);
		/*
		printf("dataMuliti1[%d][%d] : %lf\n", idx[0], idx[1], dataMuliti1[idx[0]][idx[1]]);
		printf("dataMuliti2[%d][%d] : %lf\n", idx[0], idx[1], dataMuliti2[idx[0]][idx[1]]);
		printf("dataMuliti3[%d][%d] : %lf\n", idx[0], idx[1], dataMuliti3[idx[0]][idx[1]]);
		printf("dataMuliti4_1[%d][%d] : %lf\n", idx[0], idx[1], dataMuliti4_1[idx[0]][idx[1]]);
		printf("dataMuliti4_2[%d][%d] : %lf\n", idx[0], idx[1], dataMuliti4_2[idx[0]][idx[1]]);
		*/
	}

	MatrixXd A(n, n);
	MatrixXd B(n, n);

	for (int i = 0; i < 8; i++)
	{
		file_name.pop_back();
	}

	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A(i, j) = dataSingle[i][j]; 
			B(i, j) = dataMuliti1[i][j];
		}
	}

	std::string version[] = { "SingleMethod", "MultiMethod1", "MultiMethod2", "MultiMethod3", "MultiMethod4_1", "MultiMethod4_2" };

	#pragma omp parallel for num_threads(2)
	for (int tid = 0; tid < 2; tid++) {
		int* results = new int[n];

		MatrixXd eigenVectors(n, n);
		VectorXd eigenValues(n);
		EigenSolver<MatrixXd> solver;
		if (tid == 0)
			solver = EigenSolver<MatrixXd>(A);
		else
			solver = EigenSolver<MatrixXd>(B);

		MatrixXcd solveEigenVector = solver.eigenvectors();
		VectorXcd solveEigenValue = solver.eigenvalues();

		for (int i = 0; i < n; i++) {
			eigenValues[i] = solveEigenValue[i].real();
			for (int j = 0; j < n; j++) {
				eigenVectors(i, j) = solveEigenVector(i, j).real();
			}
		}

		find_2nd_Min(eigenValues, eigenVectors);
		labelDecomposition(eigenVectors, results, n);
		
		std::string saveName;
		saveName = file_name + version[tid] + "_result.txt";
		saveData(saveName.c_str(), points, results, n);

		delete[] results;
	}
	
	timer.printTimer();

	/*
	for (int i = 0; i < n; i++) {
		delete[] dataSingle[i], dataMuliti1[i], dataMuliti2[i], dataMuliti3[i], dataMuliti4_1[i], dataMuliti4_2[i];
	}
	*/

	for (int i = 0; i < n; i++) {
		delete[] dataSingle[i], dataMuliti1[i];
	}
	delete[] points;

	return 0;
}