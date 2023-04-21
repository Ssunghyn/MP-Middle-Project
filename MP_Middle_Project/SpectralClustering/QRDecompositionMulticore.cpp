#include "QRDecomposition.h"
#include "QRDecompositionMulticore.h"
#include "DS_definitions.h"
#include "DS_timer.h"
#include "jacobi_pd.hpp"
#include "matrix_alloc_jpd.hpp"
#include "omp.h"
#include "SpectralClustring.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#define GenFloat rand()%100 + ((rand()%100)/100.0)

void MultiTransform(double** A, double** AT, int n)
{
	#pragma omp parallel for num_threads(THREADS)
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			AT[i][j] = A[j][i];
}
void Multicopy(double** input, double** output, int n)
{
	#pragma omp parallel for num_threads(THREADS)
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			output[i][j] = input[i][j];
}
double Multinorm(double* v, int n)
{
	double result = 0.0;
	#pragma omp parallel for num_threads(THREADS) reduction(+: result)
	for (int i = 0; i < n; i++) {
		result += pow(v[i], 2);
	}
	return sqrt(result);
}
void MultimatrixMul(double** A, double** B, double** output, int n)
{
	#pragma omp parallel for num_threads(THREADS)
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double sum = 0.0;
			for (int k = 0; k < n; k++) {
				sum += A[i][k] * B[k][j];
			}
			output[i][j] = sum;
		}
	}
}
double Multidots(double* A, double* B, int n)
{
	double C = 0.0;
	#pragma omp parallel for num_threads(THREADS)
	for (int i = 0; i < n; i++) {
		C += A[i] * B[i];
	}
	return C;
}
void MultiQR_Factorization(double** A, double** Q, double** R, double** H, int n) {

	// step 1. set vector c, e, v
	double* c = new double[n];
	double e = 0;
	double* v = new double[n];
	for (int i = 0; i < n; i++)
		c[i] = A[i][0];

	e = c[0] > 0 ? 1 : -1;

	double normC = Multinorm(c, n);
	v[0] = c[0] + normC * e;
	for (int i = 1; i < n; i++) {
		v[i] = c[i];
	}

	double factor = Multidots(v, v, n);
	/*
	double factor = 0;
	for (int i = 0; i < n; i++)
		factor += v[i] * v[i];
	*/
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			H[i][j] = -2.0 / factor * v[i] * v[j];
		}
		H[i][i] += 1;
	}

	// Q = H, R = HA
	Multicopy(H, Q, n);
	/*
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Q[i][j] = H[i][j];
			R[i][j] = 0;
		}
	}
	*/
	MultimatrixMul(H, A, R, n);
	/*
	for (int i = 0; i < n; i++) {
		double sum = 0;
		for (int k = 0; k < n; k++) {
			for (int j = 0; j < n; j++) {
				R[i][j] += H[i][k] * A[k][j];
			}
		}
	}
	*/
	//double** RT = Transform(R, n);  // R.T

	for (int i = 1; i < n - 1; i++) {

		for (int j = 0; j < n; j++) {
			if (j < i)
				c[j] = 0;
			else
				c[j] = R[j][i];
		}

		e = c[i] > 0 ? 1 : -1;

		normC = Multinorm(c, n);
		for (int j = 0; j < n; j++) {
			v[j] = c[j];
		}
		v[i] += normC * e;

		factor = Multidots(v, v, n);
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				H[j][k] = -2.0 / factor * v[j] * v[k];
			}
			H[j][j] += 1;
		}

		MultimatrixMul(Q, H, Q, n);
		MultimatrixMul(H, R, R, n);

	}
	//R = Transform(RT, n);

	delete[] c, v;
}
void MultifindEigenQR(double** A, double* eigenValue, double** eigenVector, int maxIter, int n, double tol)
{
	double** A_new = new double* [n];
	double** Q = new double* [n];
	double** R = new double* [n];
	double** H = new double* [n];

	for (int i = 0; i < n; i++) {
		eigenVector[i][i] = 1.0;
		H[i] = new double[n];
		A_new[i] = new double[n];
		R[i] = new double[n];
		Q[i] = new double[n];
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A_new[i][j] = A[i][j];
			R[i][j] = 0;
		}
	}

	int count = 0;

	while (!isDiag(A_new, n, tol) && (count < maxIter))
	{
		MultiQR_Factorization(A_new, Q, R, H, n);
		MultimatrixMul(eigenVector, Q, eigenVector, n);
		MultimatrixMul(R, Q, A_new, n);
		if (count == 0)
			printMat(A_new, n);
		count++;
	}

	

	for (int i = 0; i < n; i++) {
		eigenValue[i] = A_new[i][i];
	}

	for (int i = 0; i < n; i++) {
		delete[] A_new[i], Q[i], R[i], H[i];
	}
}

using namespace Eigen;
void qrDecompositionMain()
//void main()
{
	DS_timer timer(2);
	std::string name[2] = { "QR Decomposition", "Jacobi Method" };
	timer.setTimerName(0, name[0]);
	timer.setTimerName(1, name[1]);
	int size = 1000;
	MatrixXd B(size, size);
	double** A = new double* [size];
	for (int i = 0; i < size; i++)
		A[i] = new double[size];

	double** eigVecQR = new double* [size];
	double** eigVecJacobi = new double* [size];
	for (int i = 0; i < size; i++) {
		eigVecQR[i] = new double[size];
		eigVecJacobi[i] = new double[size];
	}

	double* eigValQR = new double[size];
	double* eigValJacobi = new double[size];
	

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i <= j)
				A[i][j] = GenFloat;
			else
				A[i][j] = A[j][i];
		}
	}

	for (int i = 0; i < size; i++) {
		double row_sum = 0;
		for (int j = 0; j < size; j++)
			row_sum += A[i][j];
		A[i][i] = row_sum - A[i][i];
		for (int j = 0; j < size; j++) {
			if (i != j) {
				A[i][j] *= -1;
			}
			//B(i, j) = A[i][j];
		}
	}

	timer.setTimerName(0, name[0]);
	timer.onTimer(0);
	// 고유벡터를 저장할 행렬 생성
	MatrixXd eigenvectors(B.rows(), B.cols());

	// 행렬을 LU 분해하여 저장
	PartialPivLU<MatrixXd> lu(B);

	// 고윳값 계산
	VectorXd eigenvalues = lu.matrixLU().diagonal();

	// 고유벡터 계산
	for (int i = 0; i < B.cols(); i++) {
		VectorXd b = VectorXd::Zero(B.cols());
		b(i) = 1.0;
		eigenvectors.col(i) = lu.solve(b);
	}
	timer.offTimer(0);


	
	initParallel();
	timer.onTimer(1);
	MatrixXd eigenvectorsMul(B.rows(), B.cols());

	// 행렬을 LU 분해하여 저장
	PartialPivLU<MatrixXd> luMul(B);

	// 고윳값 계산
	VectorXd eigenvaluesMul = luMul.matrixLU().diagonal();

	// 고유벡터 계산
	for (int i = 0; i < B.cols(); i++) {
		VectorXd b = VectorXd::Zero(B.cols());
		b(i) = 1.0;
		eigenvaluesMul.col(i) = luMul.solve(b);
	}
	timer.offTimer(1);
	

	/*
	sort(eigValQR, size);
	sort(eigValJacobi, size);



	for (int i = 0; i < size; i++) {
		
		if (eigValQR[i] * eigValJacobi[i] < 0) {
			printf("[%d]  %11.6lf  %11.6lf\n", i, eigValQR[i], eigValJacobi[i]);
		}
		
		printf("[%d]  %11.6lf  %11.6lf\n", i, eigValQR[i], eigValJacobi[i]);
	}
	*/
	timer.printTimer();

	for (int i = 0; i < size; i++) {
		delete[] A[i], eigVecQR[i], eigVecJacobi[i];
	}
	delete[] eigValQR, eigValJacobi;
}