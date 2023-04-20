#include "QRDecomposition.h"
#include "DS_definitions.h"
#include "DS_timer.h"
#include "jacobi_pd.hpp"
#include "matrix_alloc_jpd.hpp"

#define GenFloat rand()%100 + ((rand()%100)/100.0)

void Transform(double** A, double** AT, int n) 
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			AT[i][j] = A[j][i];
}
void copy(double** input, double** output, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			output[i][j] = input[i][j];
}
void swap(double& x1, double& x2)
{
	double temp = x1;
	x1 = x2;
	x2 = temp;
}
void sort(double* A, int n)
{
	for (int i = 0; i < n - 1; i++) {
		int min = i;
		for (int j = i + 1; j < n; j++) {
			if (A[min] > A[j])
				min = j;
		}
		swap(A[i], A[min]);
	}
}
double norm(double* v, int n)
{
	double result = 0.0;
	for (int i = 0; i < n; i++) {
		result += pow(v[i], 2);
	}
	return sqrt(result);
}
void matrixMul(double** A, double** B, double** output, int n)
{
	double** C = new double*[n];
	for (int i = 0; i < n; i++)
		C[i] = new double[n];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			C[i][j] = 0;

	for (int i = 0; i < n; i++) {
		for (int k = 0; k < n; k++) {
			for (int j = 0; j < n; j++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	copy(C, output, n);

	for (int i = 0; i < n; i++)
		delete[] C[i];
}
void matrixMul(double** A, double** output, double lambda, int n)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			output[i][j] = A[i][j] * lambda;
	}
}
double dots(double* A, double* B, int n)
{
	double C = 0.0;
	for (int i = 0; i < n; i++) {
		C += A[i] * B[i];
	}
	return C;
}
void QR_Decomposition(double** A, double** Q, double** R, int n)
{
	double** uT = new double* [n];
	for (int i = 0; i < n; i++)
		uT[i] = new double[n];

	double** AT = new double*[n];
	for (int i = 0; i < n; i++)
		AT[i] = new double[n];
	Transform(A, AT, n);

	for (int i = 0; i < n; i++)
		uT[0][i] = AT[0][i];


	for (int i = 0; i < n; i++) {
		Q[0][i] = uT[0][i] / norm(uT[0], n);
	}

	for (int i = 1; i < n; i++) {
		// u_{i th} 열 벡터에 A_{i th} 열벡터 값 대입
		for (int j = 0; j < n; j++)
			uT[i][j] = AT[i][j];

		// 각 u 벡터 구하기
		for (int j = 0; j < i; j++) {
			double mulAQ = dots(AT[i], Q[j], n);
			for (int k = 0; k < n; k++) {
				uT[i][k] -= mulAQ * Q[j][k];
			}
		}

		// 각 정규화된 e 벡터 계산
		double norm_u_i = norm(uT[i], n);

		for (int j = 0; j < n; j++) {
			Q[i][j] = uT[i][j] / norm_u_i;
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			R[i][j] = dots(AT[j], Q[i], n);
		}
	}

	Transform(Q, AT, n);
	copy(AT, Q, n);

	for (int i = 0; i < n; i++)
		delete[] uT[i], AT[i];
}
double getDifference(double** A, double** B, int n)
{
	double max = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double diff = abs(A[i][j] - B[i][j]);
			if (diff > max) {
				max = diff;
			}
		}
	}
	return max;
}
void findEigenQR(double** A, double* eigenValue, double** eigenVector, int maxIter, int n, double tol)
{
	for (int i = 0; i < n; i++) {
		eigenVector[i][i] = 1.0;
	}

	double** A_old = new double*[n];
	for (int i = 0; i < n; i++)
		A_old[i] = new double[n];

	double** A_new = new double*[n];
	for (int i = 0; i < n; i++)
		A_new[i] = new double[n];

	double** Q = new double*[n];
	for (int i = 0; i < n; i++)
		Q[i] = new double[n];

	double** R = new double*[n];
	for (int i = 0; i < n; i++)
		R[i] = new double[n];

	copy(A, A_old, n);
	copy(A, A_new, n);
	int count = 0;
	double diff = INT_MAX;

	while ((diff > tol) && (count < maxIter))
		//while ((diff > tol))
		//while (count < maxIter)
	{
		copy(A_new, A_old, n);
		QR_Decomposition(A_old, Q, R, n);
		matrixMul(eigenVector, Q, eigenVector, n);
		matrixMul(R, Q, A_new, n);
		diff = getDifference(A_new, A_old, n);
		count++;
	}

	for (int i = 0; i < n; i++) {
		eigenValue[i] = A_new[i][i];
	}

	//printf("%d 번 돌았어요!\n", count);
	//printf("%lf 의 오차가 있어요\n\n", diff);

	for (int i = 0; i < n; i++) {
		delete[] A_old[i], A_new[i], Q[i], R[i];
	}
}

//void qrDecompositionMain()
void main()
{
	DS_timer timer(2);
	std::string name[2] = { "QR Decomposition", "Jacobi Method" };
	timer.setTimerName(0, name[0]);
	timer.setTimerName(1, name[1]);
	int size = 500;

	double** A = new double*[size];
	for (int i = 0; i < size; i++)
		A[i] = new double[size];

	double** eigVecQR = new double*[size];
	double** eigVecJacobi = new double*[size];
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

	timer.onTimer(0);
	findEigenQR(A, eigValQR, eigVecQR, 250, size);
	timer.offTimer(0);

	jacobi_pd::Jacobi<double, double*, double**> eigen_cal(size);
	timer.onTimer(1);
	eigen_cal.Diagonalize(A, eigValJacobi, eigVecJacobi);
	timer.offTimer(1);

	sort(eigValQR, size);
	sort(eigValJacobi, size);

	for (int i = 0; i < size; i++) {
		if (eigValQR[i] * eigValJacobi[i] < 0) {
			printf("[%d]  %11.6lf  %11.6lf\n", i, eigValQR[i], eigValJacobi[i]);
		}
	}
	
	timer.printTimer();

	for (int i = 0; i < size; i++) {
		delete[] A[i], eigVecQR[i], eigVecJacobi[i];
	}
	delete[] eigValQR, eigValJacobi;
}