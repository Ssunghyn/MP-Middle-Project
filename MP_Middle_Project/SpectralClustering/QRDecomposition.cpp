#include "QRDecomposition.h"
#include "DS_definitions.h"
#include "DS_timer.h"
#include "jacobi_pd.hpp"
#include "matrix_alloc_jpd.hpp"

#define GenFloat rand()%100 + ((rand()%100)/100.0)

void printMat(double** A, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			/*
			if (A[i][j] < 0)
				printf("%9.6lf   ", A[i][j]);
			else
				printf("%10.6lf   ", A[i][j]);
			*/
			printf("%11.6lf   ", A[i][j]);
		}
		printf("\n");
	}
}
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
void QR_Factorization(double** A, double** Q, double** R, double** H, int n) {

	// step 1. set vector c, e, v
	double* c = new double[n];
	double e = 0;
	double* v = new double[n];
	for (int i = 0; i < n; i++)
		c[i] = A[i][0];

	e = c[0] > 0 ? 1 : -1;

	double normC = norm(c, n);
	v[0] = c[0] + normC * e;
	for (int i = 1; i < n; i++) {
		v[i] = c[i];
	}

	//double factor = dots(v, v, n);
	double factor = 0;
	for (int i = 0; i < n; i++)
		factor += v[i] * v[i];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			H[i][j] = -2.0 / factor * v[i] * v[j];
		}
		H[i][i] += 1;
	}

	// Q = H, R = HA
	//copy(H, Q, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Q[i][j] = H[i][j];
			R[i][j] = 0;
		}
	}
	//matMul(H, A, R, n);
	for (int i = 0; i < n; i++) {
		double sum = 0;
		for (int k = 0; k < n; k++) {
			for (int j = 0; j < n; j++) {
				R[i][j] += H[i][k] * A[k][j];
			}
		}
	}
	//double** RT = Transform(R, n);  // R.T

	for (int i = 1; i < n - 1; i++) {

		for (int j = 0; j < n; j++) {
			if (j < i)
				c[j] = 0;
			else
				c[j] = R[j][i];
		}

		e = c[i] > 0 ? 1 : -1;

		normC = norm(c, n);
		for (int j = 0; j < n; j++) {
			v[j] = c[j];
		}
		v[i] += normC * e;

		factor = dots(v, v, n);
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				H[j][k] = -2.0 / factor * v[j] * v[k];
			}
			H[j][j] += 1;
		}

		matrixMul(Q, H, Q, n);
		matrixMul(H, R, R, n);

	}
	//R = Transform(RT, n);

	delete[] c, v;
}
bool isDiag(double** A, int n, double tol) {
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (abs(A[i][j]) > tol) {
				return false;
			}
		}
	}
	return true;
}
void findEigenQR(double** A, double* eigenValue, double** eigenVector, int maxIter, int n, double tol)
{
	double** A_new = new double*[n];
	double** Q = new double*[n];
	double** R = new double*[n];
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
		}
	}
	
	int count = 0;

	while (!isDiag(A_new, n, tol) && (count < maxIter))
	{
		QR_Factorization(A_new, Q, R, H, n);
		matrixMul(eigenVector, Q, eigenVector, n);
		matrixMul(R, Q, A_new, n);
		count++;
	}

	for (int i = 0; i < n; i++) {
		eigenValue[i] = A_new[i][i];
	}

	for (int i = 0; i < n; i++) {
		delete[] A_new[i], Q[i], R[i], H[i];
	}
}

void qrDecompositionMain()
//void main()
{
	DS_timer timer(2);
	std::string name[2] = { "QR Decomposition", "Jacobi Method" };
	timer.setTimerName(0, name[0]);
	timer.setTimerName(1, name[1]);
	int size = 10;

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

	for (int i = 0; i < size; i++) {
		double row_sum = 0;
		for(int j = 0; j < size; j++)
			row_sum += A[i][j];
		A[i][i] = row_sum - A[i][i];
		for (int j = 0; j < size; j++) {
			if (i != j) {
				A[i][j] *= -1;
			}
			//B(i, j) = A[i][j];
		}
	}


	printf("[");
	for (int i = 0; i < size - 1; i++) {
		printf("[%lf, ", A[i][0]);
		for (int j = 1; j < size - 1; j++) {
			printf("%lf, ", A[i][j]);
		}
		printf("%lf],\n", A[i][size - 1]);
	}
	printf("[");
	for (int j = 0; j < size - 1; j++)
		printf("%lf, ", A[size - 1][j]);
	printf("%lf]]\n\n", A[size - 1][size - 1]);

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
		/*
		if (eigValQR[i] * eigValJacobi[i] < 0) {
			printf("[%d]  %11.6lf  %11.6lf\n", i, eigValQR[i], eigValJacobi[i]);
		}
		*/
		printf("[%d]  %11.6lf  %11.6lf\n", i, eigValQR[i], eigValJacobi[i]);
	}
	
	timer.printTimer();

	for (int i = 0; i < size; i++) {
		delete[] A[i], eigVecQR[i], eigVecJacobi[i];
	}
	delete[] eigValQR, eigValJacobi;
}