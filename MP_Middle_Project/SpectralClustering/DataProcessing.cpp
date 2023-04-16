#include "DataProcessing.h"

void getData(FILE* fp, std::vector<Point>& points) {
	if (fp == NULL) {
		printf("FILE NOT OPEN!\n");
		return;
	}

	while (!feof(fp)) {
		double x = 0;
		double y = 0;
		
		fscanf_s(fp, "%lf %lf", &x, &y);
		points.push_back(Point(x, y));
	}
}


#include <iostream>
#include <cmath>

using namespace std;


void qr_algorithm(double** A, int n, double** Q, double** R) {
    // Initialize Q to the identity matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                Q[i][j] = 1.0;
            }
            else {
                Q[i][j] = 0.0;
            }
        }
    }

    double epsilon = 1e-6;  // stopping criterion
    int max_iterations = 100;  // maximum number of iterations

    for (int k = 0; k < max_iterations; k++) {
        // Compute the QR decomposition of A
        for (int j = 0; j < n; j++) {
            // Compute the j-th column of Q and the j-th row of R
            double* qj = new double[n];
            double* rj = new double[n];
            for (int i = 0; i < n; i++) {
                qj[i] = A[i][j];
            }
            for (int i = 0; i < j; i++) {
                double dot_product = 0.0;
                for (int l = 0; l < n; l++) {
                    dot_product += Q[l][i] * A[l][j];
                }
                rj[i] = dot_product;
                for (int l = 0; l < n; l++) {
                    qj[l] -= Q[l][i] * rj[i];
                }
            }
            double norm_qj = 0.0;
            for (int i = 0; i < n; i++) {
                norm_qj += qj[i] * qj[i];
            }
            norm_qj = sqrt(norm_qj);
            for (int i = 0; i < n; i++) {
                Q[i][j] = qj[i] / norm_qj;
            }
            rj[j] = norm_qj;
            for (int i = j; i < n; i++) {
                R[j][i] = 0.0;
            }
            for (int i = j + 1; i < n; i++) {
                double dot_product = 0.0;
                for (int l = 0; l < n; l++) {
                    dot_product += Q[l][j] * A[l][i];
                }
                R[j][i] = dot_product;
                for (int l = 0; l < n; l++) {
                    A[l][i] -= Q[l][j] * R[j][i];
                }
            }

            delete[] qj;
            delete[] rj;
        }

        // Compute the norm of the subdiagonal elements of R
        double norm_subdiag = 0.0;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                norm_subdiag += R[j][i] * R[j][i];
            }
        }
        norm_subdiag = sqrt(norm_subdiag);

        // Check if the norm of the subdiagonal elements of R is small enough
        if (norm_subdiag < epsilon) {
            break;
        }
    }
}

int main() {
    int n = 3;

    double** A = new double* [n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
    }

    A[0][0] = 1.0; A[0][1] = 2.0; A[0][2] = 3.0;
    A[1][0] = 4.0; A[1][1] = 5.0; A[1][2] = 6.0;
    A[2][0] = 7.0; A[2][1] = 8.0; A[2][2] = 9.0;

    double** Q = new double* [n];
    for (int i = 0; i < n; i++) {
        Q[i] = new double[n];
        for (int j = 0; j < n; j++) {
            Q[i][j] = 0.0;
        }
    }

    double** R = new double* [n];
    for (int i = 0; i < n; i++) {
        R[i] = new double[n];
        for (int j = 0; j < n; j++) {
            R[i][j] = 0.0;
        }
    }

    qr_algorithm(Q, n, Q, R);
     
    cout << "Q mat" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%02.6lf ", Q[i][j]);
        }
        cout << endl;
    }
    cout << endl;
    cout << "R mat" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%02.6lf ", R[i][j]);
        }
        cout << endl;
    }
}


int data_main()
//int main() 
{
	input_path = "circle_data.txt";

	fopen_s(&fp, input_path.c_str(), "r");
	std::vector<Point> points;
	getData(fp, points);
	for (int i = 0; i < points.size(); i++)
		printf("%lf %lf\n", points[i].x, points[i].y);
}