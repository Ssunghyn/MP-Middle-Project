#include <cmath>
#include "DS_definitions.h"
#include "DS_timer.h"
#include "AffinityMatrix.h"
#include <omp.h>

#define SIGMA_DIMENSION 5
#define TEST_DATA_COUNT 10000
#define OMP_OFFSET 16
#define NUM_THREADS 8
#define PRINT_RESULT false
#define USE_V2 false

#define GenDouble (rand() % 4 + ((double)(rand() % 100) / 100.0))

inline double getDistance(Point p1, Point p2) {
    double x = pow(p2.x - p1.x, 2);
    double y = pow(p2.y - p1.y, 2);
    return sqrt(x + y);
}

double quickSelection(double* distances, double* pivot_left, double* pivot_right, int size, int k) {
    double pivot = distances[0];
    int same = 0;
    int pivot_left_count = 0;
    int pivot_right_count = 0;

    for(int i = 0; i < size; i++) {
        if(distances[i] < pivot) {
            pivot_left[pivot_left_count++] = distances[i];
        } else if(distances[i] > pivot) {
            pivot_right[pivot_right_count++] = distances[i];
        } else {
            same++;
        }
    }

    if(k <= pivot_left_count) {
        return quickSelection(pivot_left, pivot_left, pivot_right, pivot_left_count, k);
    } else if(k >= pivot_left_count + 1 && k <= pivot_left_count + same) {
        return pivot; 
    } else {
        return quickSelection(pivot_right, pivot_left, pivot_right, pivot_right_count, k - pivot_left_count - same);
    }
}

double** generateAffinityMatrix(Point* points, int point_count) {
    double** result = new double*[point_count];
    double** distance = new double*[point_count];
    double* deltas = new double[point_count];

    double* pivot_left = new double[point_count];
    double* pivot_right = new double[point_count];

    // Get distance of 2 points
    for(int p1 = 0; p1 < point_count; p1++) {
        double* partial_distance = new double[point_count];

        for(int p2 = 0; p2 < point_count; p2++) {
            #if USE_V2
            if(p1 > p2) {
                partial_distance[p2] = distance[p2][p1];
            } else {
                partial_distance[p2] = getDistance(points[p1], points[p2]);
            }
            #else
            partial_distance[p2] = getDistance(points[p1], points[p2]);
            #endif
        }

        distance[p1] = partial_distance;
        deltas[p1] = quickSelection(partial_distance, pivot_left, pivot_right, point_count, SIGMA_DIMENSION);
    }

    // Make affinity matrix
    for(int p1 = 0; p1 < point_count; p1++) {
        double* partial_result = new double[point_count];

        for(int p2 = 0; p2 < point_count; p2++) {
            double affinity;

            #if USE_V2 
            if(p1 > p2) {
                affinity = result[p2][p1];
            } else {
                affinity = exp(-distance[p1][p2] / (2 * deltas[p1] * deltas[p2]));
            }
            #else
                affinity = exp(-distance[p1][p2] / (2 * deltas[p1] * deltas[p2]));
            #endif
            partial_result[p2] = affinity;
        }

        result[p1] = partial_result;
    }
    
    // Free memory
    for(int i = 0; i < point_count; i++) {
        delete[] distance[i];
    }

    delete[] distance;
    delete[] deltas;
    delete[] pivot_left;
    delete[] pivot_right;

    return result;
}

double** generateAffinityMatrix_parallel(Point* points, int point_count) {
    double** result = new double*[point_count];
    double** distance = new double*[point_count];
    double* deltas = new double[point_count];

    #pragma omp parallel num_threads(NUM_THREADS) shared(result, distance, deltas)
    {
        double* local_pivot_left = new double[point_count];
        double* local_pivot_right = new double[point_count];
        // Get distance of 2 points
        #pragma omp for
        for(int p1 = 0; p1 < point_count; p1++) {
            int tId = omp_get_thread_num();
            double* partial_distance = new double[point_count];

            #if false
            for(int p2 = p1; p2 < point_count; p2++) {
            #else
            for(int p2 = 0; p2 < point_count; p2++) {
            #endif
                partial_distance[p2] = getDistance(points[p1], points[p2]);
            }

            distance[p1] = partial_distance;
            #if !USE_V2
            deltas[p1] = quickSelection(partial_distance, local_pivot_left, local_pivot_right, point_count, SIGMA_DIMENSION);
            #endif
        }

        #if USE_V2
        #pragma omp for
        for(int p1 = 0; p1 < point_count; p1++) {
            for(int p2 = 0; p2 < p1; p2++) {
                distance[p1][p2] = distance[p2][p1];
            }

            deltas[p1] = quickSelection(distance[p1], local_pivot_left, local_pivot_right, point_count, SIGMA_DIMENSION);
        }
        #endif

        // Make affinity matrix
        #pragma omp for
        for(int p1 = 0; p1 < point_count; p1++) {
            double* partial_result = new double[point_count];

            #if USE_V2
            for(int p2 = p1; p2 < point_count; p2++) {
            #else
            for(int p2 = 0; p2 < point_count; p2++) {
            #endif
                partial_result[p2] = exp(-distance[p1][p2] / (2 * deltas[p1] * deltas[p2]));
            }

            result[p1] = partial_result;
        }

        #if USE_V2
        #pragma omp for
        for(int p1 = 0; p1 < point_count; p1++) {
            for(int p2 = 0; p2 < p1; p2++) {
                result[p1][p2] = result[p2][p1];
            }
        }
        #endif

        delete[] local_pivot_left;
        delete[] local_pivot_right;
    }
    
    // Free memory
    for(int i = 0; i < point_count; i++) {
        delete[] distance[i];
    }

    delete[] distance;
    delete[] deltas;

    return result;
}

void AfffinityMain()
//int main() 
{
    Point points[TEST_DATA_COUNT];

    printf("Data size: %d\n", TEST_DATA_COUNT);

    DS_timer timer(2);
	timer.setTimerName(0, (char*)"Serial");
	timer.setTimerName(1, (char*)"Parallel");
    for(int i = 0; i < TEST_DATA_COUNT; i++) {
        points[i] = Point(GenDouble, GenDouble);
    }

    timer.onTimer(0);
    double** result = generateAffinityMatrix(points, TEST_DATA_COUNT);
    timer.offTimer(0);

    timer.onTimer(1);
    double** result2 = generateAffinityMatrix_parallel(points, TEST_DATA_COUNT);
    timer.offTimer(1);


    bool same = true;

    #if PRINT_RESULT
        for(int i = 0; i < TEST_DATA_COUNT; i++) {
            printf("(%lf, %lf)\n", points[i].x, points[i].y);
        }

        for(int i = 0; i < TEST_DATA_COUNT; i++) {
            for(int j = 0; j < TEST_DATA_COUNT; j++) {
                printf("%.2lf ", result[i][j]);
            }
            printf("\n");
        }

        printf("---------------------------------------------------------------------\n");

        for(int i = 0; i < TEST_DATA_COUNT; i++) {
            for(int j = 0; j < TEST_DATA_COUNT; j++) {
                printf("%lf ", result2[i][j]);
            }
            printf("\n");
        }

        printf("---------------------------------------------------------------------\n");

        for(int i = 0; i < TEST_DATA_COUNT; i++) {
            for(int j = 0; j < TEST_DATA_COUNT; j++) {
                printf("%lf ", result[i][j] - result2[i][j]);
                if(abs(result[i][j] - result2[i][j]) > 0.0000000001) {
                    same = false;
                }
            }
            printf("\n");
        }
    #endif

    for(int i = 0; i < TEST_DATA_COUNT; i++) {
        for(int j = 0; j < TEST_DATA_COUNT; j++) {
            if(abs(result[i][j] - result2[i][j]) > 0.0000000001) {
                same = false;
            }
        }
    }

    if(same) {
        printf("RESULT IS SAME");
    } else {
        printf("RESULT IS NOT SAME");
    }

    timer.printTimer();
}
