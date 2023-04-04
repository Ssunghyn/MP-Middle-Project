#include <vector>
#include <cmath>
#include "SpectralClustring.h"
#include "DS_definitions.h"
#include "DS_timer.h"
#include <omp.h>

#define SIGMA_DIMENSION 5
#define TEST_DATA_COUNT 9000
#define THREADS 8
#define OMP_OFFSET 16

#define GenDouble (rand() % 4 + ((double)(rand() % 100) / 100.0))

inline double getDistance(Point p1, Point p2) {
    double x = pow(p2.x - p1.x, 2);
    double y = pow(p2.y - p1.y, 2);
    return sqrt(x + y);
}

double quickSelection(std::vector<double> others, int k) {
    double pivot = others[0];
    int same = 0;

    std::vector<double> pivot_left, pivot_right;

    /*for(double x : others) {
        printf("%3lf ", x);
    }
    printf("\npivot = %3lf | ", pivot);*/

    for(double x : others) {
        if(x < pivot) {
            pivot_left.push_back(x);
        } else if(x > pivot) {
            pivot_right.push_back(x);
        } else {
            same++;
        }
    }

    /*for(double x : pivot_left) {
        printf("%3lf ", x);
    }
    printf(" /  ");
    for(double x : pivot_right) {
        printf("%3lf ", x);
    }
    printf("\n");*/

    if(k <= pivot_left.size()) {
        return quickSelection(pivot_left, k);
    } else if(k >= pivot_left.size() + 1 && k <= pivot_left.size() + same) {
        return pivot; 
    } else {
        return quickSelection(pivot_right, k - pivot_left.size() - same);
    }
}

std::vector<std::vector<double> > generateAffinityMatrix(std::vector<Point> *points) {
    std::vector<std::vector<double> > result(points->size());
    std::vector<std::vector<double> > distance(points->size());
    std::vector<double> deltas(points->size());

    for(int p1 = 0; p1 < points->size(); p1++) {
        std::vector<double> partialDistance(points->size());

        for(int p2 = 0; p2 < points->size(); p2++) {
            partialDistance[p2] = getDistance((*points)[p1], (*points)[p2]);
        }

        distance[p1] = partialDistance;
        deltas[p1] = quickSelection(partialDistance, SIGMA_DIMENSION);
    }

    for(int p1 = 0; p1 < points->size(); p1++) {
        std::vector<double> partialResult(points->size());

        for(int p2 = 0; p2 < points->size(); p2++) {
            double affinity = exp(-distance[p1][p2] / (2 * deltas[p1] * deltas[p2]));
            partialResult[p2] = affinity;
        }

        result[p1] = partialResult;
    }

    return result;
}

std::vector<std::vector<double> > generateAffinityMatrix_parallel(std::vector<Point> *points) {
    std::vector<std::vector<double> > result(points->size());
    std::vector<std::vector<double> > distance(points->size());
    std::vector<double> deltas(points->size());

    #pragma omp parallel num_threads(THREADS)
    {
        #pragma omp for schedule(static, points->size() / THREADS)
        for(int p1 = 0; p1 < points->size(); p1++) {
            std::vector<double> partialDistance(points->size());

            for(int p2 = 0; p2 < points->size(); p2++) {
                partialDistance[p2] = getDistance((*points)[p1], (*points)[p2]);
            }

            distance[p1] = partialDistance;
            deltas[p1] = quickSelection(partialDistance, SIGMA_DIMENSION);
        }

        #pragma omp for schedule(static, points->size() / THREADS)
        for(int p1 = 0; p1 < points->size(); p1++) {
            std::vector<double> partialResult(points->size());

            for(int p2 = 0; p2 < points->size(); p2++) {
                double affinity = exp(-distance[p1][p2] / (2 * deltas[p1] * deltas[p2]));
                partialResult[p2] = affinity;
            }

            result[p1] = partialResult;
        }
    }

    return result;
}

int main() {
    double arr[] = {5.974688, 9.825365, 10.688012, 7.398466, 9.35243, 5.811162, 9.177718, 3.696228, 11.196732, 0.000000};
    int n = sizeof(arr) / sizeof(arr[0]);
    std::vector<double> vec(arr, arr + n);

    //printf("%3lf\n", quickSelection(vec, 5));

    /*for(int i = 0; i < n; i++) {
        printf("%3lf\n", quickSelection(vec, i + 1));
    }*/

    printf("Data size: %d", TEST_DATA_COUNT);

    DS_timer timer(2);
	timer.setTimerName(0, (char*)"Serial");
	timer.setTimerName(1, (char*)"Parallel");

    std::vector<Point> points;
    for(int i = 0; i < TEST_DATA_COUNT / 3; i++) {
        Point p = Point(GenDouble, GenDouble);
        points.push_back(p);
    }

    for(int i = 0; i < TEST_DATA_COUNT / 3; i++) {
        Point p = Point(GenDouble + 20, GenDouble + 10);
        points.push_back(p);
    }

    for(int i = 0; i < TEST_DATA_COUNT / 3; i++) {
        Point p = Point(GenDouble - 20 , GenDouble - 20);
        points.push_back(p);
    }

    timer.onTimer(0);
    std::vector<std::vector<double> > result = generateAffinityMatrix(&points);
    timer.offTimer(0);

    timer.onTimer(1);
    std::vector<std::vector<double> > result2 = generateAffinityMatrix_parallel(&points);
    timer.offTimer(1);

    timer.printTimer();
}

