#include <stdio.h>
#include "SpectralClustring.h"


bool getData(FILE* fp, Point* points, int n);
void saveData(const char* fileName, Point* points, double* results, int n);