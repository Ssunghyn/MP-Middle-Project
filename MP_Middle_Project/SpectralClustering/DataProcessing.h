#include <stdio.h>
#include "SpectralClustring.h"


bool getData(FILE* fp, Point* points);
void saveData(const char* fileName, double* points, int n);