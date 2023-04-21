#include "DataProcessing.h"
#include <stdio.h>

using namespace std;


bool getData(FILE* fp, Point* points, int n) {
	if (fp == NULL) {
		printf("FILE NOT OPEN!\n");
		return;
	}

    int i = 0;
	while (!feof(fp)) {
		double x = 0;
		double y = 0;
		
		fscanf_s(fp, "%lf %lf", &x, &y);
		points[i] = Point(x, y);
        i++;
	}

	fclose(fp);

	if (i == n)
		return true;
	
	return false;
}

void saveData(const char* fileName, double* points, int n) {
	FILE* fp;
	fopen_s(&fp,fileName, "w");

	if (fp == NULL) {
		printf("FILE NOT OPEN!\n");
		return;
	}

	for(int i = 0; i < n; i++)
		fprintf(fp, "%lf\n", points[i]);

	fclose(fp);
}


void data_main()
//int main() 
{
	const char* input_path = "circle_data.txt";
    int n = 10000;
	FILE* fp;
	fopen_s(&fp, input_path, "r");
	Point* points = new Point[n];
	getData(fp, points);
	for (int i = 0; i < n; i++)
		printf("%lf %lf\n", points[i].x, points[i].y);
}