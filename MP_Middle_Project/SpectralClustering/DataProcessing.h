// #pragma once
#include <stdio.h>
#include "SpectralClustring.h"
#include <vector>
#include <string>

FILE* fp;
std::string input_path;  // 불러올 파일 경로

void getData(FILE* fp, std::vector<Point>& points);