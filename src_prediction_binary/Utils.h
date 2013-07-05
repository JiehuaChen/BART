#ifndef GUARD_UTILS_H
#define GUARD_UTILS_H

#include <vector>
#include <string>
#include <cstring>

using namespace std;

std::vector< std::vector<double>* > * GetTestData(string infile);
std::vector<double>* ParseOneInput(string &oneline);
std::vector<double>* GetRangeData(string infile);

#endif
