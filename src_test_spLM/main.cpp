#include <iostream>
#include <stdio.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <time.h>
#include <math.h>
#include <glob.h>
#include <string>
#include <numeric>
#include <algorithm>
#include "testSPLM.h"

using namespace std;

int main(int argc, char **argv)
{	
	testSPLM* tester = new testSPLM();
	tester->start("spLM_input.txt");
	return 0;
}
