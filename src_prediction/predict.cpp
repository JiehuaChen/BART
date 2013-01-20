#include <iostream>
#include <stdio.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <fstream>

#include "Node.h"
#include "Trees.h"
#include "Utils.h"


using namespace std;
int NumX;

int main(int argc, const char *argv[])
{
	NumX = 20;
	double * p_OneX = NULL;
	double OneY = 0.0;
	string testdatafile("../MCMCresults/xdat.txt");
	string treefile;
	std::vector< std::vector<double>* > * p_vv_testdata;

	treefile = "../MCMCresults/MCMC1100.txt";
	Trees *p_tree = new Trees;
	p_tree->ParseTreesFile(treefile);
	p_tree->BuildTrees();
//	p_tree->PrintTrees();

	string predictedfileName = "predictedY.txt";
	FILE*  predictedfile= fopen(predictedfileName.c_str(),"w+t"); // append mode

	p_vv_testdata = GetTestData(testdatafile);

	for(unsigned int i=1;i<p_vv_testdata->size();i++)
	{
//		printf("%u:%u\n",p_vv_testdata->size(),p_vv_testdata->at(i)->size());
		p_OneX = new double[p_vv_testdata->at(i)->size()];
		for(unsigned int j=1; j<p_vv_testdata->at(i)->size();j++)
		{
			p_OneX[j] = p_vv_testdata->at(i)->at(j);
		}
		OneY = p_tree->GetPredictedResult(p_vv_testdata->at(i)->size(), p_OneX);
		fprintf(predictedfile,"%0.17f\n",OneY);
		printf("PredictedY: %f\n",OneY);
		delete p_OneX;
	}
	
	fclose(predictedfile);
	delete p_tree;
	return 1;
}
