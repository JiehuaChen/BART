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

void
usage (void) {
	printf ("NAME:\n");
	printf ("\tpredict - calculates BART predictions, given explanatory variables and trees\n");
	printf ("SYNOPSIS:\n");
	printf ("\t./predict -x [explanatory_data_file] -t [trees_file] -o [predictions_output_file]\n");
	printf ("OPTIONS:\n");
	printf ("\t-x [filename]\n");
	printf ("\t\tSpecify input data file of explanatory variables.\n");
	printf ("\t-t [filename]\n");
	printf ("\t\tSpecify input data file of trees.\n");
	printf ("\t-o [filename]\n");
	printf ("\t\tSpecify output data file of predictions.\n");
	printf ("\t-v\n");
	printf ("\t\tPrints verbose output.\n");
	printf ("\t-h\n");
	printf ("\t\tDisplays this help message.\n");
	printf ("Primary Contact:\n");
	printf ("\tJiehua Chen <Jiehua.Chen@themathpath.com>, 2012\n");
}

bool file_read_check(const char *filename) {
	std::ifstream ifile(filename);
	return ifile;
}
bool file_write_check(const char *filename) {
	std::ofstream ofile(filename);
	return ofile;
}

int main(int argc, char **argv)
{
	/**********************************************************************
    COMMAND-LINE ARGUMENTS
	**********************************************************************/
	int g;
	int xflag = 0, tflag = 0, oflag = 0, vflag = 0;
    // defaults
    const char *xpath = "../MCMCresults/xdat.txt";
    const char *tpath = "../MCMCresults/MCMC1100.txt";
	const char *opath = "predictedY.txt";
	// options
	while ((g = getopt(argc,argv,"x:t:vh")) != -1) {
		switch (g) {
			case 'x': {
				xflag=1;
				xpath = optarg;
				break;
			}
			case 't': {
				tflag=1;
				tpath = optarg;
				break;
			}
			case 'o': {
				oflag=1;
				opath = optarg;
				break;
			}
			case 'v': {
				vflag=1;
				break;
			}
			case 'h': {
				usage ();
				return 1;
				break;
			}
			case '?': {
				usage ();
				return 1;
				break;
			}
			default: {
				usage ();
				return 1;
			}			
		}
	}
	if (0==xflag) {
		printf ("No explanatory variables file specified as input -- defaulting to %s\n",xpath);
	} else {
		printf ("Explanatory variables file (input): %s\n",xpath);
	}
	if (0==tflag) {
		printf ("No trees file specified as input -- defaulting to %s\n",tpath);	
	} else {
		printf ("Trees file (input): %s\n",tpath);	
	}
	if (0==oflag) {
		printf ("No predictions file (output) specified -- defaulting to %s.\n",opath);
	} else {
		printf ("Predictions file (output): %s\n",opath);	
	}
	// argument checking
	if (!(file_read_check(xpath))) { 
		printf ("Error: Cannot open %s for reading.", xpath);   
		return -1;
	}
	if (!(file_read_check(tpath))) { 
		printf ("Error: Cannot open %s for reading.", tpath);   
		return -2;
	}
	if (!(file_write_check(opath))) { 
		printf ("Error: Cannot open %s for writing.", opath);   
		return -3;
	}

	/**********************************************************************
    READ EXPLANATORY VARIABLES AND TREES
	**********************************************************************/
	NumX = 20;
	double * p_OneX = NULL;
	double OneY = 0.0;

	string testdatafile(xpath);
	string treefile(tpath);
	std::vector< std::vector<double>* > * p_vv_testdata;

	Trees *p_tree = new Trees;
	p_tree->ParseTreesFile(treefile);
	p_tree->BuildTrees();
//	p_tree->PrintTrees();

	string predictedfileName(opath);
	FILE* predictedfile= fopen(predictedfileName.c_str(),"w+t"); // overwrite mode

	p_vv_testdata = GetTestData(testdatafile);

	/**********************************************************************
    COMPUTE PREDICTIONS
	**********************************************************************/
	printf("Number of predicted values: %d\n", (int) (p_vv_testdata->size() - 1));
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
		printf("Predicted Y[%d]: %f\n",i,OneY);
		delete p_OneX;
	}
	
	fclose(predictedfile);
	delete p_tree;
	return 1;
}