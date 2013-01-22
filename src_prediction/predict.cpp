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
#include "Node.h"
#include "Trees.h"
#include "Utils.h"

using namespace std;
int NumX;

void
usage (void) {
	printf ("\nNAME:\n");
	printf ("\tpredict - calculates BART predictions, given explanatory variables and trees\n");
	printf ("\nSYNOPSIS:\n");
	printf ("\t./predict -x [explanatory_data_file] -f [forests_glob_pattern] -r [range_file] -o [predictions_output_file]\n");
	printf ("\nOPTIONS:\n");
	printf ("\t-x [filename]\n");
	printf ("\t\tSpecify input data file of explanatory variables.\n");
	printf ("\t-f [glob pattern]\n");
	printf ("\t\tSpecify input glob pattern matching all files of forests. Surround glob pattern with quotation marks.\n");
	printf ("\t-r [filename]\n");
	printf ("\t\tSpecify input range file of range bounds.\n");
	printf ("\t-o [filename]\n");
	printf ("\t\tSpecify output data file of predictions.\n");
	printf ("\t-s [separator string]\n");
	printf ("\t\tSpecify seperator string to use in output data file. Defaults to ',' for CSV format.\n");
	printf ("\t\tNote: For tab-separated values, type -s \"	\", where (Control+V followed by TAB) is supplied in between the quotation marks.\n");
	printf ("\t-v\n");
	printf ("\t\tPrints verbose output.\n");
	printf ("\t-h\n");
	printf ("\t\tDisplays this help message.\n");
	printf ("\nEXAMPLE:\n");
	printf ("\t./predict -x ../MCMCresults/xdat.txt -f \"../MCMCresults/MCMC*.txt\" -s \" \" -r ../MCMCresults/rgy.txt -o predictedY.txt\n");
	printf ("\nCONTACT:\n");
	printf ("\tJiehua Chen <jc3288@columbia.edu>, 2013-01-20\n");
}

bool file_read_check(const char *filename) {
	std::ifstream ifile(filename);
	return ifile;
}
bool file_write_check(const char *filename) {
	std::ofstream ofile(filename);
	return ofile;
}
inline double range_correction(double y, double lower, double upper) {
	return ((upper-lower)*(y + .5) + lower);
}
inline std::vector<std::string> glob(const std::string& pat){
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}

int main(int argc, char **argv)
{	
	// timing variables
	double t1_read = 0.0, t2_read = 0.0;
	double t1_prediction = 0.0, t2_prediction = 0.0;	
	
	/************************************************************
    COMMAND-LINE ARGUMENTS
	************************************************************/
	int g;
	int xflag = 0, fflag = 0, rflag = 0, oflag = 0, sflag = 0, vflag = 0;
	// defaults
	const char *xpath = "../MCMCresults/xdat.txt";
	const char *fpath = "../MCMCresults/MCMC*.txt";
	const char *rpath = "../MCMCresults/rgy.txt";
	const char *opath = "predictedY.txt";
	const char *sep = ","; // csv format
	// options
	while ((g = getopt(argc,argv,"x:f:r:o:s:vh")) != -1) {
		switch (g) {
			case 'x': {
				xflag=1;
				xpath = optarg;
				break;
			}
			case 'f': {
				fflag=1;
				fpath = optarg;
				break;
			}
			case 'r': {
				rflag=1;
				rpath = optarg;
				break;
			}
			case 'o': {
				oflag=1;
				opath = optarg;
				break;
			}
			case 's': {
				sflag=1;
				sep = optarg;
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
		if (vflag) printf ("Explanatory variables file (input): %s\n",xpath);
	}
	if (0==fflag) {
		printf ("No forests glob specified as input -- defaulting to %s\n",fpath);	
	} else {
		if (vflag) printf ("Forests glob (input): %s\n",fpath);	
	}
	if (0==rflag) {
		printf ("No range file (input) specified -- defaulting to %s.\n",rpath);
	} else {
		if (vflag) printf ("Range file (input): %s\n",rpath);	
	}
	if (0==oflag) {
		printf ("No predictions file (output) specified -- defaulting to %s.\n",opath);
	} else {
		if (vflag) printf ("Predictions file (output): %s\n",opath);	
	}
	// argument checking
	if (!(file_read_check(xpath))) { 
		printf ("Error: Cannot open %s for reading.", xpath);   
		return -1;
	}
	if (!(file_read_check(rpath))) { 
		printf ("Error: Cannot open %s for reading.", rpath);   
		return -2;
	}
	if (!(file_write_check(opath))) { 
		printf ("Error: Cannot open %s for writing.", opath);   
		return -3;
	}
	
	/************************************************************
    READ EXPLANATORY VARIABLES, TREES, AND RANGE BOUNDS
	************************************************************/
	t1_read = clock();

	NumX = 20;
	double * p_OneX = NULL;
	double OneY = 0.0;
	double OneY_corrected = 0.0;

	// read explanatory variables
	string testdatafile(xpath);
	std::vector< std::vector<double>* > * p_vv_testdata;
	p_vv_testdata = GetTestData(testdatafile);

	// read range 
	string rangefile(rpath);
	std::vector<double>* p_range;
	p_range = GetRangeData(rangefile);
	double range_low = p_range->at(0);
	double range_high = p_range->at(1);

	t2_read = clock();
	
	/************************************************************
    COMPUTE PREDICTIONS
	************************************************************/
	std::vector<std::string> forestsfilenames = glob(fpath); // fetch filenames of forests
	std::sort(forestsfilenames.begin(), forestsfilenames.end());

	std::vector<std::vector<double> > prediction_matrix;
	std::vector<double> prediction_times;		
	std::vector<int> forest_sizes;
	
	for (unsigned int k=0; k<forestsfilenames.size(); k++)  // iterate over forests
	{
		if (vflag) {
			printf("=============================================\n");
			printf("FOREST: %s\n", forestsfilenames[k].c_str());
			printf("=============================================\n");
			printf("Number of predicted values: %d\n", (int) (p_vv_testdata->size() - 1));
		}

		t1_prediction = clock();
		prediction_matrix.push_back(std::vector<double>());
		
		Trees *p_tree = new Trees;
		p_tree->ParseTreesFile(forestsfilenames[k]);
		p_tree->BuildTrees();		
		forest_sizes.push_back(p_tree->theTrees.size() - 1);
		
		for(unsigned int i=1;i<p_vv_testdata->size();i++) // iterate over variables
		{
	//		printf("%u:%u\n",p_vv_testdata->size(),p_vv_testdata->at(i)->size());
			p_OneX = new double[p_vv_testdata->at(i)->size()];
			for(unsigned int j=1; j<p_vv_testdata->at(i)->size();j++)
			{
				p_OneX[j] = p_vv_testdata->at(i)->at(j);
			}
			OneY = p_tree->GetPredictedResult(p_vv_testdata->at(i)->size(), p_OneX); // prediction
			OneY_corrected = range_correction(OneY, range_low, range_high); // range correction
			prediction_matrix[k].push_back(OneY_corrected);
			
			if (vflag) {
				printf("Predicted and range-corrected Y[%d]: %f --> %f\n",i,OneY,OneY_corrected); // print to screen
			}
			delete p_OneX; // free memory
		}
		delete p_tree;

		t2_prediction = clock();
		prediction_times.push_back((t2_prediction - t1_prediction) / CLOCKS_PER_SEC);
	}

	/************************************************************
    WRITE PREDICTIONS TO DISK
	************************************************************/

	// transpose prediction_matrix
	int F = (int)prediction_matrix.size(); // number of variables
	int V = (int)prediction_matrix[0].size(); // number of forests
	std::vector<std::vector<double> > A;
	A.resize(V); // V x F matrix
	for (int i = 0; i < V; i++) A[i].resize(F);
	for (int i = 0; i < F; i++) {
		for (int j = 0; j < V; j++) {
			A[j][i] = prediction_matrix[i][j];
		}
	}

	// initializations
	cout << "\nPredictions complete." << endl;
	string predictedfileName(opath);
	FILE* predictedfile = fopen(predictedfileName.c_str(),"w+t"); // overwrite mode
	double mean, deviation;	
	
	// print headers
	for (unsigned int j=0; j<A[0].size(); j++) {
		fprintf(predictedfile,"draw %d",j+1);
		fprintf(predictedfile,"%s",sep);
		if (j==A[0].size()-1) { // compute mean and standard deviation
			fprintf(predictedfile,"mean");
			fprintf(predictedfile,"%s",sep);
			fprintf(predictedfile,"std");
		}
	}
	fprintf(predictedfile,"\n");
	
	// print data
	for (unsigned int i=0; i<A.size(); i++) {
		for (unsigned int j=0; j<A[0].size(); j++) {
			fprintf(predictedfile,"%0.6f",A[i][j]);
			fprintf(predictedfile,"%s",sep);

			if (j==A[0].size()-1) { // compute mean and standard deviation
				mean = accumulate( A[i].begin(), A[i].end(), 0.0f )/ A[i].size();
				vector<double> zero_mean( A[i] );
				transform( zero_mean.begin(), zero_mean.end(), zero_mean.begin(),bind2nd( minus<float>(), mean ) );
			   	deviation = inner_product( zero_mean.begin(),zero_mean.end(), zero_mean.begin(), 0.0f );
				deviation = sqrt( deviation / ( A[i].size() - 1 ) );

				fprintf(predictedfile,"%0.6f",mean);
				fprintf(predictedfile,"%s",sep);
				fprintf(predictedfile,"%0.6f",deviation);
			}
		}
		fprintf(predictedfile,"\n");
	}
	if (vflag) {
		cout << "Range-corrected predictions using the forests in: " << endl;
		for (unsigned int k=0; k<forestsfilenames.size(); k++) {
			cout << "\t" << forestsfilenames[k] << endl;
		}
		cout << "have been written to disk as " << predictedfileName << endl;
	}

	fclose(predictedfile);

	/************************************************************
    REPORT TIMES
	************************************************************/
	printf("%35s: %15.6lf seconds\n", "Reading time", (t2_read - t1_read) / CLOCKS_PER_SEC);
	char str[200];
	for (unsigned int k=0; k<prediction_times.size(); k++) {
		sprintf(str,"Prediction time for Forest #%d",k);
		printf("%35s: %15.6lf seconds\n", str, prediction_times[k]);		
	}

	printf("Number of trees in the BART model: %d.\n", forest_sizes[0]); // assumes that all forests have same # of trees
	printf("Output is a %d-by-%d matrix, where each posterior prediction is written column by column.\n\n", (int) A.size(), (int) A[0].size());
	
	return 0;
}
