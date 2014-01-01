#ifndef GUARD_UTILS_H
#define GUARD_UTILS_H

#include <vector>
#include <string>
#include <cstring>

using namespace std;
#define DIMENSION 100
class testSPLM {
		public:
	    testSPLM();
		~testSPLM();
		int start(string inputfilename);
		int GetTestData(string infile);
		int ParseOneInput(string &oneline);
		int GetRangeData(string infile);

		double Y_r[DIMENSION];
		double X_r[DIMENSION];
		int p_r;
		int n_r;
		double *coordsD_r; 
		std::string betaPrior_r;
		double sigmaSqIG_r[2];
		double tauSqIG_r[2];
		double nuUnif_r[2]; 
		double phiUnif_r[2];
		double phiStarting_r;
		double sigmaSqStarting_r; 
		double tauSqStarting_r;
		double nuStarting_r;
		double phiTuning_r; 
		double sigmaSqTuning_r; 
		double tauSqTuning_r; 
		double nuTuning_r;
		bool nugget_r; 
		std::string covModel_r;
		bool amcmc_r; 
		int nBatch_r; 
		int batchLength_r;
		double acceptRate_r; 
		int verbose_r;
		int nReport_r;
		double **betaNorm_r;
		double **result;

};
#endif
