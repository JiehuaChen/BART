#ifndef HEADER_SPLM_H
#define HEADER_SPLM_H

int spLM(double* Y_r, double* X_r, int p_r, int n_r, double* coordsD_r, 
	std::string betaPrior_r, double** betaNorm_r, double* sigmaSqIG_r, double* tauSqIG_r, double* nuUnif_r, double* phiUnif_r,
	double phiStarting_r, double sigmaSqStarting_r, double tauSqStarting_r, double nuStarting_r,
	double phiTuning_r, double sigmaSqTuning_r, double tauSqTuning_r, double nuTuning_r, 
	bool nugget_r, std::string covModel_r, bool amcmc_r, int nBatch_r, int batchLength_r, double acceptRate_r, 
	int verbose_r, int nReport_r, double **result);

#endif
