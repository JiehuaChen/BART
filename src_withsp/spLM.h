#ifndef HEADER_SPLM_H
#define HEADER_SPLM_H

void spLM(double* Y_r, int n_r, double* coordsD_r, 
        double* sigmaSqIG_r, double* tauSqIG_r, double* nuUnif_r, 
        double* phiUnif_r, double* phiStarting_r, double* sigmaSqStarting_o, 
        double tauSqStarting_r, double* nuStarting_r, double* phiTuning_r, 
        double* sigmaSqTuning_r, double* tauSqTuning_r, double* nuTuning_r, 
        bool nugget_r, std::string covModel_r, bool amcmc_r, 
        int nBatch_r, int batchLength_r, double acceptRate_r, 
		int nParams, int sigmaSqIndx, int tauSqIndx, int phiIndx, int nuIndx,
        int verbose, int nReport_r, double *spdraw, double *params);

#endif
