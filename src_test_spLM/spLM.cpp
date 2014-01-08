#include <algorithm>
#include <string>
extern "C" {
#include <R_ext/Print.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R.h>
#include <Rmath.h>
#include <S.h>
};
#include "util.h"
#include "spLM.h"


#define Rprintf printf
#define R_alloc calloc

void spLM(double* Y_r, int n_r, double* coordsD_r, 
        double* sigmaSqIG_r, double* tauSqIG_r, double* nuUnif_r, 
        double* phiUnif_r, double phiStarting_r, double sigmaSqStarting_r, 
        double tauSqStarting_r, double nuStarting_r, double* phiTuning_r, 
        double* sigmaSqTuning_r, double* tauSqTuning_r, double* nuTuning_r, 
        bool nugget_r, std::string covModel_r, bool amcmc_r, 
        int nBatch_r, int batchLength_r, double acceptRate_r, 
		int nParams, int sigmaSqIndx, int tauSqIndx, int phiIndx, int nuIndx,
        int verbose, int nReport_r, double* spdraw, double* params) {

		/*****************************************
		  Common variables
		 *****************************************/
		int i, j, k, b, s, info, n_accept_r = 0;
		char const *lower = "L";
		char const *upper = "U";
		char const *nUnit = "N";
		char const *yUnit = "U";
		char const *ntran = "N";
		char const *ytran = "T";
		char const *rside = "R";
		char const *lside = "L";
		const double one = 1.0;
		const double negOne = -1.0;
		const double zero = 0.0;
		const int incOne = 1;

		/*****************************************
		  Set-up
		 *****************************************/
		double *Y = Y_r;
		int n = n_r;
		int nn = n*n;

		//priors
		double sigmaSqIGa = sigmaSqIG_r[0]; double sigmaSqIGb = sigmaSqIG_r[1];
		double phiUnifa = phiUnif_r[0]; double phiUnifb = phiUnif_r[1];

		bool nugget = nugget_r;
		double tauSqIGa = 0, tauSqIGb = 0;
		if(nugget){
			tauSqIGa = tauSqIG_r[0]; tauSqIGb = tauSqIG_r[1]; 
		}

		//matern
		double nuUnifa = 0, nuUnifb = 0;
		if(covModel == "matern"){
			nuUnifa = nuUnif_r[0]; nuUnifb = nuUnif_r[1]; 
		}

        // adaptive mcmc
		bool amcmc = amcmc_r;
		int nBatch = nBatch_r;
		// number of MCMC samplings
		int batchLength = batchLength_r;
		double acceptRate = acceptRate_r;
		int nSamples = nBatch*batchLength;
		int nReport = nReport_r;

		// if(verbose){
		// 	Rprintf("----------------------------------------\n");
		// 	Rprintf("\tGeneral model description\n");
		// 	Rprintf("----------------------------------------\n");
		// 	Rprintf("Model fit with %i observations.\n\n", n);
		// 	Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());

		// 	if(amcmc){
		// 		Rprintf("Using adaptive MCMC.\n\n");
		// 		Rprintf("\tNumber of batches %i.\n", nBatch);
		// 		Rprintf("\tBatch length %i.\n", batchLength);
		// 		Rprintf("\tTarget acceptance rate %.5f.\n", acceptRate);
		// 		Rprintf("\n");
		// 	}else{
		// 		Rprintf("Number of MCMC samples %i.\n\n", nSamples);
		// 	}

		// 	if(!nugget){
		// 		Rprintf("tau.sq not included in the model (i.e., no nugget model).\n\n");
		// 	}

		// 	Rprintf("Priors and hyperpriors:\n");
		// 	Rprintf("\tsigma.sq IG hyperpriors shape=%.5f and scale=%.5f\n", sigmaSqIGa, sigmaSqIGb);
		// 	if(nugget){
		// 		Rprintf("\ttau.sq IG hyperpriors shape=%.5f and scale=%.5f\n", tauSqIGa, tauSqIGb); 
		// 	}
		// 	Rprintf("\tphi Unif hyperpriors a=%.5f and b=%.5f\n", phiUnifa, phiUnifb);
		// 	if(covModel == "matern"){
		// 		Rprintf("\tnu Unif hyperpriors a=%.5f and b=%.5f\n", nuUnifa, nuUnifb);	  
		// 	}
		// } 

		/*****************************************
		  Set-up MCMC sample matrices etc.
		 *****************************************/ 
		//parameters: nParams number of parameters


        // the vector of parameters of the spatial model
        double  *params = new double[nParam];

		//starting
		params[sigmaSqIndx] = log(sigmaSqStarting_r);

		if(nugget){
			params[tauSqIndx] = log(tauSqStarting_r);
		}

		params[phiIndx] = logit(phiStarting_r, phiUnifa, phiUnifb);

		if(covModel == "matern"){
			params[nuIndx] = logit(nuStarting_r, nuUnifa, nuUnifb);
		}

		//tuning and fixed
		double *tuning = new double[nParams];
		int *fixed = new int[nParams]; zeros(fixed, nParams);

		tuning[sigmaSqIndx] = *sigmaSqTuning_r;
		if(tuning[sigmaSqIndx] == 0){
			fixed[sigmaSqIndx] = 1;
		}

		if(nugget){
			tuning[tauSqIndx] = *tauSqTuning_r;
			if(tuning[tauSqIndx] == 0){
				fixed[tauSqIndx] = 1;
			}
		}

		tuning[phiIndx] = *phiTuning_r;
		if(tuning[phiIndx] == 0){
			fixed[phiIndx] = 1;
		}

		if(covModel == "matern"){
			tuning[nuIndx] = *nuTuning_r;
			if(tuning[nuIndx] == 0){
				fixed[nuIndx] = 1;
			}
		}

		for(i = 0; i < nParams; i++){
			tuning[i] = log(sqrt(tuning[i]));
		}

		//return values: 
		//samples_r: estimated parameters;
		//accept_r: acceptance rate;
		//tuning_r: jumping distance for Methoplis steps
		double *samples_r = new double[nParams * nSamples];
		
		if(amcmc){
			double *accept_r = new double[nParams * nBatch];
			double *tuning_r = new double[nParams * nBatch];
		}else{
			n_accept_r = floor(static_cast<double>(nSamples/nReport));
			double *accept_r = new double[floor(static_cast<double>(nSamples/nReport))];
		}

		/*****************************************
		  Set-up MCMC alg. vars. matrices etc.
		 *****************************************/

		int status=1, batchAccept=0, reportCnt=0;
		double logMHRatio =0, logPostCurrent = R_NegInf, logPostCand = 0, det = 0, paramsjCurrent = 0;

		double *paramsCurrent = new double[nParams];
		double *accept = new double[nParams]; zeros(accept, nParams);

        // C: covariance matrix
		double *C = new double[nn]; zeros(C, nn);
		double sigmaSq, phi, tauSq, nu, Q;
		double *theta = new double[3]; //phi, nu, and perhaps more in the future

		int p1 = 1;
		double *vU = new double[n*p1];


		if(verbose){
			Rprintf("-------------------------------------------------\n");
			Rprintf("\t\tSampling\n");
			Rprintf("-------------------------------------------------\n");
#ifdef Win32
			R_FlushConsole();
#endif
		}

		// GetRNGstate();

		for(b = 0, s = 0; b < nBatch; b++){    
			for(i = 0; i < batchLength; i++, s++){
				for(j = 0; j < nParams; j++){
					//propose
					if(amcmc){
						if(fixed[j] == 1){
							paramsjCurrent = params[j];
						}else{
							paramsjCurrent = params[j];
							params[j] = rnorm(paramsjCurrent, exp(tuning[j]));
						}
					}else{
						F77_NAME(dcopy)(&nParams, params, &incOne, paramsCurrent, &incOne);
						for(j = 0; j < nParams; j++){
							if(fixed[j] == 1){
								params[j] = params[j];
							}else{
								params[j] = rnorm(params[j], exp(tuning[j]));
							}
						}
					}
					//extract and transform
					sigmaSq = theta[0] = exp(params[sigmaSqIndx]);
					phi = theta[1] = logitInv(params[phiIndx], phiUnifa, phiUnifb);

					if(nugget){
						tauSq = exp(params[tauSqIndx]);
					}

					if(covModel == "matern"){
						nu = theta[2] = logitInv(params[nuIndx], nuUnifa, nuUnifb);
					}

					//construct covariance matrix: create the covariance matrix C
					// C is a n(n-1)/2 by 1 array
					spCovLT(coordsD, n, theta, covModel, C);

					if(nugget){
						for(k = 0; k < n; k++){
							C[k*n+k] += tauSq;
						}
					}

					det = 0;
					// dpotrf - computes the Cholesky factorization of a real symmetric positive definite matrix C to be stored in C
					F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
					for(k = 0; k < n; k++) det += 2*log(C[k*n+k]);

					//Q = Y^T C^{-1} Y
					F77_NAME(dcopy)(&n, Y, &incOne, vU, &incOne);
					F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &p1, &one, C, &n, vU, &n);//Lv = y
					Q = pow(F77_NAME(dnrm2)(&n, vU, &incOne),2) ;

					//
					//priors, jacobian adjustments, and likelihood
					//
					logPostCand = -1.0*(1.0+sigmaSqIGa)*log(sigmaSq)-sigmaSqIGb/sigmaSq+log(sigmaSq);

					if(nugget){
						logPostCand += -1.0*(1.0+tauSqIGa)*log(tauSq)-tauSqIGb/tauSq+log(tauSq);
					}

					logPostCand += log(phi - phiUnifa) + log(phiUnifb - phi); 

					if(covModel == "matern"){
						logPostCand += log(nu - nuUnifa) + log(nuUnifb - nu);   
					}

					logPostCand += -0.5*det-0.5*Q;

					//
					//MH accept/reject	
					//      
					logMHRatio = logPostCand - logPostCurrent;

					if(runif(0.0,1.0) <= exp(logMHRatio)){
						logPostCurrent = logPostCand;
						if(amcmc){
							accept[j]++;
						}else{
							accept[0]++;
							batchAccept++;
						}
					}else{
						if(amcmc){
							params[j] = paramsjCurrent;
						}else{
							F77_NAME(dcopy)(&nParams, paramsCurrent, &incOne, params, &incOne);
						}
					}

					if(!amcmc){
						break;
					}
				}//end params

				/******************************
				  Save samples
				 *******************************/
				F77_NAME(dcopy)(&nParams, params, &incOne, &samples_r[s*nParams], &incOne);

				R_CheckUserInterrupt();
			}//end batch

			//adjust tuning
			if(amcmc){
				for(j = 0; j < nParams; j++){
					accept_r[b*nParams+j] = accept[j]/batchLength;
					tuning_r[b*nParams+j] = exp(tuning[j]);

					if(accept[j]/batchLength > acceptRate){
						tuning[j] += std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
					}else{
						tuning[j] -= std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
					}
					accept[j] = 0.0;
				}
			}
			//report
			if(status == nReport){
				if(verbose){
					if(amcmc){
						Rprintf("Batch: %i of %i, %3.2f%%\n", b+1, nBatch, 100.0*(b+1)/nBatch);
						Rprintf("\tparameter\tacceptance\ttuning\n");	  
						Rprintf("\tsigma.sq\t%3.1f\t\t%1.5f\n", 100.0*accept_r[b*nParams+sigmaSqIndx], exp(tuning[sigmaSqIndx]));
						if(nugget){
							Rprintf("\ttau.sq\t\t%3.1f\t\t%1.5f\n", 100.0*accept_r[b*nParams+tauSqIndx], exp(tuning[tauSqIndx]));
						}
						Rprintf("\tphi\t\t%3.1f\t\t%1.5f\n", 100.0*accept_r[b*nParams+phiIndx], exp(tuning[phiIndx]));
						if(covModel == "matern"){
							Rprintf("\tnu\t\t%3.1f\t\t%1.5f\n", 100.0*accept_r[b*nParams+nuIndx], exp(tuning[nuIndx]));
						}
					}else{
						Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
						Rprintf("Report interval Metrop. Acceptance rate: %3.2f%%\n", 100.0*batchAccept/nReport);
						Rprintf("Overall Metrop. Acceptance rate: %3.2f%%\n", 100.0*accept[0]/s);
					}
					Rprintf("-------------------------------------------------\n");
#ifdef Win32
					R_FlushConsole();
#endif
				}
				if(!amcmc){
					accept_r[reportCnt] = 100.0*batchAccept/nReport;
					reportCnt++;
				}
				batchAccept = 0;
				status = 0;
			}
			status++;
		}//end sample loop

		// PutRNGstate();

		//untransform variance variables
		for(s = 0; s < nSamples; s++){
			samples_r[s*nParams+sigmaSqIndx] = exp(samples_r[s*nParams+sigmaSqIndx]);
			if(nugget){
				samples_r[s*nParams+tauSqIndx] = exp(samples_r[s*nParams+tauSqIndx]);
			}
			samples_r[s*nParams+phiIndx] = logitInv(samples_r[s*nParams+phiIndx], phiUnifa, phiUnifb);
			if(covModel == "matern")
				samples_r[s*nParams+nuIndx] = logitInv(samples_r[s*nParams+nuIndx], nuUnifa, nuUnifb);
		}

        //make return object

		//samples
        // params

		//simulated spatial process: last draw of spatial parameters
        sigmaSq = theta[0] = exp(params[sigmaSqIndx]);
        phi = theta[1] = logitInv(params[phiIndx], phiUnifa, phiUnifb);

        if(nugget){
            tauSq = exp(params[tauSqIndx]);
        }

        if(covModel == "matern"){
            nu = theta[2] = logitInv(params[nuIndx], nuUnifa, nuUnifb);
        }

        //construct covariance matrix: create the covariance matrix C
        // C is a n(n-1)/2 by 1 array

        spCovLT(coordsD, n, theta, covModel, C);
        if(nugget){
      	    for(k = 0; k < n; k++){
                C[k*n+k] += tauSq[s];
            }
        } 
        // cholesky decomposition of C
        F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}//L_1

		double *spdraw = new double[n];
		int *betamu = new int[n]; zeros(betamu, n);
        mvrnorm(spdraw, betamu, C, n);
        
		if(amcmc){
            for(j = 0;j < nParams;j++){
                tuning[j] = exp(tuning[j]);
            }
			
            sigmaSqTuning_r = &tuning[sigmaSqIndx];

            if(nugget){
                tauSqTuning_r = &tuning[tauSqIndx];
            }

            tuning[phiIndx] = *phiTuning_r;
            phiTuning_r = &tuning[phiIndx];

            if(covModel == "matern"){
                tuning[nuIndx] = *nuTuning_r;
                nuTuning_r = &tuning[nuIndx];
            }
		}
		//delete
        delete [] tuning;
        delete [] fixed;
        delete [] samples_r;
        delete [] accept_r;
        delete [] turning_r;
        delete [] paramsCurrent;
        delete [] accept;
        delete [] C;
        delete [] theta;
        delete [] vU;
        delete [] betamu;
}

