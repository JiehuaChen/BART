#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include "testSPLM.h"

//#ifdef __cplusplus
//extern "C"{
//#endif

#include "spLM.h"

//#ifdef __cplusplus
//}
//#endif

testSPLM::testSPLM()
{
//		Y_r = new double[DIMENSION];
//		X_r = new double[DIMENSION];
//		p_r = 0;
//		n_r = 0;
//		coordsD_r = new double[DIMENSION][DIMENSION]; 
		coordsD_r = new double[DIMENSION*DIMENSION];
//		betaPrior_r = "";
//		sigmaSqIG_r = new double[2];
//		tauSqIG_r = new double[2];
//		nuUnif_r = new double[2]; 
//		phiUnif_r = new double[2];
//		phiStarting_r = 0;
//		sigmaSqStarting_r = 0; 
//		tauSqStarting_r = 0;
//		nuStarting_r = 0;
//		phiTuning_r = 0; 
//		sigmaSqTuning_r = 0; 
//		tauSqTuning_r = 0; 
//		nuTuning_r = 0;
//		nugget_r = true; 
//		covModel_r = "";
//		amcmc_r = true; 
//		nBatch_r = 0; 
//		batchLength_r = 0;
//		acceptRate_r = 0; 
//		verbose_r = 0;
//		nReport_r = 0;
		result = new double*[3];
		result[0] = new double[DIMENSION];
		result[1] = new double[DIMENSION];
		result[2] = new double[DIMENSION];
		betaNorm_r = NULL;
		return;
}

testSPLM::~testSPLM()
{
		return;
}

int testSPLM::start(string inputfilename)
{
	GetTestData(inputfilename);
    spLM( Y_r, n_r,  coordsD_r, 
	      sigmaSqIG_r,  tauSqIG_r,  nuUnif_r, phiUnif_r,
		  phiStarting_r,  sigmaSqStarting_r, tauSqStarting_r,  nuStarting_r,
		  phiTuning_r, sigmaSqTuning_r,  tauSqTuning_r,  nuTuning_r, 
	      nugget_r,  covModel_r,  amcmc_r, nBatch_r,  batchLength_r,  acceptRate_r, 
	      verbose_r,  nReport_r, result);	
		return 0;
}

int testSPLM::GetTestData(string infile)
{
		string s_oneline;
		ifstream ifile;
//		int ret = 0;
		string param;
		int i = 0, j = 0;

		ifile.open(infile.c_str(), ios::in);
		if(ifile.fail())
		{
				printf("Cannot open file: %s\n",infile.c_str());
				exit(1);
		}
		getline(ifile, s_oneline);
		while(!ifile.eof())
		{
				//printf("Oneline: %s\n",s_oneline.c_str());
		    	string::size_type l_pos = 0, r_pos = 0;
		    	r_pos = s_oneline.find_first_of(':',l_pos);
		    	if(r_pos != string::npos)
		    	{
						param = s_oneline.substr(l_pos, r_pos-l_pos);
						if(0 == param.compare("Y_r")){
								i = 0;
								getline(ifile, s_oneline);
								while(!ifile.eof()){
                                        l_pos = 0; r_pos = 0;
										r_pos = s_oneline.find_first_of(':',l_pos);
										if(r_pos == string::npos) {
												Y_r[i] = atof(s_oneline.c_str());
												i++;
										}
										else{
												break;
										}
										getline(ifile, s_oneline);
								}
								printf("Y_r:\n");
								for(i = 0;i < DIMENSION; i++){
										printf("%0.12lf\n",Y_r[i]);
								}
								continue;
						}else if(0 == param.compare("n_r")){
								n_r = atoi(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("n_r: %d\n",n_r);
						}else if(0 == param.compare("coordsD_r")){
								i = 0;j = 0;
								getline(ifile, s_oneline);
								while(!ifile.eof()){
                                        l_pos = 0; r_pos = 0;
										r_pos = s_oneline.find_first_of(':',l_pos);
										if(r_pos == string::npos) {
												l_pos = 0; r_pos = 0; j = 0;
												r_pos = s_oneline.find_first_of(' ',l_pos);
												while(r_pos != string::npos){
														coordsD_r[i * DIMENSION + j] = atof(s_oneline.substr(l_pos, r_pos-l_pos).c_str());
														j++;
														l_pos = r_pos + 1;
														r_pos = s_oneline.find_first_of(' ',l_pos);
												}
												i++;
										}
										else{
												break;
										}
										getline(ifile, s_oneline);
								}
//								printf("coordsD_r:\n");
//								for(i = 0;i < DIMENSION; i++){
//										for(j = 0;j < DIMENSION; j++){
//												printf("%0.12lf ",coordsD_r[i * DIMENSION + j]);
//										}
//										printf("\n");
//								}
								continue;
						}else if(0 == param.compare("sigmaSqIG_ra")){
								sigmaSqIG_r[0] = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("sigmaSqIG_ra: %0.12lf\n", sigmaSqIG_r[0]);
						}else if(0 == param.compare("sigmaSqIG_rb")){
								sigmaSqIG_r[1] = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("sigmaSqIG_rb: %0.12lf\n", sigmaSqIG_r[1]);
						}else if(0 == param.compare("tauSqIG_ra")){
								tauSqIG_r[0] = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("tauSqIG_ra: %0.12lf\n", tauSqIG_r[0]);
						}else if(0 == param.compare("tauSqIG_rb")){
								tauSqIG_r[1] = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("tauSqIG_rb: %0.12lf\n", tauSqIG_r[1]);
						}else if(0 == param.compare("nuUnif_ra")){
								nuUnif_r[0] = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("nuUnif_ra: %0.12lf\n", nuUnif_r[0]);
						}else if(0 == param.compare("nuUnif_rb")){
								nuUnif_r[1] = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("nuUnif_rb: %0.12lf\n", nuUnif_r[1]);
						}else if(0 == param.compare("phiUnif_ra")){
								phiUnif_r[0] = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("phiUnif_ra: %0.12lf\n", phiUnif_r[0]);
						}else if(0 == param.compare("phiUnif_rb")){
								phiUnif_r[1] = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("phiUnif_rb: %0.12lf\n", phiUnif_r[1]);
						}else if(0 == param.compare("phiStarting_r")){
								phiStarting_r = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("phiStarting_r: %0.12lf\n", phiStarting_r);
						}else if(0 == param.compare("sigmaSqStarting_r")){
								sigmaSqStarting_r = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("sigmaSqStarting_r: %0.12lf\n", sigmaSqStarting_r);
						}else if(0 == param.compare("tauSqStarting_r")){
								tauSqStarting_r = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("tauSqStarting_r: %0.12lf\n", tauSqStarting_r);
						}else if(0 == param.compare("nuStarting_r")){
								nuStarting_r = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("nuStarting_r: %0.12lf\n", nuStarting_r);
						}else if(0 == param.compare("phiTuning_r")){
								phiTuning_r = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("phiTuning_r: %0.12lf\n", phiTuning_r);
						}else if(0 == param.compare("sigmaSqTuning_r")){
								sigmaSqTuning_r = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("sigmaSqTuning_r: %0.12lf\n",sigmaSqTuning_r );
						}else if(0 == param.compare("tauSqTuning_r")){
								tauSqTuning_r = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("tauSqTuning_r: %0.12lf\n",tauSqTuning_r );
						}else if(0 == param.compare("nuTuning_r")){
								nuTuning_r = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("nuTuning_r: %0.12lf\n", nuTuning_r);
						}else if(0 == param.compare("nugget_r")){
							    nugget_r = atoi(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("nugget_r: %d\n",nugget_r );
						}else if(0 == param.compare("covModel_r")){
							    covModel_r = s_oneline.substr(r_pos+2,s_oneline.size()).c_str();
								printf("covModel_r: %s\n",covModel_r.c_str() );
						}else if(0 == param.compare("amcmc_r")){
							    amcmc_r = atoi(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("amcmc_r: %d\n",amcmc_r );
						}else if(0 == param.compare("nBatch_r")){
							    nBatch_r = atoi(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("nBatch_r: %d\n",nBatch_r );
						}else if(0 == param.compare("batchLength_r")){
							    batchLength_r = atoi(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("batchLength_r: %d\n",batchLength_r );
						}else if(0 == param.compare("acceptRate_r")){
								acceptRate_r = atof(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("acceptRate_r: %0.12lf\n",acceptRate_r );
						}else if(0 == param.compare("verbose_r")){
							    verbose_r= atoi(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("verbose_r: %d\n",verbose_r );
						}else if(0 == param.compare("nReport_r")){
							    nReport_r = atoi(s_oneline.substr(r_pos+2,s_oneline.size()).c_str());
								printf("nReport_r: %d\n",nReport_r );
                        }
				}
				getline(ifile, s_oneline);
		}
		return 0;
}

int testSPLM::ParseOneInput(string &oneline)
{
	string data;
	double d_data;
	string::size_type l_pos = 0, r_pos = 0;

	r_pos = oneline.find_first_of(' ',l_pos+1);
	while(r_pos != string::npos)
	{
		data = oneline.substr(l_pos, r_pos-l_pos);
		d_data = atof(data.c_str());
		l_pos = r_pos + 1;
		r_pos = oneline.find_first_of(' ',l_pos+1);
	}
	return int(d_data);
}

int testSPLM::GetRangeData(string infile)
{
	std::vector<double>* p_v = new std::vector<double>;
	ifstream range_f(infile.c_str()); 	
	int numbers_in_row;
	double num;
	string str;
			
	while (!range_f.eof()) {
		getline(range_f,str); // read a line
		stringstream ss(str); // convert to string stream
		numbers_in_row = 0;
		while (ss >> num) {
			if (numbers_in_row > 1) {
				cerr << "Error: range file has a row with more than one number." << endl;
				continue;
			} 
			(*p_v).push_back(num); // store words from stream into vector
			numbers_in_row++;
		}
	}

	return 0;
}
