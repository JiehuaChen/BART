#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include "Utils.h"

std::vector< std::vector<double>* > * GetTestData(string infile)
{
	std::vector< std::vector<double>* > * p_v_v = NULL;
	std::vector<double>* p_v = NULL;
	string s_oneline;
	ifstream ifile;

	p_v_v = new std::vector< std::vector<double>* >;
	p_v_v->push_back(NULL);
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
                p_v = ParseOneInput(s_oneline);
                if(p_v != NULL)
                {
			p_v_v->push_back(p_v);
/*
			for(unsigned int i = 1;i<p_v->size();i++)
			{
				printf("%f ",p_v->at(i));
			}
			printf("\n");
*/
                }
                getline(ifile, s_oneline);
        }
	return p_v_v;
}

std::vector<double>* ParseOneInput(string &oneline)
{
	std::vector<double>* p_vd = NULL;
	p_vd = new std::vector<double>;
	
	string data;
	double d_data;
	string::size_type l_pos = 0, r_pos = 0;

	p_vd->push_back(0.0);

	r_pos = oneline.find_first_of(' ',l_pos+1);
	while(r_pos != string::npos)
	{
		data = oneline.substr(l_pos, r_pos-l_pos);
		d_data = atof(data.c_str());
		p_vd->push_back(d_data);
		l_pos = r_pos + 1;
		r_pos = oneline.find_first_of(' ',l_pos+1);
	}
	if(p_vd->size() > 1)
	{
		return p_vd;
	}
	else
	{
		delete p_vd;
		return NULL;
	}
}

std::vector<double>* GetRangeData(string infile)
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

	return p_v;
}