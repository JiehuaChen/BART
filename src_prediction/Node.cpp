#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

#include "Node.h"
#include "List.h"

typedef double *dp;


Node::Node()
{
	Top = 1;
	Bot = 1;
	Nog = 0;
	TrainedMu = 0.0;
}

Node::~Node()
{

}

bool Node::GoRight(double *x)
{
	if(x[Var] > SplitVal)
	{
		return true;
	}
	else
	{
		return false;
	}
}



double Node::currentFits(unsigned int n, double* x)
{
	double fittedMu = 0.0;
	if(Bot)
	{
		fittedMu = TrainedMu;
	}
	else
	{
		if(GoRight(x))
		{
			fittedMu = RightC->currentFits(n, x);
		}
		else
		{
			fittedMu = LeftC->currentFits(n, x);
		}
	}
	return fittedMu;
}



