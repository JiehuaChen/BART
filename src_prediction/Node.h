#ifndef GUARD_NODE_H
#define GUARD_NODE_H

#include <stdio.h>
#include <vector>
#include <cmath>

#include "List.h"

extern int NumX;
class Node {
public:
	Node();
	~Node();

	// next three are indicators for 
	// whether the node is a top, bottom, or nogrand
	int Top;
	int Bot;
	int Nog; // node with children but no grandchildren -- second-to-last to the end, so to speak
	double TrainedMu; //storing the Mu value of each bottom node trained by MCMC
	int Depth; //The depth of this node in this tree
	int Var; //The Index of the vector variable
	int OrdRule; //The Index of the split value in RuleMat
	double SplitVal; //Split value of this node

	// pointers for tree structure
	Node *Parent;
	Node *LeftC;
	Node *RightC;


	//functions
	bool GoRight(double *x);
        //void  currentFits(MuS* mod,int nTrain,double** xTrain,double* yTrain,int nTest,double** xTest,double* w, double** fits);
        double currentFits(unsigned int n, double* x);

};


#endif

