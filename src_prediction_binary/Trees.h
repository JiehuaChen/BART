#ifndef GUARD_TREES_H
#define GUARD_TREES_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

#include "List.h"
#include "Node.h"

using namespace std;
extern int NumX;

class Trees {
public:
	Trees();
	~Trees();

	ifstream ifile;
	ofstream ofile;
	std::vector<Node*> theTrees;
	std::vector<std::vector<Node*> *> NodesofTheTrees;
	int CurTree;
	
	Node* ParseOneLine(string &oneline);
	int ParseTreesFile(string filename);
	Node* Brother(Node *n);
	int BuildTrees();
	int BuildTree(int start, int end, Node* &proot, Node* parent);
	int FindNode(int start, int end, int dep);
	void PrintTrees();
	void PrintPre(Node* node);
	double GetPredictedResult(unsigned int n, double* x);
};

#endif

