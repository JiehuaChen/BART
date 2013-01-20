#include "Trees.h"

Trees::Trees()
{
	CurTree = 0;
	NodesofTheTrees.push_back(new std::vector<Node*>);
	theTrees.push_back(NULL);
}

Trees::~Trees()
{
	for(unsigned int i = 0; i < NodesofTheTrees.size(); i++)
	{
		for(unsigned int j = 0; j < NodesofTheTrees[i]->size(); j++)
		{
			if(NodesofTheTrees[i]->at(j) != NULL)
			{
				delete NodesofTheTrees[i]->at(j);
			}
		}
		delete NodesofTheTrees[i];
	}
	theTrees.clear();
}

Node* Trees::ParseOneLine(string &oneline)
{
	Node* p_node = NULL;
	int Top = 0, Bot = 0, Nog = 0;
	int Depth = 0, Var = 0;
//	int OrdRule = 0;
	double TrainedMu = 0.0, SplitVal = 0.0;

	string data;
	string::size_type r_pos = 0, l_pos = 0;

	l_pos = oneline.find("Tree");
	if(l_pos != string::npos)
	{
		CurTree++;
		NodesofTheTrees.push_back(new std::vector<Node*>);
		NodesofTheTrees[CurTree]->push_back(NULL);
		return NULL;
	}

	l_pos = oneline.find("Depth:");
	if(l_pos == string::npos)
	{
		return NULL;
	}
	else
	{
		l_pos += strlen("Depth:");
		r_pos = oneline.find_first_of(' ',l_pos);
		data = oneline.substr(l_pos, r_pos-l_pos);
		Depth = atoi(data.c_str());
		
		l_pos = oneline.find("TBN: ",r_pos);
		l_pos += strlen("TBN: ");
		Top = oneline.at(l_pos) - '0';
		Bot = oneline.at(l_pos+1) - '0';
		Nog = oneline.at(l_pos+2) - '0';

		l_pos = oneline.find("Var:");
		l_pos += strlen("Var:");
		r_pos = oneline.find_first_of(' ',l_pos);
		data = oneline.substr(l_pos, r_pos-l_pos);
		Var = atoi(data.c_str());

		l_pos = oneline.find("ORDRule:",r_pos);
		if(l_pos == string::npos)
		{
			l_pos = oneline.find("Mu:",r_pos);
			if(l_pos == string::npos)
			{
				printf("Mu Error when Parsing Tree %d: Depth=%d\n",CurTree,Depth);
				printf("Oneline: %s\n",oneline.c_str());
				exit(1);
			}
			l_pos += strlen("Mu:");
			data = oneline.substr(l_pos, oneline.length()-l_pos);
			TrainedMu = atof(data.c_str());
		}
		else
		{
			l_pos = oneline.find_first_of('=',l_pos);
			if(l_pos == string::npos)
			{
				printf("ORDRule Error when Parsing Tree %d: Depth=%d\n",CurTree,Depth);
				printf("Oneline: %s\n",oneline.c_str());
				exit(1);
			}
			l_pos++;
			data = oneline.substr(l_pos, oneline.length()-l_pos);
			SplitVal = atof(data.c_str());
		}
		p_node = new Node();
		p_node->Top = Top;
		p_node->Bot = Bot;
		p_node->Nog = Nog;
		p_node->Depth = Depth;
		p_node->Var = Var;
		p_node->TrainedMu = TrainedMu;
		p_node->SplitVal = SplitVal;
		p_node->LeftC = NULL;
		p_node->RightC = NULL;
		p_node->Parent = NULL;
		//printf("Depth=%d,TBN:%d%d%d,Var:%d,Mu:%f,SplitVal:%f\n",Depth,Top,Bot,Nog,Var,TrainedMu,SplitVal);
		return p_node;
	}
}
int Trees::ParseTreesFile(string filename)
{
	Node* p_node;
	string s_oneline;
	ifile.open(filename.c_str(), ios::in);
	if(ifile.fail())
	{
		printf("Can not open file: %s\n",filename.c_str());
		exit(1);
	}
	getline(ifile, s_oneline);
	while(!ifile.eof())
	{
		//printf("Oneline: %s\n",s_oneline.c_str());
		p_node = ParseOneLine(s_oneline);
		if(p_node != NULL)
		{
			NodesofTheTrees[CurTree]->push_back(p_node);
		}
		getline(ifile, s_oneline);
	}
	return 0;
}

Node* Trees::Brother(Node *n)
{
        if(n->Top) return 0;

        if(n==((n->Parent)->LeftC)) {
                return (n->Parent)->RightC;
        } else {
                return (n->Parent)->LeftC;
        }
}

int Trees::BuildTrees()
{
	CurTree = 1;
	for(unsigned int i = 1; i < NodesofTheTrees.size();i++)
	{
		theTrees.push_back(NULL);
		BuildTree(1, int(NodesofTheTrees[i]->size() - 1),theTrees.at(CurTree), NULL);
		CurTree++;
	}
	return 0;
}

int Trees::BuildTree(int start, int end, Node* &proot, Node* parent)
{
	int tmp = 0;
	proot = NodesofTheTrees[CurTree]->at(start);
	proot->Parent = parent;
	tmp = FindNode(start+2,end,proot->Depth+1);
	if(tmp != -1)
	{
		BuildTree(start+1, tmp-1, proot->LeftC, proot);
		BuildTree(tmp, end, proot->RightC, proot);
	}	
	return 0;
}

int Trees::FindNode(int start, int end, int dep)
{
	int i = 0;
	if(start > end)
	{
		return -1;
	}
	for(i=start; i<=end; i++)
	{
		if(NodesofTheTrees[CurTree]->at(i)->Depth == dep)
		{
			return i;
		}
	}
	return -1;
}

void Trees::PrintTrees()
{
	for(unsigned int i = 1;i < theTrees.size(); i++)
	{
		printf("Tree%u\n",i);
		PrintPre(theTrees[i]);
	}
}

void Trees::PrintPre(Node * node)
{
	if(node == NULL)
	{
		return;
	}
	printf("Depth=%d,TBN:%d%d%d,Var:%d,Mu:%0.17f,SplitVal:%0.17f\n",
		node->Depth,node->Top,node->Bot,node->Nog,node->Var,node->TrainedMu,node->SplitVal);
	PrintPre(node->LeftC);
	PrintPre(node->RightC);

}

double Trees::GetPredictedResult(unsigned int n, double* x)
{
	double result = 0.0;
	for(unsigned int i = 1; i < theTrees.size();i++)
	{
		result += theTrees[i]->currentFits(n, x);
	}
	return result;
}