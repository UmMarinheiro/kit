#include <cstdio>
#include <iostream>

#include <limits>
#include <queue>
#include <vector>
#include <list>
#include <set>

using namespace std;

#include "data.h"
#include "hungarian.h"

#define MODE HUNGARIAN_MODE_MINIMIZE_COST
enum searchType {DFS, BFS, LOWEST};
searchType searchMethod = DFS;

double **cost;
Data * dt;
void initializeData(char* file)
{
    cout << "Reading " << file << " ..." << endl;
    dt = new Data(2, file);
	dt->readData();
    cout << "Succesfully read " << file << " !" << endl;
}
void freeData() {delete dt;}
void initializeCost()
{
	cost = new double*[dt->getDimension()];
	for (int i = 0; i < dt->getDimension(); i++){
		cost[i] = new double[dt->getDimension()];
		for (int j = 0; j < dt->getDimension(); j++){
			cost[i][j] = dt->getDistance(i,j);
		}
	}
}
void freeCost()
{
	for (int i = 0; i < dt->getDimension(); i++) delete [] cost[i];
	delete [] cost;
}

typedef struct Solution
{
    vector<int> sequence;
    double cost;
    void print() const
    {
        for(int i = 0; i < sequence.size() - 1; i++)
            cout << sequence.at(i) << " -> ";
        cout << sequence.back() << std::endl;

        cout << "Cost: " << cost << std::endl;
    }
    void print(string name) const {cout<<name<<"="<<endl; print(); }
}Solution;

typedef struct Node
{
	vector<pair<int,int>> forbidden_arcs;
	vector<vector<int>> subtours;

	double lower_bound;
	int chosen;
	bool feasible;

	bool operator<(const Node& other) const
    {
		return lower_bound < other.lower_bound;   
    } 
} Node;

int getAssingment(hungarian_problem_t *p, int i)
{
	for(int j = 1; j <= dt->getDimension(); j++) 
		if(p->assignment[i-1][j-1])return j;
	return -1;
}
vector<vector<int>> getSolutionHungarian(Node *node)
{
	hungarian_problem_t p;

	vector<int> forbiden_values(node->forbidden_arcs.size());
	for(int i = 0; i < node->forbidden_arcs.size(); i++)
	{
		pair<int,int> &current = node->forbidden_arcs[i];
		forbiden_values[i] = cost[current.first-1][current.second-1];

		cost[current.first-1][current.second-1] = INFINITE;
	}

	hungarian_init(&p, cost, dt->getDimension(), dt->getDimension(), MODE); 

	node->lower_bound = hungarian_solve(&p);

	set<int> nonVisited;
	for(int i = 1; i <= dt->getDimension(); i++) nonVisited.insert(i);

	int n = 0;
	int minIndex = 0, minSize = INFINITE;

	node->subtours = {};

	while(!nonVisited.empty())
	{
		int toVisit = *nonVisited.begin();
		node->subtours.push_back({});

		while(toVisit)
		{
			int current = toVisit;

			if(nonVisited.find(current) == nonVisited.end()) break;

			node->subtours[n].push_back(current);
			
			toVisit = getAssingment(&p, current);
			nonVisited.erase(current);
		}
		if(node->subtours[n].size() < minSize)
		{
			minSize = node->subtours[n].size();
			minIndex = n;
		}
		n++;
	}
	hungarian_free(&p);

	node->feasible = node->subtours.size() == 1;
	node->chosen = minIndex;

	for(int i = 0; i < node->forbidden_arcs.size(); i++)
	{
		pair<int,int> &current = node->forbidden_arcs[i];
		
		cost[current.first-1][current.second-1] = forbiden_values[i];
	}
	return node->subtours;
}
void addSubdivisionsOfNodeToList(Node *parent, list<Node> &tree)
{
	Node n;
	n.forbidden_arcs = parent->forbidden_arcs;

	pair<int,int> forbidden_arc = {
		parent->subtours[parent->chosen][parent->subtours[parent->chosen].size()-1], 
		parent->subtours[parent->chosen][0]
	};
	n.forbidden_arcs.push_back(forbidden_arc);
	
	tree.push_back(n);

	for (int i = 0; i < parent->subtours[parent->chosen].size() - 1; i++)
	{
		n.forbidden_arcs = parent->forbidden_arcs;

		forbidden_arc = {
			parent->subtours[parent->chosen][i], 
			parent->subtours[parent->chosen][i+1]
		};
		n.forbidden_arcs.push_back(forbidden_arc);
		
		tree.push_back(n);
	}
}

_List_iterator<Node> DFSFunction(list<Node> &tree) {return prev(tree.end());}
_List_iterator<Node> BFSFunction(list<Node> &tree) {return tree.begin();}

Solution BranchNBoundList()
{
	Solution best = {{}, numeric_limits<double>::infinity()};
	double &upper_bound = best.cost; 

	Node root; 
	list<Node> tree = {root};
	_List_iterator<Node> (*searchFunction)(list<Node> &tree) = searchMethod==DFS?&DFSFunction:&BFSFunction;

	while(!tree.empty())
	{
		auto node = searchFunction(tree);
		vector<vector<int>> subtours = getSolutionHungarian(&*node);

		if(node->lower_bound > upper_bound) 
		{
			tree.erase(node);
			continue;
		}

		if(node->feasible) 
		{
			if(node->lower_bound < upper_bound)
			{
				upper_bound = node->lower_bound;
				best.sequence = subtours[0];
				best.sequence.push_back(subtours[0][0]);
			}
		}
		else addSubdivisionsOfNodeToList(&*node, tree);
		tree.erase(node);
	}
	return best;
}

void addSubdivisionsOfNodeToPriorityQueue(Node *parent, priority_queue<Node> &tree)
{
	Node n;
	n.forbidden_arcs = parent->forbidden_arcs;

	pair<int,int> forbidden_arc = {
		parent->subtours[parent->chosen][parent->subtours[parent->chosen].size()-1], 
		parent->subtours[parent->chosen][0]
	};
	n.forbidden_arcs.push_back(forbidden_arc);
	
	tree.push(n);

	for (int i = 0; i < parent->subtours[parent->chosen].size() - 1; i++)
	{
		n.forbidden_arcs = parent->forbidden_arcs;

		forbidden_arc = {
			parent->subtours[parent->chosen][i], 
			parent->subtours[parent->chosen][i+1]
		};
		n.forbidden_arcs.push_back(forbidden_arc);
		
		tree.push(n);
	}
}
Solution BranchNBoundPriorityQueue()
{
	Solution best = {{}, numeric_limits<double>::infinity()}; // TODO Testar com construcao
	double &upper_bound = best.cost; 

	Node root; 
	priority_queue<Node> tree;

	tree.push(root);

	while(!tree.empty())
	{
		auto node = tree.top(); 
		tree.pop();
		vector<vector<int>> subtours = getSolutionHungarian(&node);

		if(node.lower_bound > upper_bound) continue;

		if(node.feasible) 
		{
			if(node.lower_bound < upper_bound)
			{
				upper_bound = node.lower_bound;
				best.sequence = subtours[0];
				best.sequence.push_back(subtours[0][0]);
			}
		}
		else addSubdivisionsOfNodeToPriorityQueue(&node, tree);
	}
	return best;
}
int main(int argc, char** argv) 
{
	initializeData(argv[1]);
	initializeCost();

	if(argc >= 3 && ( (string)argv[2]=="BFS" || (string)argv[2]=="LOWEST"))
	{
		if((string)argv[2]=="BFS") searchMethod = BFS;
		else searchMethod = LOWEST;
		cout<<"Searching method: "<<argv[2]<<endl;
	}
	else cout<<"Searching method: DFS"<<endl;

	cout << "Running BnB..." << endl;
    clock_t before = clock();
	Solution s = searchMethod==LOWEST?BranchNBoundPriorityQueue():BranchNBoundList();

    float duration = (clock()-before);
    s.print("BnB Solution");
    cout << "Took " << (float)duration/CLOCKS_PER_SEC << " seconds" << endl;
	
	freeData();
	freeCost();

	return 0;
}
