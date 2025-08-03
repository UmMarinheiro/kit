#include <cstdio>
#include <iostream>

#include <limits>
#include <vector>
#include <list>
#include <set>

using namespace std;

#include "data.h"
#include "hungarian.h"


#define MODE HUNGARIAN_MODE_MINIMIZE_COST

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

typedef struct
{
	vector<pair<int,int>> forbidden_arcs;
	vector<vector<int>> subtours;

	double lower_bound;
	int chosen;
	bool feasible;
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

	cout<<"Updating costMatrix..."<<endl;
	vector<int> forbiden_values(node->forbidden_arcs.size());
	cout<<"Blocking: ";
	for(int i = 0; i < node->forbidden_arcs.size(); i++)
	{
		pair<int,int> &current = node->forbidden_arcs[i];
		forbiden_values[i] = cost[current.first-1][current.second-1];

		cost[current.first-1][current.second-1] = INFINITE;
		cout<<current.first<<"-"<<current.second<<"  ";
	}
	cout<<endl;

	cout<<"Running hungarian..."<<endl;
	hungarian_init(&p, cost, dt->getDimension(), dt->getDimension(), MODE); // Carregando o problema

	node->lower_bound = hungarian_solve(&p);

	set<int> nonVisited;
	for(int i = 1; i <= dt->getDimension(); i++) nonVisited.insert(i);

	int n = 0;
	int minIndex = 0, minSize = INFINITE;

	cout<<"--Verifying subtours--"<<endl;
	node->subtours = {};

	while(!nonVisited.empty())
	{
		cout<<"Subtour "<<n<<": "<<endl;
		int toVisit = *nonVisited.begin();
		node->subtours.push_back({});

		cout<<"       ";
		while(toVisit)
		{
			int current = toVisit;

			if(nonVisited.find(current) == nonVisited.end()) break;

			node->subtours[n].push_back(current);
			cout<<current<<" ";

			toVisit = getAssingment(&p, current);
			nonVisited.erase(current);
		}
		cout<<endl;
		cout<<"Subtour "<<n<<" ended"<<endl;
		if(node->subtours[n].size() < minSize)
		{
			minSize = node->subtours[n].size();
			minIndex = n;
		}
		n++;
	}
	cout<<"--Subtours veryfied--"<<endl;
	hungarian_free(&p);

	cout<<"Hungarian concluded"<<endl;

	node->feasible = node->subtours.size() == 1;
	node->chosen = minIndex;

	cout<<"Cleaning costMatrix"<<endl;
	for(int i = 0; i < node->forbidden_arcs.size(); i++)
	{
		pair<int,int> &current = node->forbidden_arcs[i];
		
		cost[current.first-1][current.second-1] = forbiden_values[i];
		//cost[current.second-1][current.first-1] = forbiden_values[i];
	}
	cout<<"Solution complete"<<endl;
	return node->subtours;
}
double BranchNBound()
{
	double upper_bound = numeric_limits<double>::infinity(); // TODO Testar com construcao

	Node root; 
	list<Node> tree = {root};

	while(!tree.empty())
	{
		auto node = prev(tree.end());
		getSolutionHungarian(&*node);

		if(node->lower_bound > upper_bound) 
		{
			tree.erase(node);
			continue;
		}
		else if(node->feasible) upper_bound = min(upper_bound, node->lower_bound);
		else
		{
			Node n;
			n.forbidden_arcs = node->forbidden_arcs;

			pair<int,int> forbidden_arc = {
				node->subtours[node->chosen][node->subtours[node->chosen].size()-1], 
				node->subtours[node->chosen][0]
			};
			n.forbidden_arcs.push_back(forbidden_arc);
			
			tree.push_back(n);

			for (int i = 0; i < node->subtours[node->chosen].size() - 1; i++)
			{
				n.forbidden_arcs = node->forbidden_arcs;

				forbidden_arc = {
					node->subtours[node->chosen][i], 
					node->subtours[node->chosen][i+1]
				};
				n.forbidden_arcs.push_back(forbidden_arc);
				
				tree.push_back(n);
			}
		}
		tree.erase(node);
	}
	return upper_bound;
}
int main(int argc, char** argv) 
{
	initializeData(argv[1]);
	initializeCost();

	cout << "Running BnB..." << endl;
	double cost = BranchNBound();
	cout << "Finished Bnb with cost: " << cost << endl;
	
	freeData();
	freeCost();

	return 0;
}
