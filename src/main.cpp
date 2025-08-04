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

void updateCost(Solution &toUpdate)
{
    toUpdate.cost = 0;
    for(int i = 0; i < toUpdate.sequence.size()-1; i++) 
        toUpdate.cost += dt->getDistance(toUpdate.sequence[i], toUpdate.sequence[i+1]);
}

typedef struct InsertionInfo
{
    int insertedNode;
    int removedEdge;
    double cost;
} InsertionInfo;
vector<int> choseFromInterval(int n, int first, int limit)
{
    vector<int> chosen;
    for(int i = 0; i < n; i++)
    {
        int toAdd = (rand() % ((limit-first)-i)) + first;
        bool inserted = false;
        for(int j = 0; j < chosen.size(); j++) 
        {
            if(toAdd < chosen[j]) 
            {
                chosen.insert(chosen.begin() + j,toAdd);
                inserted = true;
                break;
            }
            toAdd++;
        }
        if(!inserted) chosen.push_back(toAdd);
    }
    return chosen;
}
vector<int> choseRandom3NodeSequence()
{
    vector<int> chosen = {1};
    vector<int> sorted3 = choseFromInterval(3, 2, dt->getDimension());

    int i = rand()%3;
    chosen.push_back(sorted3[i]);
    sorted3.erase(sorted3.begin() + i);

    i = rand()%2;
    chosen.push_back(sorted3[i]);
    chosen.push_back(sorted3[!i]);

    chosen.push_back(1);
    
    return chosen;
}
vector<int> getUnusedNodes(const vector<int> & used)
{
    vector<int> unused;
    
    for(int i = 1; i < dt->getDimension() + 1; i++)
    {
        bool isUsed = false;
        for(int node : used)
        {
            if(i == node) 
            {
                isUsed = true;
                break;
            }   
        }
        if(!isUsed) unused.push_back(i);
    }

    return unused;
}
InsertionInfo calculateBestInsertion
    (const vector<int> & sequence, const vector<int> & inserting)
{
    InsertionInfo best;
	best.cost = numeric_limits<double>::infinity();

    for(int a = 0; a < sequence.size() - 1; a++)
    {
        int predecessor = sequence[a];
        int successor = sequence[a+1];
        for(int inserted : inserting) 
        {
			double cost = 
                - dt->getDistance(predecessor-1, successor-1)
                + dt->getDistance(predecessor-1, inserted-1)
                + dt->getDistance(inserted-1, successor-1);

			if(cost < best.cost)
			{
				best.cost = cost;
				best.insertedNode = inserted;
				best.removedEdge = a;
			}
        }
    }
    return best;
}

Solution Construct()
{
    Solution s = Solution();
    
    s.sequence = choseRandom3NodeSequence();
    
    updateCost(s);
    std::vector<int> CL = getUnusedNodes(s.sequence);
    
    while(!CL.empty())
    {
        InsertionInfo insertion = 
            calculateBestInsertion(s.sequence, CL);
        
        s.sequence.insert(
            s.sequence.begin() + insertion.removedEdge + 1,
            insertion.insertedNode);
        s.cost += insertion.cost;

        CL = getUnusedNodes(s.sequence);
    }

    return s;
}

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
	Solution best = Construct();
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
