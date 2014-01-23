#ifndef PHYLOGENETIC_TREE_HEADER
#define PHYLOGENETIC_TREE_HEADER 

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace std;


#include "dm.h"

/*
	class node 
*/

class node {
public:
	int id;
	double distance;
	vector< int > children;
	node(int, double);

};


/*
	class phylogenetic tree
*/

class pTree {
private:
	vector< node* > nodes;
	int curId;
	int firstNode;
	vector< int > mapNode;

	void setDistance(int, double);
	void addChild(int, int);
	void remapNodes(int, int, int);
	void updateNodes(vector<double>&);
	void traversal(node*);

public:
	pTree(int);
	~pTree();

	void addLeaf(int, double);
	void addNode(int, int, double, double, DistanceMat*);
	void addLastNode(int, DistanceMat*);
	void print();
	void print2();
};

 

#endif