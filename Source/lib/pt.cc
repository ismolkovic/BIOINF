#include "pt.h"


/******
		PHYLOGENETIC TREE
******/

/*
	node constructor 
	params: node.id, node.distance to parent
*/
node::node(int _id, double _dist) : id(_id), distance(_dist) { } 


/*
	pTree constructor
	param: number of taxons
*/
pTree::pTree(int n) : curId(n), firstNode(n) {
	for(int i = 0; i < n; i++) {
		addLeaf(i, 0);
		mapNode.push_back(i);
	}
}


/*
	pTree desctructor
*/
pTree::~pTree() {
	for(int i = 0; i < (int)nodes.size(); i++)
		delete nodes.at(i);
}


/*	
	Add taxon to tree 
	prams: id - taxon, distane from parent
*/
void pTree::addLeaf(int id, double dist) {
	nodes.push_back(new node(id, dist));
}


/*
	Add new node to tree
	params: i = f, j = g, d(f,u), d(g,u), distance matrix
*/
void pTree::addNode(int i, int j, double dfu, double dgu, DistanceMat* mat) {
	int f = mapNode.at(i);
	int g = mapNode.at(j);

	node *n = new node(curId, 0);
	n->children.push_back(f); n->children.push_back(g);
	nodes.push_back(n);

	setDistance(f, dfu);
	setDistance(g, dgu);

	remapNodes(i, j, curId++);
	updateNodes(mat->getRow(i));
}


/*
	Set node distance from parent
	params: id node, distance
*/
void pTree::setDistance(int id, double dist) {
	nodes.at(id)->distance = dist;
}


/*
	Add child to node
	params: id node, id child
*/
void pTree::addChild(int id, int c) {
	nodes.at(id)->children.push_back(c);
}


/*
	Remap nodes
	mapNode store real node id, i,j are id's from current distance matrix
	param: first node id, second node i, real new node id
*/
void pTree::remapNodes(int i, int j, int u) {
	mapNode.at(i) = u;
	mapNode.erase(mapNode.begin() + j);
}


/*
	Update distance of other taxa to new node
	params: vector distances; distance=0 -> new node
*/
void pTree::updateNodes(vector<double>& duk) {
	for(int i = 0; i < (int)mapNode.size(); i++) {
		if(duk.at(i) == 0) continue; // it's a new node u 
		int index = mapNode.at(i);
		setDistance(index, duk.at(i));
	}
}


/*
	Add Last 3 nodes to tree
	params: id new node, distance matrix
*/
void pTree::addLastNode(int u, DistanceMat* mat) {	
	int c = 0, d = 0;
	if(u == 0) { c = 1; d = 2; }
	else if(u == 1 ) { c = 0; d = 2;}
	else { u = 0; d = 1;}

	// c - first node, d - second node, u - node id from last iteration

	double Dcd = mat->get(c,d); // d(c,d)
	double Dcu = mat->get(u,c); // d(c,u)
	double Ddu = mat->get(u,d); // d(d,u)

	// v - new node 
	double Dvu = (double)0.5*(Dcu + Ddu - Dcd); //d(v,u)
	double Dvd = (double)0.5*(Dcd + Ddu - Dcu); //d(v,d)
	double Dvc = (double)0.5*(Dcd + Dcu - Ddu); //d(v,c)

	// create new node, add children
	node *v = new node(curId, 0);
	v->children.push_back(mapNode.at(c));
	v->children.push_back(mapNode.at(d));
	v->children.push_back(mapNode.at(u));

	// update distances to last node v
	setDistance(mapNode.at(c), Dvc); // v -> c
	setDistance(mapNode.at(d), Dvd); // v -> d
	setDistance(mapNode.at(u), Dvu); // u - > v
	nodes.push_back(v);
}	


/*
	print phylogenetic tree
*/
void pTree::print() {
	for(unsigned int i = firstNode; i < nodes.size(); i++) {
		node *root = nodes.at(i);

		for(unsigned int j = 0; j < root->children.size(); j++) {
			node *child = nodes.at(root->children.at(j));
			printf("%d %d %.6lf\n", root->id, child->id, child->distance);
		}
	}

}

/*
	Print phygenetic tree using tree traversal
*/
void pTree::print2() {
	node *root = nodes.at(nodes.size()-1);
	traversal(root);
}


/*
	tree traversal 
	params: root - current node
*/
void pTree::traversal(node *root) {
	for(unsigned int i = 0; i < root->children.size(); i++) {
		node *nroot = nodes.at(root->children.at(i));
		//cout << root->id << " " << nroot->id << " " << nroot->distance << endl;
		printf("%d %d %.6lf\n", root->id, nroot->id, nroot->distance);

		if(nroot->children.size() > 0) { 
			traversal(nroot); 
		}
	}
}