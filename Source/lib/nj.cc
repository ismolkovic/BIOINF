#include "nj.h"


/******
			Neighbour joining
******/

/*
	NJ constructor
	params: distance matrix, number of taxons
*/
NJ::NJ(DistanceMat *_mat, int _r) : dMat(_mat), r(_r) { }


/*
	NJ destructor
	delete distance matrix
*/
NJ::~NJ() { delete dMat; }	



/*
	NJ main method
	return phylogenetic tree
*/
pTree* NJ::start() {
	pTree *pt = new pTree(r);
	int i, j, u = 0;
	int n = r;

	// r-3 iterations
	for(int iter = 0; iter < n - 3; iter++) {
		dMat->calcAllSums();
		findMin(i, j);
		double _dfu = d_fu(i, j); // u -> f
		double _dgu = d_gu(i, j, _dfu); // u -> g
		u = d_uk(i, j);  // u -> k

		// add new node
		pt->addNode(i, j, _dfu, _dgu, dMat); 
		r--;
	}
	
	// last step
	pt->addLastNode(u, dMat);

	return pt;
}


/*
	Find the pair of taxa for which Q(i,j) has its lowest value
	params: reference to Qi, reference to Qj
*/
void NJ::findMin(int &Qi, int &Qj) {
	double Qmin = (r -2 ) * dMat->get(0,1);
	Qmin= Qmin - dMat->sumRow(0) - dMat->sumRow(1);
	Qi = 0; Qj = 1;

	for(int i = 0; i < r - 1; i++) {
		for(int j = i + 1; j < r; j++) {
			double Qij = (r-2) * dMat->get(i, j) - dMat->sumRow(i) - dMat->sumRow(j);
			if( Qij < Qmin) {
				Qi = i; Qj = j; Qmin = Qij;
			}   
		}
	}
}


/*	
	Calculate distance between f and new node u -> d(f,u)
	params: first taxon, second taxon
*/
double NJ::d_fu(int f, int g) {
	double val = dMat->get(f, g) * 0.5;
	val += (double)1/(2*(r-2)) * (dMat->sumRow(f) - dMat->sumRow(g));
	return val;
}


/*	
	Calculate distance between g and new node u   -> d(g,u)
	params: fist taxon, second taxon, distance(f,u)
*/
double NJ::d_gu(int f, int g, double d_fu) {
	return dMat->get(f, g) - d_fu;
}


/*
	Calculate distance of the other taxa to the new node  -> d(u,k)
	params: fist taxon, second taxon 
*/
int NJ::d_uk(int f, int g) {
	int u = f; // position of new node in distance matrix

	// d(u,k) //
	dMat->insert(u, u, 0);
	for(int k = 0; k < r; k++) {
		if(k == f || k == g) continue;
		double val = (double)0.5*(dMat->get(f, k) + dMat->get(g, k) - dMat->get(f,g));
		
		dMat->insert(u, k, val);
		dMat->insert(k, u, val);
	}

	// remove row and column g
	dMat->removeRow(g);
	dMat->removeCol(g);
	return u;
}
