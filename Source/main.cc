/*
	Copyright (c) 2014. Igor Smolkovič 0036453288

	Fakultet elektrotehnike i računarstva
	Bioinformatika projekt
	
	Neighbor­ joining
*/



#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
using namespace std;

#include "nj.h"



int main(int argc, char const **argv)
{
	try {
		int N, a, b,ret; double d;
		ret = fscanf(stdin, "%d", &N); // number of taxons

		DistanceMat *dMat = new DistanceMat(N);

		//  read N(N-1)/2 distances between taxons
		for(int i = 0; i < N *(N-1)/2; ++i) {
			ret = fscanf(stdin, "%d %d %lf", &a, &b, &d);
			
			if( a >= N || b >= N)
				throw string("Index out of range!!");
			if (ret > 0) {
				dMat->insert(a, b, d);
				dMat->insert(b, a, d);
			}
		}

	
		NJ *nj = new NJ(dMat, N);
		pTree *pt = nj->start();

		pt->print(); // print phylogenetic tree

		delete nj; delete pt;
	}
	
	catch(string e) { 
		cout << e << endl;
	}

	return 0;
}
