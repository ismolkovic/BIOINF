#ifndef NJ_HEADER
#define NJ_HEADER

#include <vector>
#include <string>
#include <iostream>
#include <map>
using namespace std;

#include "dm.h"
#include "pt.h"


/*	
	Neighbor joining - main class
		- find pair of taxa with lowest Q value
		- calculate distances fu, gu, uk
*/

class NJ {	
private:
	DistanceMat *dMat;
	int r;

	void findMin(int&, int&);
	double d_fu(int, int);
	double d_gu(int, int, double);
	int d_uk(int i, int j);

public:
	NJ(DistanceMat*, int);
	~NJ();
	pTree* start();
};

#endif