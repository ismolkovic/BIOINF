#ifndef DISTANCE_MATRIX
#define DISTANCE_MATRIX 

#include <vector>
#include <string>
#include <iostream>
#include <map>
using namespace std;


/*
	Distance Matrix
		- get, insert distance value
		- get row, get sum of row
		- calculate sum row, calc all sums
		- remove row, column
		- print distance matrix
*/

class DistanceMat {
private:
	vector< vector< double > > d;
	vector< double > sumRowHash;


public:
	DistanceMat(int);
	~DistanceMat() { }

	void insert(int, int, double);
	double get(int, int);
	int getSize();
	vector<double>& getRow(int);
	double sumRow(int);
	double calcSumRow(int);
	void calcAllSums();


	void removeRow(int);
	void removeCol(int);

	void print();
};


#endif