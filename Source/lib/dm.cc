
#include "dm.h"

/******
			DistanceMat
******/


/*
	Distance matrix constructor
	param: n - matrix dimension (n x n matrix)
*/
DistanceMat::DistanceMat(int n) {
	d = vector< vector< double > >(n, vector<double>(n, 0)); 
	sumRowHash = vector< double >(n, 0);
}

/*
	Insert distance value between i and j
	params: i - first node, j - second node
*/
void DistanceMat::insert(int i , int j, double val) {
	if( i >= getSize() || j >= getSize() ) 
		throw string("Index out of range!!");
	d[ i ][ j ] = val;
}


/*
	Return distance value between i and j
	params: i - first node, j - second node
*/
double DistanceMat::get(int i, int j) {
	if( i >= getSize() || j >= getSize() ) 
		throw string("Index out of range!!");
	return d[ i ][ j ];	
}


/*
	Get sum of row
	param: i - row number
*/
double DistanceMat::sumRow(int i) {
	if( i >= getSize()) throw string("Index out of range!!");
	return sumRowHash.at(i);
}

/*
	Calculate sum of row
	param: i - row number
*/
double DistanceMat::calcSumRow(int i) {
	double sum = 0;
	for(int j = 0; j < getSize(); j++) {
		sum += d[ i ][ j ];
	}

	return sum;
}

/*
	Calculate sums of rows for current distance matrix
*/
void DistanceMat::calcAllSums() {
	for(int i = 0; i < getSize(); i++) {
		sumRowHash.at(i) = calcSumRow(i);
	}
}

/*
	Return number of rows/columns
*/
int DistanceMat::getSize() {
	return d.size();
}


/*
	Remove row
	param: row number
*/
void DistanceMat::removeRow(int row) {
	if(row >= getSize()) throw string("Index out of range!!");
	d.erase(d.begin()+row);
}


/*
	Remove Column
	param: column number
*/
void DistanceMat::removeCol(int col) {
	for(int i = 0; i < (int)d.size(); i++) {
		if( col >= (int)d.at(i).size()) 
			throw string("Index out of range!!");
		d.at(i).erase(d.at(i).begin() + col);
	}
}


/*
	Return row
	param: row number
*/
vector<double>& DistanceMat::getRow(int row) {
	if(row >= getSize()) throw string("Index out of range!!");
	return d.at(row);
}


/*
	print distance matrix
*/
void DistanceMat::print() {
	for(int i = 0; i < getSize(); i++) {
		for(int j = 0; j < i; j++) {
			cout << d[ i ][ j ] << " ";
		}
		cout << endl;
	}
}
