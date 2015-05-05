#include <armadillo>
#include <iostream>
#include <vector>

using namespace arma;
using namespace std;

int main(){

	mat A = randu<mat>(5,5);
	A.print("A");
	cout << endl;
	A.col(0).print();
	cout << endl;

	rowvec B = A.row(0).subvec(0,3);
	B.print("B");
	rowvec C = ones<rowvec>(4);
	C.print("C");
	A.row(2).subvec(0,3) = C;
	A.print("A");


	return 0;
}
