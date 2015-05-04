#include "koch.h"
#include <iostream>


using namespace std;

int main(){

	Line test;
	test.i1 = 1;
	test.i2 = 2;

	Koch k;
	cout << "delta " << k.delta << endl;
	cout << k.lines.size() << endl;
	cout << k.grid_max << endl;
	cout << k.N << endl;
	k.fill_x();
	k.x.print("x: ");
	cout << endl << endl;
	cout << k.x(k.origin) << endl;
	cout << k.intpower(2,4) << endl;
	k.initialize_line();


	return 0;
}
