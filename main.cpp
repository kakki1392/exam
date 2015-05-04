#include "koch.h"
#include <iostream>
#include <vector>


using namespace std;

int main(){

	Koch k;
	cout << k.lines.size() << endl;
	cout << k.grid_max << endl;
	cout << k.N << endl;
	k.fill_x();
	k.x.print("x: ");
	cout << endl << endl;
	//cout << k.x(k.origin) << endl;
	k.initialize_line();
	k.draw_lines();
	k.update_l();
	/*
	cout << k.lines.size() << endl;
	for(size_t i=0; i<k.lines.size(); i++){
		k.lines[i].print();
	}
	*/
	k.draw_lines();
	k.update_l();
	k.draw_lines();
	k.update_l();
	k.draw_lines();
	k.update_l();
	k.draw_lines();



	return 0;
}
