#include "koch.h"
#include <iostream>
#include <vector>
#include <unistd.h>
using namespace arma;
using namespace std;

int main(){
	wall_clock timer;
	timer.tic();

	Koch k;
	k.gplt.cmd("set term pdfcairo");
	k.gplt.cmd("set output 'data/koch3.pdf'");
	k.gplt.cmd("set size square");
	k.gplt.cmd("set xlabel 'x'");
	k.gplt.cmd("set ylabel 'y'");
	k.gplt.cmd("set title 'l = 3, L = 1'");


	cout << "time spent creating object: " << timer.toc() << endl;
	cout << k.lines.size() << endl;
	cout << k.grid_max << endl;
	cout << k.N << endl;
	timer.tic();
	k.fill_x();
	cout << "time spent filling grid: " << timer.toc();
	//k.x.print("x: ");
	cout << endl << endl;
	//cout << k.x(k.origin) << endl;
	k.initialize_line();
	//k.plot_lines();
	//k.draw_lines();
	//k.plot_boundary();
	k.update_l();
	//k.plot_lines();
	//k.draw_lines();
	//sleep(4);
	//k.plot_boundary();
	k.update_l();
	//k.plot_lines();
	//k.draw_lines();
	//sleep(4);
	//k.plot_boundary();
	k.update_l();
	k.plot_lines();
	k.gplt.cmd("unset output");
	/*
	k.update_l();
	sleep(2);
	k.plot_lines();
	k.update_l();
	sleep(2);
	k.plot_lines();
	*/
	//k.draw_lines();
	//sleep(4);
	//k.plot_lines();
	//k.plot_single_line_corners(3);
	/*
	k.plot_lines();
//	k.draw_lines();
	k.update_l();
	sleep(4);
	k.plot_lines();
	k.update_l();
	sleep(4);
	k.plot_lines();
	sleep(4);
	k.update_l();
	k.plot_lines();
	cout << k.lines.size() << endl;
	for(size_t i=0; i<k.lines.size(); i++){
		k.lines[i].print();
	}
	*/



	return 0;
}
