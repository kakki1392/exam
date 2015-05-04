#include "koch.h"
#include <cmath>
#include <iostream>

Line::Line(){
	i1 = 0;
	i2 = 0;
	j1 = 0;
	j2 = 0;
	dir = 0;
}

Line::~Line(){}


Koch::Koch(){
	n = 4;
	l_max = 2;
	l = 0;
	L = 1.0;
	s = 1.0;
	m = 1.0;
	delta_min = pow(0.25,l_max);
	delta = delta_min;
	grid_max = L/2.0;
	for(int i = 1; i<(l_max+1);i++){
		grid_max = grid_max + L*pow(0.25,i);
	}
	lines = vector<Line>(n);
	double N_double = grid_max/delta;
	N = (size_t) N_double;
	N = 2*N + 3;
	origin = N/2;
	
	boundary = zeros<umat>(N,N);
	interior = zeros<umat>(N,N);
	grid = zeros<mat>(N,N);


	/*
	int factor = 4;
	int factor2 = 4;
	for(int i=1; i<l_max; i++){
		factor = factor*factor;
	}
	N = factor/2;
	for(int i=1; i<(l_max+1); i++){
		for(int j=1; j<i; j++){
			factor2 = factor2*factor2;
		}
		N = N + factor/factor2;
	}
	*/
	
}

Koch::~Koch(){}

void Koch::fill_x(){
	x = zeros<vec>(N);
	for(size_t i=(origin+1); i<N; i++){
		x(i) = x(i-1) + delta;
	}
	for(size_t i=0; i<(origin); i++){
		x(i) = -x(N-1-i);
	}
}

size_t Koch::intpower(size_t base, size_t exponent){
	size_t out = base;
	for(size_t i=0; i<(exponent-1); i++){
		out = out*base;
	}
	return out;
}

void Koch::initialize_line(){
	size_t steps = intpower(4,l_max)/4;
	cout << "steps: " << steps << endl;
	//Fill interior
	size_t longstep = 2*steps;
	cout << "longstep: " << longstep << endl;
	interior.submat(origin-longstep, origin-longstep, origin+longstep, origin+longstep) = ones<umat>(2*longstep+1,2*longstep+1);
	boundary = interior;
	boundary.submat(origin-longstep+1, origin-longstep+1, origin+longstep-1, origin+longstep-1) = zeros<umat>(2*(longstep-1)+1,2*(longstep-1)+1);
	interior.save("interior.dat", arma_ascii);
	boundary.save("boundary.dat", arma_ascii);

	//line 1
	
	



}


