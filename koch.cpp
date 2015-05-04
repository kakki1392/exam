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

Line::Line(size_t i_1, size_t j_1, size_t i_2, size_t j_2, int d){
	i1 = i_1;
	j1 = j_1;
	i2 = i_2;
	j2 = j_2;
	dir = d;
}

Line::~Line(){}

void Line::print(){
	cout << "i1 :" << i1 << endl
	     << "j1 :" << j1 << endl
	     << "i2 :" << i2 << endl
	     << "j2 :" << j2 << endl
	     << "dir: " << dir << endl << endl;
}




Koch::Koch(){
	n = 4;
	l_max = 4;
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
	if(exponent==0){
		return base;
	}
	for(size_t i=0; i<(exponent-1); i++){
		out = out*base;
	}
	return out;
}

void Koch::initialize_line(){
	size_t steps = intpower(4,l_max)/4;
	s_index = 4*steps + 1;
	cout << "steps: " << steps << endl;
	cout << "s_index: " << s_index << endl;
	size_t longstep = 2*steps;
	/*
	stp = steps;
	lstp = longstep;
	cout << "longstep: " << longstep << endl;
	interior.submat(origin-longstep, origin-longstep, origin+longstep, origin+longstep) = ones<umat>(2*longstep+1,2*longstep+1);
	boundary = interior;
	boundary.submat(origin-longstep+1, origin-longstep+1, origin+longstep-1, origin+longstep-1) = zeros<umat>(2*(longstep-1)+1,2*(longstep-1)+1);
	interior.save("interior.dat", arma_ascii);
	boundary.save("boundary.dat", arma_ascii);
	*/
	
	//line 1
	Line line1;
	line1.i1 = origin - longstep;
	line1.i2 = origin + longstep;
	line1.j1 = origin + longstep;
	line1.j2 = line1.j1;
	line1.dir = 1;
	lines[0] = line1;
	
	//line 2
	Line line2;
	line2.i1 = origin + longstep;
	line2.i2 = line2.i1;
	line2.j1 = origin + longstep;
	line2.j2 = origin - longstep;
	line2.dir = 2;
	lines[1] = line2;
	
	//line 3
	Line line3;
	line3.i1 = origin + longstep;
	line3.i2 = origin - longstep;
	line3.j1 = origin - longstep;
	line3.j2 = line3.j1;
	line3.dir = 3;
	lines[2] = line3;
	
	//line 4
	Line line4;
	line4.i1 = origin - longstep;
	line4.i2 = line2.i1;
	line4.j1 = origin - longstep;
	line4.j2 = origin + longstep;
	line4.dir = 4;
	lines[3] = line4;
}


void Koch::draw_lines(){
	boundary = ones<umat>(N,N);
	uvec ins_col = zeros<uvec>(s_index);
	urowvec ins_row = zeros<urowvec>(s_index);
	cout << "s_index: " << s_index << endl;
	for(size_t it=0; it<lines.size(); it++){
		if(lines[it].dir == 1){
			boundary.col(lines[it].j1).subvec(lines[it].i1, lines[it].i2) = ins_col;
		}
		else if(lines[it].dir == 2){
			boundary.row(lines[it].i1).subvec(lines[it].j2, lines[it].j1) = ins_row;
		}
		else if(lines[it].dir == 3){
			boundary.col(lines[it].j1).subvec(lines[it].i2, lines[it].i1) = ins_col;
		}
		else{
			boundary.row(lines[it].i1).subvec(lines[it].j1, lines[it].j2) = ins_row;
		}
	}
	boundary(0,0) = 0;
	boundary(0,N-1) = 0;
	boundary.save("data/boundaryBigBig5.dat", raw_ascii);
}


void Koch::test_lines(){
	umat test = zeros<umat>(N,N);

	uvec del_col = zeros<uvec>(stp-1);
	uvec ins_col = ones<uvec>(2*lstp+1);
	urowvec del_row = zeros<urowvec>(stp-1);
	urowvec ins_row = ones<urowvec>(2*lstp+1);

	ins_col.print();
	cout << endl;
	ins_row.print();
	cout << endl;

	for(size_t it=0; it<lines.size(); it++){
		if(lines[it].dir == 1){
			test.col(lines[it].j1).subvec(lines[it].i1, lines[it].i2) = ins_col;
		}
		else if(lines[it].dir == 2){
			test.row(lines[it].i1).subvec(lines[it].j1, lines[it].j2) = ins_row;
		}
		else if(lines[it].dir == 3){
			test.col(lines[it].j1).subvec(lines[it].i2, lines[it].i1) = ins_col;
		}
		else{
			test.row(lines[it].i1).subvec(lines[it].j2, lines[it].j1) = ins_row;
		}
	}
	test.save("test_lines.dat",arma_ascii);
}

void Koch::update_l(){
	if((l+1) > l_max){
		return;
	}
	cout << "entering update_l()" << endl;
	int new_l = l+1;
	size_t new_n = 8*n;
	size_t s = intpower(4,l_max)/intpower(4,new_l); //CAREFUL
	size_t longstep = 2*s;
	cout << "new_l: " << new_l << endl;
	cout << "steps: " << s << endl;
	cout << "longstep: " << longstep << endl;
	vector<Line> new_lines;
	s_index = s + 1;
	l = new_l;
	n = new_n;
	cout << "s_index:" << s_index << endl;

	for(size_t it=0; it<lines.size(); it++){
		if(lines[it].dir == 1){ size_t i1 = lines[it].i1;
			size_t j1 = lines[it].j1;
			Line line1(i1,j1,i1+s,j1,1);
			new_lines.push_back(line1);
			Line line2(i1+s,j1,i1+s,j1+s,4);
			new_lines.push_back(line2);
			Line line3(i1+s,j1+s,i1+2*s,j1+s,1);
			new_lines.push_back(line3);
			Line line4(i1+2*s,j1+s,i1+2*s,j1,2);
			new_lines.push_back(line4);
			Line line5(i1+2*s,j1,i1+2*s,j1-s,2);
			new_lines.push_back(line5);
			Line line6(i1+2*s,j1-s,i1+3*s,j1-s,1);
			new_lines.push_back(line6);
			Line line7(i1+3*s,j1-s,i1+3*s,j1,4);
			new_lines.push_back(line7);
			Line line8(i1+3*s,j1,i1+4*s,j1,1);
			new_lines.push_back(line8);
		}
		else if(lines[it].dir == 2){
			size_t i1 = lines[it].i1;
			size_t j1 = lines[it].j1;
			Line line1(i1,j1,i1,j1-s,2);
			new_lines.push_back(line1);
			Line line2(i1,j1-s,i1+s,j1-s,1);
			new_lines.push_back(line2);
			Line line3(i1+s,j1-s,i1+s,j1-2*s,2);
			new_lines.push_back(line3);
			Line line4(i1+s,j1-2*s,i1,j1-2*s,3);
			new_lines.push_back(line4);
			Line line5(i1,j1-2*s,i1-s,j1-2*s,3);
			new_lines.push_back(line5);
			Line line6(i1-s,j1-2*s,i1-s,j1-3*s,2);
			new_lines.push_back(line6);
			Line line7(i1-s,j1-3*s,i1,j1-3*s,1);
			new_lines.push_back(line7);
			Line line8(i1,j1-3*s,i1,j1-4*s,2);
			new_lines.push_back(line8);
		}
		else if(lines[it].dir == 3){
			size_t i1 = lines[it].i1;
			size_t j1 = lines[it].j1;
			Line line1(i1,j1,i1-s,j1,3);
			new_lines.push_back(line1);
			Line line2(i1-s,j1,i1-s,j1-s,2);
			new_lines.push_back(line2);
			Line line3(i1-s,j1-s,i1-2*s,j1-s,3);
			new_lines.push_back(line3);
			Line line4(i1-2*s,j1-s,i1-2*s,j1,4);
			new_lines.push_back(line4);
			Line line5(i1-2*s,j1,i1-2*s,j1+s,4);
			new_lines.push_back(line5);
			Line line6(i1-2*s,j1+s,i1-3*s,j1+s,3);
			new_lines.push_back(line6);
			Line line7(i1-3*s,j1+s,i1-3*s,j1,2);
			new_lines.push_back(line7);
			Line line8(i1-3*s,j1,i1-4*s,j1,3);
			new_lines.push_back(line8);
		}
		else{
			size_t i1 = lines[it].i1;
			size_t j1 = lines[it].j1;
			Line line1(i1,j1,i1,j1+s,4);
			new_lines.push_back(line1);
			Line line2(i1,j1+s,i1-s,j1+s,3);
			new_lines.push_back(line2);
			Line line3(i1-s,j1+s,i1-s,j1+2*s,4);
			new_lines.push_back(line3);
			Line line4(i1-s,j1+2*s,i1,j1+2*s,1);
			new_lines.push_back(line4);
			Line line5(i1,j1+2*s,i1+s,j1+2*s,1);
			new_lines.push_back(line5);
			Line line6(i1+s,j1+2*s,i1+s,j1+3*s,4);
			new_lines.push_back(line6);
			Line line7(i1+s,j1+3*s,i1,j1+3*s,3);
			new_lines.push_back(line7);
			Line line8(i1,j1+3*s,i1,j1+4*s,4);
			new_lines.push_back(line8);
		}
	}
	lines = new_lines;
	cout << "leaving update_l()" << endl;
}









