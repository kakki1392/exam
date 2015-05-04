#include <armadillo>
#include <cstdio>
#include <vector>

using namespace arma;
using namespace std;

class Line{
	public:
		Line();
		~Line();
		size_t i1;
		size_t j1;
		size_t i2;
		size_t j2;
		int dir;
};

class Koch{
	public:
		Koch();
		~Koch();
		vector<Line> lines;
		size_t n;         //number of lines
		int l_max;        //maximum allowed l
		int l;            //current l
		double s;         //length of line
		double L;         //length of line at l=0
		int m;            //
		double delta_min; //minimum grid step
		double delta;     //used grid step
		double grid_max;

		size_t N;
		size_t origin;
		umat boundary;
		umat interior;
		mat grid;
		vec x;

		void fill_x();
		void initialize_line();
		size_t intpower(size_t base, size_t exponent);

};
