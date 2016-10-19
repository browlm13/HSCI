#include <functional>	//must compile with -std=c++11
#include <iostream>
#include <vector>
#include <math.h>       /* sin, pow, sqrt */

#include "newton.h"		/* newton */
#include "matrix.h"		/* Linspace, write */

using namespace std;

//residual function
double kepler(double t, double e, double x){return (e * sin(x)) -x -t;}

//derivitive of residual function
double dkepler(double e, double x){return (e * cos(x)) - 1;}

int main(){

	//settings for kepler's equation - kepler(t, e, x)
	double a = 2.0;
	double b = 1.25;
	double epsilon = sqrt(1 - (pow(b,2)/pow(a,2)));

	//settings for newton's method - newton(f, df, x, mi, tol, show)
	double xn = 0;
	int max_iter = 6;				//6
	double tolerance = pow(10,-5);	//-5
	bool show_iterate = false;

	//times {0, 0.001, . . . , 10}
	vector<double> ts = Linspace(0,10,10000);

	//vector of omegas corresponding to ts
	vector<double> ws;


	//run trials for all values of t
	double t;
	for (int i=0; i<ts.size(); i++){
		t = ts[i];

		//bind std::ref of t and e to the functions, use ref to pass by refrence
		const auto f = std::bind(&kepler, ref(t), ref(epsilon), std::placeholders::_1);
		const auto df = std::bind(&dkepler, ref(epsilon), std::placeholders::_1);

		//call newton
		ws.push_back(newton(f, df, xn, max_iter, tolerance, show_iterate));
		
		//the intial guess for each subsequent solve should be the returned value w
		xn = ws[i];
	}

	//x and y radial coordiantes coorispondint to omega/t
	vector<double> xs;
	vector<double> ys;

	//radial coordinates
	double rx;
	double ry;
	double r;
	for (int i=0; i<ws.size(); i++){
		r = (a*b) / (sqrt( pow(b * cos(ws[i]),2) + pow(a * sin(ws[i]),2) ));
		
		rx =  r * cos(ws[i]);
		xs.push_back(rx);

		ry = r * sin(ws[i]);
		ys.push_back(ry);
	}

	//write vectos ts, xs, ys, and ws
	write(ts, "t.txt");
	write(ws, "w.txt");
	write(xs, "x.txt");
	write(ys, "y.txt");

	return 0;
}
