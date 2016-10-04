#include <cmath> /* abs */
#include <functional> /* function type */
#include <iostream>

using namespace std;

double newton(function<double(double)> func, function<double(double)> dfunc, double x, int maxit, double tol, bool show_iterates){

	//display
	if(show_iterates){
		cout << "x0 = " << x << " , tolerence = " << tol << endl;
		cout << "iteration | current guess | residual" << endl;
	}

	//intial guess
	double xn = x;
	double xn_1 = x;
	for(int i=1; i <= maxit; i++)
	{

		//Newton Method Formula
		//xn_1 = xn - (func(xn))/(dfunc(xn));
		double h = func(xn)/dfunc(xn);
		xn_1 = xn - h;

		//check if current guess is within tolernce range
		//double residual = abs(xn-xn_1);
		double residual = abs(h);

		//display
		if(show_iterates)
			cout << endl << "\t" << i << "\t"<< xn_1 << "\t"<< residual << endl;

		//break condition
		if (residual < tol){
		  cout << "\n" << "success, root found at x = " << xn_1 << endl <<endl; 	//display root
		  //return 0.0;
		  return xn_1;
		}

		//update
		double tmp = xn_1;
		xn = tmp;
	}

	cerr << "error, no convergence." <<endl << endl;			//no convergence
	
	//return 0;
	return xn_1;	//if no converges should not return this value then***

}
