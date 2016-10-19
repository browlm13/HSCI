#include "matrix.h"
#include <vector>
#include <stdlib.h> /* abs */
#include <cmath> /* pow */
#include <iostream>
#include <random>

using namespace std;

//main function
int main(){

	vector<double> ns = {5, 9, 17, 33, 65};
	vector<double> e_2norm;		//e = x - xhat
	vector<double> r_2norm; 	//r = Axhat -b, Ae = r

	vector<double> e;
	vector<double> r;
	int n;
	for (int i=0; i<ns.size(); i++){

		//create Vandermonde matrix of size nxn
		n = ns[i];
		vector<double> v = Linspace(0,1,n);

		Matrix M(n,n);
		M.Vandermonde_Matrix(v,n);

		//create random vector x
		vector<double> x = random(n);

		//solve for b
		vector<double> b = M.dot(x);


		//linear solve
		vector<double> computed_x = M.solve(b);

		//calculate error vector
		e = sub(x, computed_x);
		e_2norm.push_back(norm(e, 2));

		//calculate residual vector
		vector<double> r = M.dot(e);
		r_2norm.push_back(norm(r, 2));

	}

	//write vectors, 2-norm(residual) 2-norm(error)
	//write(e_2norm, "e_2norm.txt");
	//write(r_2norm, "r_2norm.txt");

	return 0;
}	