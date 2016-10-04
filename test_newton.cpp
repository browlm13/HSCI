#include <cmath> /* pow */
#include <iostream>
//#include <vector>
#include "newton.h"

using namespace std;

//test function and its dervitive
double f(double x){return (pow(x,2)) * (x-3) * (x+2);}
double df(double x){return (4*pow(x,3)) - (3*pow(x,2)) - (12*x);}

int main(){

	//For your tests, start with initial guesses of x0 = {−3, 1, 2},
	double x0s[] = {-3,1,2};
	int x0s_size = 3;

	//tolerances of ε = {10^−1, 10^−5, 10−9}
	double tols[] = {-1,-5,-9};
	int tols_size = 3;

	//max iterations
	int max = 50;

	//call newton
	for (int i=0; i<x0s_size; i++){
		for (int j=0; j<tols_size; j++){
			double x0 = x0s[i];
			double epsilon = pow(10,tols[j]);
		
			newton(f, df, x0, max, epsilon, true);
		}
	}
	
	return 0;
}