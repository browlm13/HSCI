#include "matrix.h"
#include <vector>
#include <stdlib.h> /* abs */
#include <cmath> /* pow */
#include <iostream>
#include <random>

using namespace std;


//main function
int main(){

	//testing all matrix library functions

	//size of square matrices to be multiplied
	int n = 4;

	//create identity
	Matrix I;
	I.Identity_Matrix(n);

	cout << endl << "Identity: " << endl;
	I.print();

	//create vandermonde
	//vector<double> v = Linspace(n,1,n);		//works
	vector<double> v = Linspace(n,1,n);			//doesn't work

	//		*** row swaps make it fail ***

	Matrix V;
	V.Vandermonde_Matrix(v,n);

	cout << endl << "Vandermonde: " << endl;
	V.print();

	//test multiplication
	Matrix R;
	R = V.dot(I);

	cout << endl << "Result V*I: " << endl;
	R.print();

	R = I.dot(V);

	cout << endl << "Result I*V: " << endl;
	R.print();

	//test transpose
	cout << endl << "I: " << endl;
	I.print();

	cout << endl << "I Transpose: " << endl;
	I.T();
	I.print();

	cout << endl << "V: " << endl;
	V.print();

	cout << endl << "V Transpose: " << endl;
	V.T();
	V.print();

	//revert
	V.T();


	//test row operations
	cout << endl << "I: " << endl;
	I.print();

	//swap
	cout << endl << "I swap rows 0 and n-1: " << endl;
	I.swap_rows(0,n-1);
	I.print();

	//revert
	I.Identity_Matrix(n);

	//scale
	cout << endl << "I: " << endl;
	I.print();

	cout << endl << "I scale rows n/2  9: " << endl;
	int c = int(n/2);
	I.scale_row(c,9);
	I.print();


	//revert
	I.Identity_Matrix(n);

	cout << endl << "I: " << endl;
	I.print();

	cout << endl << "add row 0 to 1: " << endl;
	I.add_row(0,1);
	I.print();

	//revert
	I.Identity_Matrix(n);

	cout << endl << "I: " << endl;
	I.print();

	cout << endl << "add scaled row cn-1 to n/2, scaler = 7: " << endl;
	I.add_scaled_row(n-1,c,7);
	I.print();

	//revert
	I.Identity_Matrix(n);

	//vector matrix multiplication
	vector<double> b = Linspace(1,1,n);

	cout << endl << "vector multiplier b: " << endl;
	print_vector(b);

	cout << endl << "I: " << endl;
	I.print();

	cout << endl << "I * b: " << endl;
	R = I.dot(b);
	R.print();

	b = Linspace(1,1,n);

	cout << endl << "vector multiplier b: " << endl;
	print_vector(b);

	cout << endl << "V: " << endl;
	V.print();

	cout << endl << "V * b: " << endl;
	R = V.dot(b);
	R.print();

	b = Linspace(1,n,n);
	cout << endl << "vector multiplier b: " << endl;
	print_vector(b);

	cout << endl << "I: " << endl;
	I.print();

	cout << endl << "I * b: " << endl;
	R = I.dot(b);
	R.print();


	b = Linspace(1,n,n);
	cout << endl << "vector multiplier b: " << endl;
	print_vector(b);

	cout << endl << "V: " << endl;
	V.print();

	cout << endl << "V * b: " << endl;
	R = V.dot(b);
	R.print();

	//test term elimination (for gaussian elimination)
	cout << endl << "V: " << endl;
	V.print();

	cout << endl << "test term elimination (for gaussian elimination), eliminate row 1 with row n-1 at column 0, print scaler d " << endl;
	double d = V.eliminate_0(n-1,1,0);
	V.print();
	cout << endl << "d: " << d << endl;

	//test LU factorization (with P). test recombination LU = PA
	Matrix L;
	Matrix U;
	Matrix P;
	V.LUP_factorization(L,U,P);

	R = P.dot(V);

	cout << endl << "Result P*V: " << endl;
	R.print();

	R = L.dot(U);

	cout << endl << "Result L*U: " << endl;
	R.print();

	//testing forward and back substitution


	//testing solver

	//create random vector x
	vector<double> x = random(n);

	cout << " vector x (multiplier) :" << endl;

	print_vector(x);

	cout << endl << " ---- " << endl;

	//Matrix b = V.dot(x);
	b = V.dot(x);


	cout << " vector b (answer) :" << endl;

	print_vector(b);

	cout << endl << " ---- " << endl;

	//testing linear solve
	V.solve(b);


	//compare with actual answer
	cout << endl << "original vector x :" << endl;
	print_vector(x);


	return 0;
}	

	
	/*
	//create Vandermonde matrix of size nxn
	int n = 5;
	vector<double> v = Linspace(5,1,n);

	Matrix M(n,n);
	M.Vandermonde_Matrix(v,n);
	M.print();

	cout << endl << " ---- " << endl;

	//create random vector x
	vector<double> x = random(n);

	cout << " vector x (multiplier) :" << endl;

	print_vector(x);

	cout << endl << " ---- " << endl;

	//Matrix b = M.dot(x);
	vector<double> b = M.dot(x);


	cout << " vector b (answer) :" << endl;

	print_vector(b);

	cout << endl << " ---- " << endl;

	//testing linear solve
	M.solve(b);


	//compare with actual answer
	cout << endl << "original vector x :" << endl;
	print_vector(x);


	//testing solve function on upper triangular matrix
	cout << endl << "Matrix::solve (test), testing solve on upper triangular matrix :" << endl;

	Matrix test_u(3,3);
	test_u.M[0][0] = 1;
	test_u.M[0][1] = 0;
	test_u.M[1][0] = 1;
	test_u.M[1][1] = 1;
	test_u.M[0][2] = 0;
	test_u.M[2][0] = 1;
	test_u.M[2][1] = 1;
	test_u.M[2][2] = 1;
	test_u.print();

	vector<double> test_b;
	test_b.push_back(1);
	test_b.push_back(2);
	test_b.push_back(5);
	print_vector(test_b);

	vector<double> test_x = test_u.solve(test_b);
	*/