/*
	create a basic matrix class
*/

/*

	Matrix Class

		functionality:
							Matrix(nrow, ncol)
							*Matrix(Matrix B)
							*Matrix(vector<double> v)
							void Vandermonde_Matrix(vector<double> v)
							operator=(Matrix)
							operator=(vector<double>)
							void transpose()
							Matrix dot(Matrix)
							vector<double> dot(vector<double>)			//should maybe be the same
							int write(string fname)

*/

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

//basic vector functions
std::vector<double> Linspace(double a, double b, int n);
std::vector<double> times(std::vector<double> a, std::vector<double> b);
std::vector<double> scale(std::vector<double> v, double s);
std::vector<double> add(std::vector<double> a, std::vector<double> b);
std::vector<double> sub(std::vector<double> a, std::vector<double> b);
double sum_vector(std::vector<double> v);
double norm(std::vector<double> v, int p);
std::vector<double> vector_pow(std::vector<double> v, double n);
std::vector<double> random(int n);
int max_index(std::vector<double> v);
void reverse(std::vector<double> & v);
void swap(std::vector<double> & v, int a, int b);	//swap row a with row b inplace
void print_vector(std::vector<double> v);
int write(std::vector<double> v, const char *outfile);

//basic matrix class
class Matrix{

	public:
		int nrow, ncol;
		std::vector<std::vector<double> > M;

		//constructors
		Matrix();
		Matrix(int r, int c);
		Matrix(std::vector<double>);

		//copy constructor
		Matrix(const Matrix & rhs);

		//assignment operators
		Matrix& operator=(const Matrix & rhs);
		Matrix& operator=(const std::vector<double> & rhs);

		//access operator
		double operator() (int i, int j);
		std::vector<double> operator() (int i);		//return col vector

		//destrutor
		//~Matrix();

		//							types of matricies

		//Identity 'constructor'
		void Identity_Matrix(int n);

		//Vandermond 'constructor'
		void Vandermonde_Matrix(std::vector<double> x, int c);

		//							basic matrix operations

		//clear
		void clear();

		//transpose
		void T();

		//dot (is this dot for real?)
		std::vector<double> dot(std::vector<double> v);
		Matrix dot(Matrix b);
		//Matrix dot(std::vector<double> v);
		//Matrix dot(Matrix B);


		//							elementary row operations

		//swap, swap row a with row b
		void swap_rows(int a, int b);

		//scale, scale row r by constant c
		void scale_row(int r, double c);

		//add, add row a to row b
		void add_row(int a, int b);

		//add row with factor, add row a to b with factor or c
		void add_scaled_row(int a, int b, double c);

		//eliminate leading non0 in row b with row a
		void eliminate_lead_non0(int a, int b);

		//eliminate leading non0 in row b with row a, returns scaler c
		double eliminate_0(int a, int b, int col);

		//							other

		//LUP factorization
		int LUP_factorization(Matrix &L, Matrix &U, Matrix &P);

		//linear Solver
		std::vector<double> solve(std::vector<double> b);

		//back substitution
		std::vector<double> back_substitution(std::vector<double> b);

		//forward substitution
		std::vector<double> forward_substitution(std::vector<double> b);


		//infinity norm (max sum of rows)
		double inifinity_norm();

		//column max index
		int col_max_index(int col);

		//write
		//int write(std::string fname);

		//print
		void print();

};

#endif