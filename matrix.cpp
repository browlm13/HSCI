#include "matrix.h"
//#include <vector>
#include <cmath> /* pow */
#include <algorithm> /* max */
#include <random>
#include <iostream>

// singularity tolerance
#define STOL 1.e-15

using namespace std;

/*
	Basic Vector Functions
*/

//create a new vector of linearly space data
vector<double> Linspace(double a, double b, int n){
	if (n<2) cerr << "Linspace::length must be > 1\n";
	vector<double> v(n);
	
	double h = (b-a)/(n-1);
	for (int i=0; i<n; i++)
	 	v[i] = a + i*h;

	return v;
}

//entry wise multiplication of vectors
vector<double> times(vector<double> a, vector<double> b){
	if (a.size() != b.size()) cerr << "times::two vectors must be of same length\n";
	vector<double> c;
	for (int i=0; i<a.size(); i++)
		c.push_back(a[i] * b[i]);
	return c;
}

//matrix scaler
vector<double> scale(vector<double> v, double s){
	vector<double> result;
	for (int i=0; i<v.size(); i++)
		result.push_back(v[i] * s);
	return result;
}

//add vectors element wise
std::vector<double> add(std::vector<double> a, std::vector<double> b){
	if (a.size() != b.size()) cerr << "times::two vectors must be of same length\n";
	vector<double> c;
	for (int i=0; i<a.size(); i++)
		c.push_back(a[i] + b[i]);
	return c;
}

//colapse vector to double, (sum entries)
double sum_vector(vector<double> v){
	double  sum = 0;
	for (int i=0; i<v.size(); i++)
		sum += v[i];
	return sum;
}

//raise each element to a power n in a vector
vector<double> vector_pow(vector<double> v, double n){
	vector<double> result;
	for (int i=0; i<v.size(); i++)
		result.push_back(pow(v[i],n));
	return result;
}

//generate random vector
vector<double> random(int n){
	//uniform distribution
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0,100);			//uniform real double distribution
	vector<double> v;
	for (int i=0; i<n; i++)
		v.push_back(distribution(generator));  							//generates number in the range(0...100)
	return v;
}

//returns max element index
int max_index(vector<double> v){
	int mx = 0;
	for (int i=1; i<v.size(); i++){
		if (v[i] > v[mx])
			mx = i;
	}
	return mx;
}

//reverses vector in place
void reverse(std::vector<double> & v){
	vector<double> tmp;

	for (int i=v.size()-1; i>=0; i--)
		tmp.push_back(v[i]);

	v = tmp;
}

//swap vector elements in place
void swap(vector<double> & v, int a, int b){
	double tmp = v[a];
	v[a] = v[b];
	v[b] = tmp;
}

//print vector
void print_vector(vector<double> v){
	for (int i=0; i<v.size(); i++)
		cout << v[i] << " ";
}

//write a vector to a file
int write(std::vector<double> v, const char *outfile) { 

  // open output file
  FILE *fptr = NULL;
  fptr = fopen(outfile, "w");
  if (fptr == NULL) {
    cerr << "Write:: error, unable to open " << outfile << " for writing\n";
    return 1;
  }

  // print data to file
  for (size_t i=0; i<v.size(); i++) {
      fprintf(fptr, "  %.16g", v[i]);
  }

  // close output file and return
  fclose(fptr);
  return 0;
}


/*

	Matrix Class Member Functions

*/

//constructors
Matrix::Matrix(){
	ncol = 1;
	nrow = 1;
	vector<double> tmp;
	M.push_back(tmp);
}

Matrix::Matrix(int r, int c){
	ncol = c;
	nrow = r;

	for (int i=0; i<c; i++){
		vector<double> tmp(r);
		M.push_back(tmp);
	}
}

Matrix::Matrix(vector<double> v){
	ncol = 1;
	nrow = v.size();
	M.push_back(v);
}

Matrix::Matrix(const Matrix & rhs){
	ncol = rhs.ncol;
	nrow = rhs.nrow;
	M = rhs.M;
}


//assigment operators
Matrix& Matrix::operator= (const Matrix & rhs){

    if (this != &rhs) // protect against invalid self-assignment
    {
    	nrow = rhs.nrow;
    	ncol = rhs.ncol;
    	M = rhs.M;
    }
    // by convention, always return *this
    return *this;
}

Matrix& Matrix::operator= (const vector<double> & rhs){
	nrow = rhs.size();
	ncol = 1;
	M[0] = rhs;

	return *this;
}

//access operators
double Matrix::operator() (int i, int j){return M[i][j];}

//return column vector
vector<double> Matrix::operator() (int i){return M[i];}

//make it an Identity of col/row size n
void Matrix::Identity_Matrix(int n){

	//create matrix of correct dimensions
	Matrix I(n,n);
	M = I.M;
	nrow = n;
	ncol = n;

	//fill diagnols with '1''s
	for(int i=0; i<ncol; i++)
		M[i][i] = 1;
}

//make it a Vandermonde
void Matrix::Vandermonde_Matrix(vector<double> v, int c){

	//clear
	M.clear();

	nrow = v.size();
	ncol = c;

	//create a Vandermonde Matrix from first column vector v
	M.push_back(v);

	//create all remaining columns
	for (int i=2; i<ncol+2; i++){
		vector<double> p = vector_pow(v, i);
		M.push_back(p);
	}
}

//clear matrix
void Matrix::clear(){
	vector<vector<double> > clr;
	M = clr;
	nrow = 0;
	ncol = 0;
}

//transpose of matrix
void Matrix::T(){
	Matrix R;
	R.clear();

	for (int r=0; r<nrow; r++)
	{
		vector<double> new_col;
		for(int c=0; c<ncol; c++){
			new_col.push_back(M[c][r]);
		}
		R.M.push_back(new_col);
	}

	int tmp = ncol;
	ncol = nrow;
	nrow = tmp;
	M = R.M;
}


//vector matrix multiplication
//Matrix Matrix::dot(vector<double> v){
vector<double> Matrix::dot(vector<double> v){

	//for each entry in the vector multiply it by corresponding column in matrix
	Matrix R(nrow, v.size());

	//scale matrix by vector entries
	for (int i=0; i<v.size(); i++){
		R.M[i] = scale(M[i], v[i]);
	}

	//collapse standing matrix (sum columns)
	//1 transpose matrix
	R.T();

	//collanpse column vectors
	vector<double> tmp;
	for(int i=0; i<ncol; i++){
		tmp.push_back(sum_vector(R.M[i]));
	}

	return tmp;
}

//Matrix Matrix multiplication
Matrix Matrix::dot(Matrix B){

	//exit condition
	if (this->ncol != B.nrow) cerr << "Matrix::dot two Matrices must be of correct dim.\n";
	
	Matrix C(this->nrow, B.ncol);

	//perform Matrix Matrix multiplication
	//perform vector matrix multiplication for each row of B (T() col), storing results in C
	//B.T();
	for (int i=0; i<B.ncol; i++)
		C.M[i] = this->dot(B.M[i]);

	return C;
}

/*
	elementary row operations
*/
//swap, swap row a with row b
void Matrix::swap_rows(int a, int b){
	if (a >= nrow || b >= nrow || a < 0 || b < 0)cerr << "error: Matrix::swap_rows - Indicies of matrix must be valid.\n";
	this->T();
	vector<double> tmp = M[a];
	M[a] = M[b];
	M[b] = tmp;
	this->T();
}

//scale, scale row r by constant c
void Matrix::scale_row(int r, double c){
	this->T();
	M[r] = scale(M[r],c);
	this->T();
}

//add, add row a to row b
void Matrix::add_row(int a, int b){
	this->T();
	M[b] = add(M[a], M[b]);
	this->T();
}

//add row with factor, add row a to b with factor or c
void Matrix::add_scaled_row(int a, int b, double c){
	this->T();
	M[b] = add( scale(M[a], c), M[b] );
	this->T();
}

//eliminate leading non zero value in row b with some scaled value of row a
void Matrix::eliminate_lead_non0(int a, int b){

	this->T();

	//find leading non0 term in row b
	int col = -1;
	for (int i=0; i< M[b].size(); i++){
		if (M[b][i] != 0){
			col = i;
			break;
		}
	}
	if (col != -1){
		double c = M[b][col] / M[a][col];
		this->T();
		this->add_scaled_row(a,b,-c);
		this->T();
	}
	this->T();
}

//eliminate col value in row b with some scaled value of row a, returns scaler c
double Matrix::eliminate_0(int a, int b, int col){

		if(M[col][b] != 0 && M[col][a] != 0){
			double c = -(M[col][b] / M[col][a]);
			this->add_scaled_row(a,b,c);
			return c;
		}
		
		return 1;
}

// LU and P Factorization, 
int Matrix::LUP_factorization(Matrix &L, Matrix &U, Matrix &P){
	if (this->nrow != this->ncol)	cerr << "(error) Matrix::LUP_factorization::Matrix mustbe square\n";

	//lower triangluar matrix base (set as identity matrix)
	L.Identity_Matrix(this->nrow);


	//upper triangular Matrix base (set this Matrix as starting place)
	U.nrow = this->nrow;
	U.ncol = this->ncol;
	U.M = this->M;

	// P matrix to keep track of row swaps (and scales?) (set identity matrix)
	P.Identity_Matrix(this->nrow);

	// perform Gaussian elimination for LU factorization
	for (int k=0; k<U.nrow -1; k++) {   // loop over diagonals

		//The error comes from row swapping 
		
		//PIVOT
		// find the pivot row p (max of col) //starting from row k down
		/*
		int p=k;
		for (int i=k; i<U.nrow; i++){
		  if( abs(U.M[k][i]) >  abs(U.M[k][p]) )
			p=i;
		}*/

		int p=k;
		for (int i=k+1; i<U.nrow; i++){
		  //if( abs(U.M[i][k]) >  abs(U.M[p][k]) )
			if( abs(U.M[k][i]) >  abs(U.M[k][p]) )
				p=i;
		}

		//swap L U and P rows around the pivot
		//U
		//L.swap_rows(p,k);		//for some matricies this makes it better and for others worse
		//U.swap_rows(p,k);

		//P
		//P.swap_rows(p,k);
		

		//ELIMINATION
		//Perform eliminations using row k
		double c;	//scaler
		for (int i=k; i<U.nrow-1; i++){
			c = U.eliminate_0(k, i+1, k);

			if (c == 1)
				L.M[k][i+1] =  0;
			else
				L.M[k][i+1] = -c;
		}

	}

	//1.  Construct the matrices  L, U and P.
	cout << endl << "U, L, P:" << endl;
	U.print();
	cout << endl;
	L.print();
	cout << endl;
	P.print();
	cout << endl;

	return 0;
}


vector<double> Matrix::solve(vector<double> b){
	if (this->nrow != b.size())	cerr << "Matrix::solve::Matrix A and vector b must be of apropriate dim\n";
	if (this->nrow != this->ncol) cerr << "Matrix::solve::Matrix must be square\n";
	
	vector<double> x;

	/*
		Theorem (PA = LU;  Factorization with Pivoting).  Given that A is nonsingular.  The solution x to the linear system  Ax = b, is found in four steps:  

	    1.  Construct the matrices  L, U and P.  
	    2.  Compute the column vector  Pb.
	    3.  Solve  Ly = Pb for  y  using forward substitution.
	    4.  Solve  Ux = y for x  using back substitution.

	    -note y = Ux
	*/

	/*
	//create lower triangluar matrix base (create identity matrix)
	Matrix L(this->nrow,this->ncol);
	for(int i=0; i<L.nrow; i++)
		L.M[i][i] = 1;

	//create upper triangular Matrix to save A
	Matrix U(this->nrow,this->ncol);
	U.M = this->M;

	//create P matrix to keep track of row swaps (create identity matrix)
	Matrix P(this->nrow,this->ncol);
	for(int i=0; i<P.nrow; i++)
		P.M[i][i] = 1;
	
	// perform Gaussian elimination to convert A,b to an upper-triangular system
	for (int k=0; k<U.nrow -1; k++) {   // loop over diagonals

		//PIVOT
		// find the pivot row p (max of col) //starting from row k down
		int p=k;
		for (int i=k; i<U.nrow; i++){
		  if( abs(U.M[k][i]) >  abs(U.M[k][p]) )
			p=i;
		}

		//swap U and b and P rows around the pivot
		//U
		U.swap_rows(p,k);

		//b
		swap(b,p,k);

		//P
		P.swap_rows(p,k);

		//ELIMINATION
		//Perform eliminations using row k
		double c;	//scaler
		for (int i=k; i<U.nrow-1; i++){
			c = U.eliminate_0(k, i+1, k);
			b[i+1] = b[i+1] * c;				//adjust b

			//this->M[k][i+1] = c;				//store multipliers in col below pivot for L
			//L.M[k][i+1] = c;					//store multipliers in Lower triangular Matrix
			//NEG***???
			if (c == 1){
				L.M[k][i+1] =  0;
			}else{
				L.M[k][i+1] = -c;
			}

			cout << "c : " << c;

		}

	}

	//1.  Construct the matrices  L, U and P.
	cout << endl << "U, L, P:" << endl;
	U.print();
	cout << endl;
	L.print();
	cout << endl;
	P.print();
	cout << endl;

	cout << "vector b (modified b) :" << endl;
	print_vector(b);
	cout << endl;




	//check that L*U = A
	cout << endl << "check that L*U = A :" << endl;
	cout << endl << "L*U :" << endl;
	Matrix Test;
	Test = L.dot(U);
	Test.print();

	cout << endl << "A :" << endl;
	this->print();
	*/

	//create lower triangluar matrix base (create identity matrix)
	Matrix L;

	//create upper triangular Matrix to save A
	Matrix U;

	//create P matrix to keep track of row swaps (create identity matrix)
	Matrix P;

	this->LUP_factorization(L,U,P);

	cout << endl << "vector b: " << endl;
	print_vector(b);

	//2.  Compute the column vector  Pb.
	vector<double> Pb = P.dot(b);

	cout << endl << "Pb vector:" << endl;
	print_vector(Pb);
	cout << endl;

	//3.  Solve  Ly = Pb for  y  using forward substitution.	(incorrect)
	vector<double> y = L.forward_substitution(Pb);

	cout << endl << "y vector:" << endl;
	print_vector(y);
	cout << endl;

	//4.  Solve  Ux = y for x  using back substitution.
	vector<double> computed_x = U.back_substitution(y);

	cout << endl << "computed vector x:" << endl;
	print_vector(computed_x);
	cout << endl;


	/*
	//Test see if modified forward substitution gives the same answer as back sub
	computed_x = U.forward_substitution(y);

	cout << endl << "computed vector x (using forward_substitution):" << endl;
	print_vector(computed_x);
	*/






	/*
	//back substitution
	// perform Backward Substitution on result
	this->print();
	U.print();
	L.print();
	P.print();

	cout << " vector b (modified b) :" << endl;
	print_vector(b);

	cout << endl << " computed vector x :" << endl;
	vector<double> computed_x =  U.back_substitution(b);
	print_vector(computed_x);

	cout << endl << " vector b (recomputed from computed x) :" << endl;
	print_vector(U.dot(computed_x));
	*/
	/*
	// perform Gaussian elimination to convert A,b to an upper-triangular system
	for (int k=0; k<this->nrow -1; k++) {   // loop over diagonals

		//PIVOT
		// find the pivot row p (max of col) //starting from row k down
		int p=k;
		for (int i=k; i<this->nrow; i++){
		  if( abs(this->M[k][i]) >  abs(this->M[k][p]) )
			p=i;
		}

		//swap A and b rows around the pivot
		//A
		this->swap_rows(p,k);

		//b
		for (int j=0; j<b.size(); j++){
			double tmp = b[p];
			b[p] = b[k];
			b[k] = tmp;
		}

		//ELIMINATION
		//Perform eliminations using row k
		double c;	//scaler
		for (int i=k; i<this->nrow-1; i++){
			c = this->eliminate_0(k, i+1, k);
			b[i+1] = b[i+1] * c;				//adjust b
			//this->M[k][i+1] = c;				//store multipliers in col below pivot for L
			L.M[k][i+1] = c;					//store multipliers in Lower triangular Matrix
		}

	}

	//back substitution
	// perform Backward Substitution on result
	this->print();
	L.print();

	cout << " vector b (modified b) :" << endl;
	print_vector(b);

	cout << endl << " computed vector x :" << endl;
	vector<double> computed_x =  this->back_substitution(b);
	print_vector(computed_x);

	cout << endl << " vector b (recomputed from computed x) :" << endl;
	print_vector(this->dot(computed_x));
	*/


	return x;
}

//back substitution (solving Ux =b)
std::vector<double> Matrix::back_substitution(std::vector<double> b){
	
	//check for correct dimensions
	if (this->nrow != b.size() || this->nrow != this->ncol){
		cerr << "Matrix::back_substitution error, illegal matrix/vector dimensions\n";
		return vector<double>(0);
	}


  // copy b into x
  vector<double> x = b;

  // analyze matrix for magnitude
  double Umax = inifinity_norm();

  // perform column-oriented Backward Substitution algorithm
  for (int j=this->nrow-1; j>=0; j--) {

    // check for nonzero matrix diagonal
    if (abs(this->M[j][j]) < STOL*Umax) {
      cerr << "BackSubstitution error: numerically singular matrix\n";
      return vector<double>(0);
    }

    x[j] /= this->M[j][j];

    for (int i=this->ncol; i>=j; i--){
    	int iter =(this->ncol)-(j+1);
      	x[iter] -= this->M[j][i]*x[j];
    }

    /*
    // solve for this row of solution
    x[j] /= this->M[j][j];

    // update all remaining rhs
    //for (int i=0; i<j; i++)
    //x[j] -= this->M[i][j]*x[j];
    for (int i=this->ncol; i>j; i--){
    	int iter =(this->ncol)-(j+1);
      x[iter] -= this->M[j][i]*x[j];
    }
    */
  }

  return x;
}

//forward substitution (solving Lx =b)
std::vector<double> Matrix::forward_substitution(std::vector<double> b){
	
	//exit condition
	if (this->nrow != b.size()) {
		cerr << "Matrix::forward_substitution error, illegal matrix/vector dimensions\n";
		return vector<double>(0);	
	}

	// copy b into x
	vector<double> x = b;

	// analyze matrix for magnitude
	double Lmax = inifinity_norm();

	// perform column-oriented Forwards Substitution algorithm
	for (int j=0; j<this->nrow; j++) {

		// check for nonzero matrix diagonal
	    if (abs(this->M[j][j]) < STOL*Lmax) {
			cerr << "ForwardSubstitution error: singular matrix\n";
			return vector<double>(0);
		}

		// solve for this row of solution
		x[j] /= this->M[j][j];

		// update all remaining rhs
		for (int i=j+1; i<this->nrow; i++)
			x[j] -= this->M[j][i]*x[i];
	}

	return x;
}

//infinity norm (max sum of rows)
double Matrix::inifinity_norm(){
	vector<double> sums;
	this->T();
	for (int i=0; i<this->ncol; i++)
		sums.push_back(sum_vector(this->M[i]));
	this->T();

	return sums[max_index(sums)];
}

//return the index of the maximum element in the column
int Matrix::col_max_index(int col){return max_index(this->M[col]);}

void Matrix::print(){
	for (int r=0; r<nrow; r++){
		for (int c=0; c<ncol; c++){
			cout<< M[c][r] << "\t";
		}
		cout << endl;
	}
}