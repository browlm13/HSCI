/*

	polynomial interpollation with newton form

*/

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;


/*
	additional functions
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

//colapse vector to double, (sum entries)
double sum_vector(vector<double> v){
	double  sum = 0;
	for (int i=0; i<v.size(); i++)
		sum += v[i];
	return sum;
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

//log vector
void print_vector(vector<double> v){
	for (int i=0; i<v.size(); i++)
		cout << v[i] << " ";
	cout << endl;
}


/*
	functions to test
*/

//actual function value
double f(double x){ return ((3.1) * pow(x,4)) + ((2.3) * pow(x,3)) - ((6.6) * pow(x,2)) + (8.7 *x) + 7.9;}

double Newton_basis(vector<double> xnodes, int n, double x){
	//phi_n+1(x) = (x - x0)(x - x1) ... (x - xn)
	double basis = 1;
	for (int i=0; i<=n; i++)
		basis *= (x - xnodes[i]);
	return basis;
}

double next(double prev, int n, double prev_a, vector<double> xnodes, double x){
	//pn(x) = pn-1(x) + an(x-x0)(x-x1)///(x-xn-1)
	double tail = prev_a * Newton_basis(xnodes, n, x);
	return prev + tail;
}

vector<double> gen_ynodes(vector<double> xnodes){
	//gernerate y nodes: f(x) = 3.1x^4 + 2.3x^3 − 6.6x^2 + 8.7x + 7.9

	vector<double> ynodes;
	for (int i=0; i<xnodes.size(); i++)
		ynodes.push_back(f(xnodes[i]));
	return ynodes;
}

//returns coefficeients for newton form
vector<double> Newton_co(vector<double> xnodes, vector<double> ynodes, double x){
	int n = xnodes.size();

	vector<double> ps;
	vector<double> as;		//coefficeients

	ps.push_back(ynodes[0]);
	as.push_back(ynodes[0]);

	for(int i=0; i<=n; i++){
		//as.push_back( (ynodes[i+1] - ps[i])  / Newton_basis(xnodes, i, x));
		as.push_back( (ynodes[i+1] - ps[i+1])  / Newton_basis(xnodes, i+1, x));

		ps.push_back( next(ps[i], i, as[i], xnodes, x) );
	}

	return as;

}

//used to obtain interpolated value after coefiencts are obtained
double Netwon_polynomial(vector<double> xnodes, vector<double> as, double x){
	//evalutes pn(x) = a0 + a1(x − x0) + a2(x − x0)(x − x1) + · · ·+an(x − x0)(x − x1)· · ·(x − xn−1).
	double pn;
	for (int i=0; i<as.size(); i++)
		pn += as[i] * Newton_basis(xnodes, i, x);
	return pn;
}

int main(){

	//data set
	vector<double> xnodes = {-2, -1, 0, 1 ,2};		//given xs
	vector<double> ynodes = gen_ynodes(xnodes);		//given ys

	//inputs to function
	vector<double> xs = Linspace(-3,3,201);

	vector<double> ps;
	vector<double> fxs;
	for (int i=0; i<xs.size(); i++){

		//compute actual value
		fxs.push_back(f(xs[i]));

		/*
			Newton interpolation
		*/

		//get coeffients
		vector<double> as = Newton_co(xnodes, ynodes, xs[i]);

		//interpolate
		ps.push_back(Netwon_polynomial(xnodes, as, xs[i]));
	}

	//write inputs to file
	write(xs, "xs.txt");

	//write actual function values
	write(fxs, "fxs.txt");

	//write Newton Ploynomial function values
	write(ps, "ps.txt");

	return 0;
}