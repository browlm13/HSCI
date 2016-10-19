/*

	polynomial interpollation with newton form

*/

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

void print_vector(vector<double> v){
	for (int i=0; i<v.size(); i++)
		cout << v[i] << " ";
	cout << endl;
}

double Newton_basis(vector<double> xnodes, int n, double x){
	//phi_n+1(x) = (x - x0)(x - x1) ... (x - xn)
	double basis = 1;
	for (int i=0; i<n; i++)
		basis *= (x - xnodes[i]);
	return basis;
}

double next(double prev, int n, double prev_a, vector<double> xnodes, double x){
	//pn(x) = pn-1(x) + an(x-x0)(x-x1)///(x-xn-1)
	double tail = prev_a * Newton_basis(xnodes, n, x);
	return prev + tail;
}


double f(double x){ return ((3.1) * pow(x,4)) + ((2.3) * pow(x,3)) - ((6.6) * pow(x,2)) + (8.7 *x) + 7.9;}

vector<double> gen_ynodes(vector<double> xnodes){
	//gernerate y nodes: f(x) = 3.1x^4 + 2.3x^3 − 6.6x^2 + 8.7x + 7.9

	vector<double> ynodes;
	for (int i=0; i<xnodes.size(); i++)
		ynodes.push_back(f(xnodes[i]));
	return ynodes;
}

vector<double> Newton_co(vector<double> xnodes, vector<double> ynodes, double x){
	//returns coefficeients
	int n = xnodes.size();

	vector<double> ps;
	vector<double> as;

	ps.push_back(ynodes[0]);
	as.push_back(ynodes[0]);

	for(int i=0; i<=n; i++){
		as.push_back( (ynodes[i+1] - ps[i])  / Newton_basis(xnodes, i, x));
		ps.push_back( next(ps[i], i, as[i], xnodes, x) );
	}

	return as;

}

int main(){

	vector<double> xnodes = {-2, -1, 0, 1 ,2};
	vector<double> ynodes = gen_ynodes(xnodes);

	print_vector(xnodes);
	print_vector(ynodes);

	double x = 201;
	vector<double> as = Newton_co(xnodes, ynodes, x);
	print_vector(as);

	cout << "tiger";
	return 0;
}