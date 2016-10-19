#include <vector>
#include <stdlib.h> /* abs */
#include <cmath> /* pow */
#include <iostream>
#include <string>
#include <fstream>

using namespace std;


//create a new vector of linearly space data
std::vector<double> Linspace(double a, double b, size_t n){
	if (n<2) std::cerr << "Linspace::length must be > 1\n";
	std::vector<double> v(n);
	
	double h = (b-a)/(n-1);
	for (size_t i=0; i<n; i++)
	 	v[i] = a + i*h;

	return v;
}

//vector element wise subtraction function
vector<double> sub(vector<double> a, vector<double> b){
	if (a.size() != b.size()) std::cerr << "sub::vectors must be same sime";

	vector<double> x;
	
	for(size_t i=0; i<a.size(); i++)
		x.push_back(a[i]-b[i]);
	
	return x;

}

//factorial function
double factorial(size_t n){
	//n = abs(n);
	if (n < 1) return 1;
	double x=n;
	for(size_t i=(n-1); i>1; i--)
		x *= i;
	return x;
}

//nest function without use of matrix class
double nest(std::vector<double> v, double x){
	if (v.size() == 0) std::cerr << "nest::vector length must be > 0";
	double p = v[0];
	for(size_t i=1; i<v.size(); i++)
	 	p += (v[i] * pow(x,i));
	return p;
}

std::vector<double> pn(size_t n, std::vector<double> z){

	//vector of answers from every trial
	vector<double> v(z);

	//generate a values
	std::vector<double> a_v(n);
	
	for(size_t i=0; i<n; i++)
		a_v[i] = (1/(factorial(i)));
	
	//fill andwers vector by running nest function on each z value
	for( size_t i=0; i<z.size(); i++)
		v[i]=(nest(a_v, z[i]));
	
	return v;
}

//c++'s exponential funciton
std::vector<double> e(std::vector<double> z){

	//vector of answers from every trial
	vector<double> v(z);

	//fill anwers vector by running exp function on each z value
	for( size_t i=0; i<z.size(); i++)
		v[i]=(double)exp(z[i]);
	
	return v;
}

//find the error between two vectors
std::vector<double> err(std::vector<double> a, std::vector<double> b){
	if (a.size() != b.size()) std::cerr << "err::vectors must be same sime";
	std::vector<double> x = sub(a,b);

	for(size_t i=0; i<a.size(); i++)
		x[i] = abs(x[i]);
	
	return x;

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


//main function
int main(){
	

	//create z
	size_t size = 201;
	std::vector<double> z = Linspace(-1.0, 1.0, size);

	//p4,p8,p12,e^x
	vector<double> p4 = pn(4,z);
	vector<double> p8 = pn(8,z);
	vector<double> p12 = pn(12,z);
	vector<double> f = e(z);

	//err4,err8,err12
	vector<double> err4 = err(f,p4);
	vector<double> err8 = err(f,p8);
	vector<double> err12 = err(f,p12);
	

	//write to files corresponding names
	write(z, "z.txt");
	write(p4, "p4.txt");
	write(p8, "p8.txt");
	write(p12, "p12.txt");
	write(f, "f.txt");
	write(err4, "err4.txt");
	write(err8, "err8.txt");
	write(err12, "err12.txt");

	return 1;
}
