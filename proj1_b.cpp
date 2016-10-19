#include <vector>
#include <stdlib.h> /* abs */
#include <cmath> /* pow */
#include <algorithm>  /* reverse */
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

//variables:
	//smallest machine number
const double error_dp = pow(2,-52);

//computer estimates for functions f(x) = x^-3, f'(x), f''(x). 
	//*used as actual values for error calculations
double f(double x)
	{	return pow(x,-3);			}

double df(double x)
	{	return -3 * pow(x,-4);		}

double ddf(double x)
	{	return 12 * pow(x,-5);		}

//computer estimates for functions f(x) = x^^-3, f'(x), f''(x). using
	//using the forward diffrence approximation method.
double df_fda(double x, double h)
	{	return (f(x + h) - f(x)) / h;	}

//computer estimates for relitive error, r ,  of df_fda(x)
double r(double x, double  h)
	{	return abs( (df(x) - df_fda(x, h)) / f(x) );	}

//c1 and c2 for next formula

double c1(double x)
	{	return abs( ddf(x) / (2.0 * df(x)) );	}

double c2(double x, double h)
{	return abs( ( f(x) * error_dp ) / df(x) );		}
	//{	return abs( df_fda(x, h) / df(x) );		}

//computer esimates for upperbound of r, R, of df_fda(x)
double R(double x, double h)
	{	return abs( (c1(x) * h) + ( c2(x,h) * (1.0/h) ) );	}

//h = 2^-n , n = {1,2,3..52}
double h(double n)
	{	return pow(2,-n);	}

//returns a vector of n values
std::vector<double> n_v(size_t max){
	std::vector<double> answers;
	for (size_t i = 1; i <= max; i++)
		answers.push_back(i);
	return answers;
}

//returns a vector of h values
std::vector<double> h_v(std::vector<double> n_v){
	std::vector<double> answers;
	for (int i=0; i<n_v.size(); i++)
		answers.push_back(h(n_v[i]));
	return answers;
}

//returns a vector of r values
std::vector<double> r_v(double x, std::vector<double> h_v){
	std::vector<double> answers;
	for (int i=0; i<h_v.size(); i++)
		answers.push_back(r(x,h_v[i]));
	return answers;
}

//returns a vector of R values
std::vector<double> R_v(double x, std::vector<double> h_v){
	std::vector<double> answers;
	for (int i=0; i<h_v.size(); i++)
		answers.push_back(R(x,h_v[i]));
	return answers;
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

	//compute n,h,r,R vectors with a =3, 
		//and write them to a file with corresponding name.txt extenstion:
	const double a = 3;

	//n
	std::vector<double> n = n_v(52);
	write(n, "n.txt");

	//h
	std::vector<double> h = h_v(n);
	write(h, "h.txt");

	//r
	std::vector<double> r = r_v(a, h);
	write(r, "r.txt");

	//R
	std::vector<double> R = R_v(a, h);
	write(R, "R_ub.txt");

	//end program
	return 1;
}

