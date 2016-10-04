#include <functional> /* function type */

double newton(std::function<double(double)> func, std::function<double(double)> dfunc, double x, int maxit, double tol, bool show_iterates);