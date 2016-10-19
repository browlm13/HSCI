# %load Makefile
#############################################################
######
# Makefile for project 2
# Laurence Brown
# SMU Mathematics
# Math 3316
# 9/28/2015
#############################################################
######
# compiler & flags
CXX = g++
CXXFLAGS = -O -std=c++0x
# makefile targets
#newton_test.exe	:	newton.cpp	test_newton.cpp
#linear_solve.exe	:	matrix.cpp	vandermonde.cpp		
#kepler.exe	:	matrix.cpp	newton.cpp	kepler.cpp

all:
	$(CXX) $(CXXFLAGS) newton.cpp test_newton.cpp -o test_newton.exe
	$(CXX) $(CXXFLAGS) matrix.cpp	vandermonde.cpp -o linear_solve.exe
	$(CXX) $(CXXFLAGS) matrix.cpp	newton.cpp kepler.cpp -o kepler.exe
	
	chmod 755 test_newton.exe
	chmod 755 linear_solve.exe
	chmod 755 kepler.exe

clean :
	rm -f *.exe *.txt
####### End of Makefile #######
