# %load Makefile
#############################################################
######
# Makefile for project 1
# Laurence Brown
# SMU Mathematics
# Math 3316
# 9/14/2015
#############################################################
######
# compiler & flags
CXX = g++
CXXFLAGS = -O -std=c++0x
# makefile targets
#proj1_a.exe	:	proj1_a.cpp
#proj1_b.exe	:	proj1_b.cpp

all:
	$(CXX) $(CXXFLAGS) proj1_a.cpp -o proj1_a.exe
	$(CXX) $(CXXFLAGS) proj1_b.cpp -o proj1_b.exe

clean :
	rm -f *.exe *.txt
####### End of Makefile #######