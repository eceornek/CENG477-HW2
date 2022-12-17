all:
	g++ *.cpp -O3 -o Main -std=c++11
debug:
	g++ *.cpp -g -O0 -o Main -std=c++11