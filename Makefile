all:rasterizer

rasterizer:
	g++ *.cpp -O3 -o rasterizer -std=c++11

clean:
	rm rasterizer