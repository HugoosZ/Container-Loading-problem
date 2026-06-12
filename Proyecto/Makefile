all: bsg_rcl
bsg_rcl: bsg_rcl.cpp
	g++ -O2 -std=c++17 -o bsg_rcl bsg_rcl.cpp
test: bsg_rcl
	./bsg_rcl --selftest
clean:
	rm -f bsg_rcl
