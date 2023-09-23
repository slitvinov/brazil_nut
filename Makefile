all: main
.cpp:
	c++ -O2 $< -lglut -lGL -o $@
