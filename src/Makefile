all:	prog

refresh: sequence.cpp globals.cpp helper.cpp msa.cpp main.cpp experiment.cpp
	g++ -o prog sequence.cpp globals.cpp helper.cpp msa.cpp main.cpp experiment.cpp -g -std=c++0x
prog:	sequence.o globals.o helper.o msa.o main.o experiment.o 
	g++ -o prog sequence.o globals.o helper.o msa.o main.o experiment.o -ggdb -pg -O0 -std=c++0x
sequence.o:	sequence.cpp 
	g++ -c sequence.cpp -g -std=c++0x
globals.o:	globals.cpp 
	g++ -c globals.cpp -g -std=c++0x
helper.o:	helper.cpp 
	g++ -c helper.cpp -g -std=c++0x
msa.o:	msa.cpp
	g++ -c msa.cpp -g -std=c++0x
main.o:	main.cpp
	g++ -c main.cpp -g -std=c++0x
experiment.o: experiment.cpp
	g++ -c experiment.cpp -g -std=c++0x