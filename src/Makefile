all:	prog

refresh: sequence.cpp globals.cpp helper.cpp msa.cpp main.cpp experiment.cpp
	g++ -o prog sequence.cpp globals.cpp helper.cpp msa.cpp main.cpp experiment.cpp -g
prog:	sequence.o globals.o helper.o msa.o main.o experiment.o
	g++ -o prog sequence.o globals.o helper.o msa.o main.o experiment.o -ggdb -pg -O0
sequence.o:	sequence.cpp
	g++ -c sequence.cpp -g
globals.o:	globals.cpp
	g++ -c globals.cpp -g
helper.o:	helper.cpp
	g++ -c helper.cpp -g
msa.o:	msa.cpp
	g++ -c msa.cpp -g
main.o:	main.cpp
	g++ -c main.cpp -g
experiment.o: experiment.cpp
	g++ -c experiment.cpp -g