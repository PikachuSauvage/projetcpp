all: main

main:  main.o Envir.o Ecoli.o
	g++ main.o Envir.o Ecoli.o -o main -O3
	
main.o: main.cpp Envir.h Ecoli.h
	g++ -g -c main.cpp -o main.o -std=c++11 -O3

Envir.o: Envir.cpp Envir.h Ecoli.h
	g++ -g -c Envir.cpp -o Envir.o -std=c++11 -O3
	
Ecoli.o: Ecoli.cpp Ecoli.h
	g++  -g -c Ecoli.cpp -o Ecoli.o -std=c++11 -O3
	
clean:
	rm -f *.o
