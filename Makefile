all: main

main:  main.o Envir.o Ecoli.o
	g++ -g main.o Envir.o Ecoli.o -o main
	
main.o: main.cpp Ecoli.h Envir.h
	g++ -g -c -Wall main.cpp -o main.o -std=c++11

Envir.o: Envir.cpp Envir.h Ecoli.h
	g++ -g -c -Wall Envir.cpp -o Envir.o -std=c++11
	
Ecoli.o: Ecoli.cpp Ecoli.h
	g++ -g -c -Wall Ecoli.cpp -o Ecoli.o -std=c++11
	
clean:
	rm -f *.o
