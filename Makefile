all: main

main:  main.o Run.o Envir.o Ecoli.o
	g++ main.o Run.o Envir.o Ecoli.o -o main -O3
	
main.o: main.cpp Envir.h Ecoli.h
	g++ -c main.cpp -o main.o -std=c++11 -O3
	
Run.o: Run.cpp Run.h Envir.h
	g++ -c Run.cpp -o Run.o -std=c++11 -O3
	
Envir.o: Envir.cpp Envir.h Ecoli.h
	g++ -c Envir.cpp -o Envir.o -std=c++11 -O3
	
Ecoli.o: Ecoli.cpp Ecoli.h
	g++ -c Ecoli.cpp -o Ecoli.o -std=c++11 -O3
	
clean:
	rm -f *.o
