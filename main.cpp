#include "Envir.h"
#include "Ecoli.h"
#include <cstdlib>
#include <iostream>
#include <ctime>

using std::cout;
using std::endl;
int main(){
	unsigned long seed = unsigned ( std::time(0) );
	std::srand (seed);
	cout << "Hello projet C++!" << endl;
	Envir test = Envir(3,3,0.1,50,0.1,0.1,0.1,0.1,0.01);
	//cout << test.getQa(0,0) << endl;
	test.init(3,3,0.5);
	test.updateMetab();
	test.print();
	for (int i=0; i<2; i++){
		cout<< "==========="<< i << "=========="<<endl;
		test.diffuse();
		test.print();
		test.updateFitness();
		test.toSurvive(0.1);
		test.plsDie(0.2);
		test.print();
	}
	return 0;
}
