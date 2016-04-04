#include "Envir.h"
#include "Ecoli.h"
#include <cstdlib>
#include <iostream>

using std::cout;
using std::endl;
int main(){
	cout << "Hello projet C++!" << endl;
	Envir test = Envir(10,10,0.1,50,0.1,0.1,0.1,0.1,0.01);
	cout << test.getQa(0,0) << endl;
	test.init(10,10,0.5);
	test.updateMetab();
	test.print();
	test.plsDie(0.1);
	test.print();
	return 0;
}
