#include "Envir.h"
#include "Ecoli.h"
#include <cstdlib>
#include <iostream>

using std::cout;
using std::endl;
int main(){
	cout << "Hello projet C++!" << endl;
	Envir test = Envir(32,32,0.1,50);
	cout << test.getQa(0,0) << endl;
	return 0;
}
