// ===========================================================================
//                                  Includes
// ===========================================================================
#include "Envir.h"
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
// =========================================================================
//                              Protected Methods
// ========================================================================= 

// =========================================================================
//                               Constructors
// =========================================================================
Envir::Envir(unsigned int W, unsigned int H, double D, double A_MAX){
	double* Qa= new double[W*H];
	double* Qb= new double[W*H];
	double* Qc= new double[W*H];
	Ecoli** indiv= new Ecoli*[W*H];
	unsigned int* reproduced= new unsigned int[W*H];
	cout.precision(3);
	for (unsigned int i=0;i<W;i++){
		for (unsigned int j=0;j<H;j++){
			Qa[i*H + j] = (double)rand()/RAND_MAX * A_MAX;
			cout<< Qa[i*H + j] <<" ";
			Qb[i*H + j] = 0;
			Qc[i*H + j] = 0;
			reproduced[i*H + j] = 0;
			indiv[i*H + j]=nullptr;
		}
		cout << endl;
	}
	Qa_=Qa;
	Qb_=Qb;
	Qc_=Qc;
	D_=D;
	reproduced_=reproduced;
	indiv_=indiv;
	W_=W;
	H_=H;
}
// ===========================================================================
//                                 Destructor
// ===========================================================================

// =========================================================================
//                              Public Methods
// =========================================================================
Envir::init(int W, int H, double percentageGA){
	char* tmpTable= new char [W*H];
	int count=0;
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			char[i*H + j]='0';
	while (count < (W * H * percentageGA)){
		x=rand()%W;
		y=rand()%H;
		if (tmpTable[x*H + y] == '0'){
			indiv_[x*H + y]=Ecoli(x,y,'A');
			count++;
		}
	}
	for (int i=0; i<w; i++)
		for (int j=0; j<w; j++)
			if (tmpTable[i*H + j] == '0')
				indiv_[i*H + j] == Ecoli(i,j,'B');
	
}
// =========================================================================
//                                  Getters
// =========================================================================
double Envir::getQa(unsigned int x, unsigned int y){
	return Qa_[x*H_ + y];
}
double Envir::getQb(unsigned int x, unsigned int y){
	return Qb_[x*H_ + y];
}
double Envir::getQc(unsigned int x, unsigned int y){
	return Qc_[x*H_ + y];
}
unsigned int Envir::getRep(unsigned int x, unsigned int y){
	return reproduced_[x*H_ + y];
}
Ecoli* Envir::getEcoli(unsigned int x, unsigned int y){
	return indiv_[x*H_ + y];
}


// =========================================================================
//                                  Setters
// =========================================================================
void Envir::setQa(unsigned int x, unsigned int y, double newValue){
	Qa_[x*H_ + y]=newValue;
}
void Envir::setQb(unsigned int x, unsigned int y, double newValue){
	Qb_[x*H_ + y]=newValue;
}
void Envir::setQc(unsigned int x, unsigned int y, double newValue){
	Qb_[x*H_ + y]=newValue;
}
