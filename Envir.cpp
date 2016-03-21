// ===========================================================================
//                                  Includes
// ===========================================================================
#include "Envir.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
using std::cout;
using std::endl;
// =========================================================================
//                              Protected Methods
// ========================================================================= 

// =========================================================================
//                               Constructors
// =========================================================================
Envir::Envir(unsigned int W, unsigned int H, double D, double A_homo,
			double RAA, double RAB, double RBB, double RBC, double Wmin){
	Qa_= new double[W*H];
	Qb_= new double[W*H];
	Qc_= new double[W*H];
	F_= new double[W*H];
	Ecoli** indiv= new Ecoli*[W*H];
	unsigned int* reproduced= new unsigned int[W*H];
	cout.precision(3);
	for (unsigned int i=0;i<W;i++){
		for (unsigned int j=0;j<H;j++){
			Qa_[i*H + j] = A_homo;
			cout<< Qa_[i*H + j] <<" ";
			Qb_[i*H + j] = 0;
			Qc_[i*H + j] = 0;
			F_[i*H + j] =0;
			reproduced[i*H + j] = 0;
			indiv[i*H + j]=nullptr;
		}
		cout << endl;
	}
	D_=D;
	reproduced_=reproduced;
	indiv_=indiv;
	W_=W;
	H_=H;
	RAA_=RAA;
	RAB_=RAB;
	RBB_=RBB;
	RBC_=RBC;
	Wmin_=Wmin;
}
// ======================================================================
//                                 Destructor
// =======================================================================

// =========================================================================
//                              Public Methods
// =========================================================================
void Envir::init(int W, int H, double percentageGA){
	char* tmpTable= new char [W*H];
	int x=0,y=0;
	int count=0;
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			tmpTable[i*H + j]='0';
	while (count < (W * H * percentageGA)){
		x=rand()%W;
		y=rand()%H;
		if (tmpTable[x*H + y] == '0'){
			Ecoli thisecoli=Ecoli(x,y,'A');
			indiv_[x*H + y]= &thisecoli;
			count++;
		}
	}
	for (int i=0; i<W; i++)
		for (int j=0; j<W; j++)
			if (tmpTable[i*H + j] == '0'){
				Ecoli thisecoli=Ecoli(x,y,'B');
				indiv_[i*H + j] = &thisecoli;
			}
}

void Envir::updateMetab(){
	double curr_qa = 0;
	double curr_qb = 0;
	double curr_qc = 0;
	for (unsigned int x=0; x<W_; x++)
		for (unsigned int y=0; y<H_; y++){
			curr_qa = indiv_[x*H_ + y]->getQa();
			curr_qb = indiv_[x*H_ + y]->getQb();
			curr_qc = indiv_[x*H_ + y]->getQc();
			if (indiv_[x*H_ + y]->getGeno()=='A'){
				Qa_[x*H_ +y] = Qa_[x*H_ + y]*(1-RAA_);//dAout
				indiv_[x*H_ + y]->setQa(curr_qa+RAA_*Qa_[x*H_ +y]
				-curr_qa*RAB_); //dA
				indiv_[x*H_ + y]->setQb(curr_qb*(1+RAB_));//dB
			}
			if (indiv_[x*H_ + y]->getGeno()=='B'){
				Qb_[x*H_ +y] = Qb_[x*H_ + y]*(1-RBB_);//dBout
				indiv_[x*H_ + y]->setQb(curr_qb+RBB_*Qb_[x*H_ + y]
				-curr_qb*RBC_);//dB
				indiv_[x*H_ + y]->setQc(curr_qc*(1+RBC_));//dC
			}
		}
}

void Envir::updateFitness(){
	for (unsigned int x=0; x<W_; x++)
		for (unsigned int y=0; y<H_; y++){
			if (indiv_[x*H_ + y]->getGeno()=='A'){
				F_[x*H_ + y]=indiv_[x*H_ + y]->getQb();
			}
			if (indiv_[x*H_ + y]->getGeno()=='B'){
				F_[x*H_ + y]=indiv_[x*H_ + y]->getQc();
			}
			if (F_[x*H_ + y]<Wmin_)
				F_[x*H_ + y]=0;
		}
}

void Envir::diffuse(){
	double* newQa= new double[W_*H_];
	double* newQb= new double[W_*H_];
	double* newQc= new double[W_*H_];
	for (int x=0; x<(int)W_; x++)
		for (int y=0; y<(int)H_; y++){
			memcpy(newQa, Qa_, W_*H_);
			memcpy(newQb, Qb_, W_*H_);
			memcpy(newQc, Qc_, W_*H_);
			int m=0, n=0;
			for (int i=-1; i<=1; i++)
				for (int j=-1; j<=1; j++){
					if ((x+i)<0)
						m=i+H_;
					if ((x+i)>= (int)H_)
						m=i-H_;
					if ((y+j)<0)
						n=j+W_;
					if ((y+j)>=(int)H_)
						n=j-W_;
					newQa[x*H_+y] += D_*Qa_[(x+m)*H_ +(y+n)];
					newQb[x*H_+y] += D_*Qb_[(x+m)*H_ +(y+n)];
					newQc[x*H_+y] += D_*Qc_[(x+m)*H_ +(y+n)];
			}
			newQa[x*H_+y]-=9*D_*Qa_[x*H_ + y];
			newQb[x*H_+y]-=9*D_*Qb_[x*H_ + y];
			newQc[x*H_+y]-=9*D_*Qc_[x*H_ + y];
		}
}

void Envir::plsDie(double prob){
	for (int x=0; x<(int)W_; x++)
		for (int y=0; y<(int)H_; y++){
			if ((double) rand()/RAND_MAX < prob)
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
