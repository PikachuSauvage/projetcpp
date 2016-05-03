// ===========================================================================
//                                  Includes
// ===========================================================================
#include "Envir.h"
#include <algorithm>
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
	reproduced_= new unsigned int[W*H];
	cout.precision(3);
	for (unsigned int i=0;i<W;i++){
		for (unsigned int j=0;j<H;j++){
			Qa_[i*H + j] = A_homo;
			//cout<< Qa_[i*H + j] <<" ";
			Qb_[i*H + j] = 0;
			Qc_[i*H + j] = 0;
			F_[i*H + j] =0;
			reproduced_[i*H + j] = 0;
			indiv[i*H + j]=nullptr;
		}
		//cout << endl;
	}
	D_=D;
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
Envir::~Envir(){
	delete[] Qa_;
	delete[] Qb_;
	delete[] Qc_;
	delete[] F_;
	delete[] reproduced_;
	for (unsigned int i=0; i<W_; i++)
		for (unsigned int j=0; j<H_; j++)
			delete indiv_[i*H_+j];
	delete[] indiv_;
	Qa_=nullptr;
	Qb_=nullptr;
	Qc_=nullptr;
	F_=nullptr;
	reproduced_=nullptr;
	indiv_=nullptr;
}
// =========================================================================
//                              Public Methods
// =========================================================================
void Envir::init(int W, int H, double percentageGA){
	char* tmpTable= new char[W*H];
	int x=0,y=0;
	int count=0;
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			tmpTable[i*H + j]='0';
	while (count < (W * H * percentageGA)){
		x=rand()%W;
		y=rand()%H;
		//cout << x <<','<< y << endl;
		if (tmpTable[x*H + y] == '0'){
			indiv_[x*H + y]= new Ecoli(x,y,'A');
			tmpTable[x*H + y] ='A';
			count++;
		}
	}
	for (int i=0; i<W; i++)
		for (int j=0; j<W; j++)
			if (tmpTable[i*H + j] == '0')
				indiv_[i*H + j] = new Ecoli(x,y,'B');
	delete[] tmpTable;
}

void Envir::updateMetab(double p, double n){
	double curr_qa = 0;
	double curr_qb = 0;
	double curr_qc = 0;
	double Qmoving = 0;
	for (unsigned int x=0; x<W_; x++)
		for (unsigned int y=0; y<H_; y++){
			if (reproduced_[x*H_ + y]==0){
				if (indiv_[x*H_ + y]->getGeno()=='A'){
					for (int i=0; i< (int)(n/p); i++){
						indiv_[x*H_ + y]->setQb(indiv_[x*H_ + y]->getQb()
						+indiv_[x*H_ + y]->getQa()* RAB_ * p);
						indiv_[x*H_ + y]->setQa(indiv_[x*H_ + y]->getQa()
						+(Qa_[x*H_ + y]*RAA_-indiv_[x*H_ + y]->getQa()*RAB_)*p);
						Qa_[x*H_+y]-=Qa_[x*H_+y]*RAA_*p;
					}
				}
				if (indiv_[x*H_ + y]->getGeno()=='B'){
					for (int i=0; i< (int)(n/p); i++){
						indiv_[x*H_ + y]->setQc(indiv_[x*H_ + y]->getQc()
						+indiv_[x*H_ + y]->getQb()*RBC_*p);
						indiv_[x*H_ + y]->setQb(indiv_[x*H_ + y]->getQb()
						+(Qb_[x*H_ + y]*RBB_-indiv_[x*H_ + y]->getQb()*RBC_)*p);
						Qb_[x*H_+y]-=Qb_[x*H_+y]*RBB_*p;
					}
				}
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
			if ((F_[x*H_ + y]<Wmin_)||(reproduced_[x*H_+y]==1))
				F_[x*H_ + y]=0;
		}
}

void Envir::diffuse(){
	double* newQa= new double[W_*H_];
	double* newQb= new double[W_*H_];
	double* newQc= new double[W_*H_];
	int m=0, n=0;
	memcpy(newQa, Qa_, W_*H_*sizeof(double));
	memcpy(newQb, Qb_, W_*H_*sizeof(double));
	memcpy(newQc, Qc_, W_*H_*sizeof(double));
	for (int x=0; x< (int)W_; x++)
		for (int y=0; y<(int)H_; y++){
			for (int i=-1; i<=1; i++)
				for (int j=-1; j<=1; j++){
					m=0;
					n=0;
					if ((x+i)<0)
						m=i+H_;
					if ((x+i)>= (int)H_)
						m=i-H_;
					if ((y+j)<0)
						n=j+W_;
					if ((y+j)>= (int)H_)
						n=j-W_;
					//cout<< x*H_+y <<','<< m << ',' << n<< endl;
					newQa[x*H_+y] += D_*Qa_[(x+m)*H_ +(y+n)];
					newQb[x*H_+y] += D_*Qb_[(x+m)*H_ +(y+n)];
					newQc[x*H_+y] += D_*Qc_[(x+m)*H_ +(y+n)];
			}
			newQa[x*H_+y]-=9*D_*Qa_[x*H_ + y];
			newQb[x*H_+y]-=9*D_*Qb_[x*H_ + y];
			newQc[x*H_+y]-=9*D_*Qc_[x*H_ + y];
		}
	delete[] Qa_;
	delete[] Qb_;
	delete[] Qc_;
	Qa_=newQa;
	Qb_=newQb;
	Qc_=newQc;
}

void Envir::plsDie(double prob){
	for (int x=0; x<(int)W_; x++)
		for (int y=0; y<(int)H_; y++)
			if ((double) rand()/RAND_MAX < prob){
				//cout<< "hi"<<endl;
				Qa_[x*H_+y] += indiv_[x*H_+y]->getQa();
				Qb_[x*H_+y] += indiv_[x*H_+y]->getQb();
				Qc_[x*H_+y] += indiv_[x*H_+y]->getQc();
				indiv_[x*H_+y]->setGeno('0');
				indiv_[x*H_+y]->setQa(0);
				indiv_[x*H_+y]->setQb(0);
				indiv_[x*H_+y]->setQc(0);
				F_[x*H_+y]=0;
			}
}
void Envir::mutation(double prob, int posi){
	if ((double) random()/RAND_MAX <= prob){
		//cout << "Congrats, your PIKACHU evolved into RAICHU" << endl;
		if (indiv_[posi]->getGeno()=='A'){
			indiv_[posi]->setGeno('B');
		} else {
			indiv_[posi]->setGeno('A');
		}
	}
}
		
void Envir::toSurvive(double prob){
	std::vector<int> v;
	int x=0;
	int y=0;
	int m=0;
	int n=0;
	double fitmax=Wmin_;
	int memmax;
	for (unsigned int i=0;i<W_;i++)
		for (unsigned int j=0;j<H_;j++){
			if (indiv_[i*H_+j]->getGeno() == '0'){
				v.push_back(i*H_+j);
			}
		}
	std::random_shuffle(v.begin(),v.end());
	//cout<< "Vector contrains:";
	//for (std::vector<int>::iterator i=v.begin();i!=v.end();++i)
		//cout<<' '<< *i;
	//cout<< '\n';
	for (std::vector<int>::iterator i=v.begin();i!=v.end();++i){
		x=*i/H_;
		y=*i-H_*x;
		//cout << "*i="<<*i<< ",X=:" << x <<",Y=:"<<y<<endl;
		for (int i=-1; i<=1; i++)
			for (int j=-1; j<=1; j++){
				m=x+i;
				n=y+j;
				if ((x+i)<0)
					m+=H_;
				if ((x+i)>= (int)H_)
					m-=H_;
				if ((y+j)<0)
					n+=W_;
				if ((y+j)>=(int)H_)
					n-=W_;
				//cout << m*H_+n <<"hi" <<endl;
				if ((reproduced_[m*H_+n]==0)&&(F_[m*H_+n] > fitmax)){
					//cout << "Mother No:"<< m*H_+n <<endl;
					memmax=m*H_+n;
					fitmax=F_[m*H_+n];
				}
			}
			//mutation
		//cout << "yo yo" << endl;
		if (fitmax>Wmin_){
			//cout << "check now" <<endl;
			mutation(prob, x*H_+y);
			mutation(prob, memmax);
			//Division
			indiv_[memmax]->setQa(indiv_[memmax]->getQa()/2);
			indiv_[memmax]->setQb(indiv_[memmax]->getQb()/2);
			indiv_[memmax]->setQc(indiv_[memmax]->getQc()/2);
			indiv_[x*H_+y]->setQa(indiv_[memmax]->getQa());
			indiv_[x*H_+y]->setQb(indiv_[memmax]->getQb());	
			indiv_[x*H_+y]->setQc(indiv_[memmax]->getQc());	
			indiv_[x*H_+y]->setGeno(indiv_[memmax]->getGeno());
			reproduced_[memmax]=1;
			reproduced_[x*H_+y]=1;
		}
	}
}

void Envir::print(){
	int A=0;
	int B=0;
	double somme=0;
	for (unsigned int i=0;i<W_;i++){
		for (unsigned int j=0;j<H_;j++){
			cout<<"No:" << i*H_+j <<
			" Geno: "<<indiv_[i*H_+j]->getGeno() <<
			//" Fit: "<<F_[i*H_+j]<<
			" QaEn: " <<Qa_[i*H_ + j] <<
			" QbEn: " <<Qb_[i*H_ + j] <<
			" QcEn: " <<Qc_[i*H_ + j] <<
			" QaIn: " <<indiv_[i*H_+j]->getQa()<<
			" QbIn: " <<indiv_[i*H_+j]->getQb()<<
			" QcIn: " <<indiv_[i*H_+j]->getQc()<<
			" rp " <<reproduced_[i*H_ + j]<<endl;
			somme=somme+Qa_[i*H_ + j]+Qb_[i*H_ + j]+Qc_[i*H_ + j];
			somme=somme+indiv_[i*H_+j]->getQa()+indiv_[i*H_+j]->getQb()+indiv_[i*H_+j]->getQc();
			if (indiv_[i*H_ + j]->getGeno()=='A')
				A++;
			else 
				if (indiv_[i*H_ + j]->getGeno()=='B')
					B++;
		}
	}
	cout<<"A: "<< A << " B: " << B << " somme:" << somme <<endl;
}

void Envir::reinit(){
	for (unsigned int i=0;i<W_;i++)
		for (unsigned int j=0;j<H_;j++)
			reproduced_[i*H_ + j] = 0;
		
}

void Envir::refeed(double A){
	for (unsigned int i=0;i<W_;i++)
		for (unsigned int j=0;j<H_;j++){
			Qa_[i*H_ + j] = A;
			Qb_[i*H_ + j] = 0;
			Qc_[i*H_ + j] = 0;
		}
}
void Envir::run(int TMAX){
	for (int t=1; t<=TMAX; t++){
		updateMetab(0.1, 1);
		plsDie(0.1);
		toSurvive(0);
	}
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
//unsigned int Envir::getRep(unsigned int x, unsigned int y){
	//return reproduced_[x*H_ + y];
//}
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
