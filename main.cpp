#include "Envir.h"
#include "Ecoli.h"
#include <cstdlib>
#include <iostream>
#include <ifstream>
#include <ctime>

using std::cout;
using std::endl;
//void dynamic_param(double tmax, double tmin, double amin, double amax){
	 

int main(){
	int refeed_count = 0;
	unsigned long seed = unsigned (std::time(0));
	std::srand (seed);
	cout << "Hello projet C++!" << endl;
	ifstream config ("config.txt");
	unsigned int W;
	unsigned int H; 
	double D=0;
	double A_homo = 0;
	double RAA=0;
	double RAB=0;
	double RBB=0;
	double RBC=0;
	double Wmin=0;
	if (config.if_open()){
		cout<<"W="<<W<<"\n";
		cout<<"H="<<H<<"\n";
		cout<<"D="<<D<<"\n";
		cout<<"A_homo="<<A_homo<<"\n";
		cout<<"RAA="<<RAA<<"\n";
		cout<<"RAB="<<RAB<<"\n";
		cout<<"RBB="<<H<<"\n";
		cout<<"H="<<H<<"\n";
		
	Envir test = Envir(32,32,0.1,50,0.1,0.1,0.1,0.1,0.001);
	//cout << test.getQa(0,0) << endl;
	test.init(32,32,0.5);
	test.updateMetab(0.1, 1);
	//test.print();
	for (int i=0; i<10000; i++){
		cout<< "==========="<< i << "=========="<<endl;
		test.diffuse();
		//test.print();
		test.updateFitness();
		test.toSurvive(0);
		test.plsDie(0.02);
		//cout << " After plsDie, after toSurvive " <<endl;
		//test.print();
		test.updateMetab(0.1, 1);
		test.reinit();
		refeed_count++;
		if (refeed_count==1300){
			refeed_count=0;
			test.refeed(50);
		}
	}
	test.print();
	return 0;
}
