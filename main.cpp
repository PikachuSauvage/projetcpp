#include "Envir.h"
#include "Ecoli.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>

using std::cout;
using std::endl;
//using std::cin;
//void dynamic_param(double tmax, double tmin, double amin, double amax){
	 

int main(){
	int refeed_count = 0;
	unsigned long seed = unsigned (std::time(0));
	std::srand (seed);
	cout << "Hello projet C++!" << endl;
	std::ifstream config;
	config.open("config.txt",std::ios::in);
	std::ofstream output;
	output.open("output.txt",std::ios::out);
	int W =0;
	int H =0; 
	double D=0;
	double A_percent = 0;
	double RAA=0;
	double RAB=0;
	double RBB=0;
	double RBC=0;
	double Wmin=0;
	std::string key;
	if (config.is_open()){
		std::getline(config, key, '=');
		config>>W;
		std::getline(config, key, '=');
		config>>H;
		std::getline(config, key, '=');
		config>>D;
		std::getline(config, key, '=');
		config>>A_percent;
		std::getline(config, key, '=');
		config>>RAA;
		std::getline(config, key, '=');
		config>>RAB;
		std::getline(config, key, '=');
		config>>RBB;
		std::getline(config, key, '=');
		config>>RBC;
		std::getline(config, key, '=');
		config>>Wmin;
	}
	//cout<< H << D << A_homo << RAA << RAB << RBB << RBC << Wmin; 
	for (int T=1; T<=1500 ;T = T+20){
		for (int A=0; A<=50; A=A+4){
			cout<<T<< " "<< A <<endl;
			Envir test = Envir(W,H,D,A,RAA,RAB,RBB,RBC,Wmin);
			test.init(W,H,A_percent);
			test.updateMetab(0.1, 1);
			for (int i=0; i<10000; i++){
				//cout<< "==========="<< i << "=========="<<endl;
				test.diffuse();
				test.updateFitness();
				test.toSurvive(0);
				if (test.getStatus()==0)
					break;
				test.plsDie(0.02);
				test.updateMetab(0.1, 1);
				test.reinit();
				refeed_count++;
				if (refeed_count==T){
					refeed_count=0;
					test.refeed(A);
				}
			}
			test.result();
			output << T << " " << A << " " <<test.getStatus()<<endl;
			//test.print();
		}
	}
	config.close();
	output.close();
	//test.print();
	return 0;
}
