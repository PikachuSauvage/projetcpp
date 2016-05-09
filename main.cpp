#include "Run.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
//#include <boost/tuple/tuple.hpp>
//#include "gnuplot-iostream/gnuplot-iostream.h"
using std::cout;
using std::endl;

//input: int Astart, Aend, double D, int Tstart, int Tend,double Apas 
//       int Tpasint rep, bool generate_squares, bool use_squares
int main(int argc, char *argv[]){
	std::string s_output;
	std::string s_squars;
	std::ofstream output;
	std::vector<std::tuple<int, int, int, int,double>> sqs;
	double A1=0;
	double A2=50;
	double Apas=5;
	int rep=5;
	int T1=1;
	int T2=1500;
	int Tpas=50;
	double D=atof(argv[3]);
	double tmpres=0;
	int generate_squares=1;
	int use_squares=1;
	unsigned long seed = unsigned (std::time(0));
	std::srand (seed);
	cout << "Hello projet C++!" << endl;
	std::string s("config.txt");
	s_output= "manualoutput"+ std::string(argv[1]) +
			"_"+std::string(argv[2])+"_"+std::string(argv[3])+".txt";
	output.open(s_output.c_str());
	if (!output){
		cout<<"File not opened!"<<endl;
	}
	if (argc>=2) A1=atof(argv[1]);
	if (argc>=3) A2=atof(argv[2]);
	if (argc>=4) D=atof(argv[3]);
	if (argc>=5) T1=atoi(argv[4]);
	if (argc>=6) T2=atoi(argv[5]);
	if (argc>=7) Apas=atof(argv[6]);
	if (argc>=8) Tpas=atoi(argv[7]);
	if (argc>=9) rep=atoi(argv[8]);
	if (argc>=10) generate_squares=atoi(argv[9]);
	if (argc>=11) use_squares=atoi(argv[10]);
	
	cout<<"Astart="<< A1 << " Aend=" << A2 << " D=" << D << " Tstart="
	<<T1<<" Tend="<<T2<<" Apas=" <<Apas << " Tpas="<<Tpas<<" rep="<<rep
	<<" generate_squares="<<generate_squares<<" use_squares="<<
	use_squares <<endl;
	cout<<"Please confirm the parameters and press enter to start"<<endl;
	cout<<"As for other parameters pls check config file"<<endl;
	std::cin.get();
	if (generate_squares){
		s_squars="squares_"+std::string(argv[3])+".txt";
	} else {
		//s_squars;
	}
	Run run(s, s_squars);
	
	if (generate_squares){
		cout<<"Generating squares"<<endl;
		double resBL=run.multip_test(T1,A1,D,rep);
		double resBR=run.multip_test(T2,A1,D,rep);
		double resTL=run.multip_test(T1,A2,D,rep);
		double resTR=run.multip_test(T2,A2,D,rep);
		output<<T1<<" "<<A1<<" "<<resTL<<endl;
		output<<T2<<" "<<A1<<" "<<resTR<<endl;
		output<<T1<<" "<<A2<<" "<<resBL<<endl;
		output<<T2<<" "<<A2<<" "<<resBR<<endl;
		run.toDet_.push_back(std::make_tuple(T1,T2,A1,A2,resTL,resTR,
		resBL,resBR,D,rep));
		output.close();
		while ((run.toDet_.size()!=0)&&(run.toDet_.size()<=300)){
			run.dynamic_analyse(run.toDet_[0],s_output);
			run.toDet_.erase(run.toDet_.begin());
			cout<<"Queue size: "<<run.toDet_.size()<<" "<<
			run.notDet_.size()<<endl;
		}
		run.squars_.close();
	}
	
	if (use_squares){//read_squares
		std::ifstream squars;
		squars.open(s_squars,std::ios::in);
		int area=0;
		int a1=0,a2=0,t1=0,t2=0;
		while(squars>>t1>>t2>>a1>>a2>>tmpres){
			sqs.push_back(std::make_tuple(t1,t2,a1,a2,tmpres));
			area+=(t2-t1)*(a2-a1);
		}
		cout<<area<<endl;
		squars.close();
	}
	
	output.clear();
	output.open(s_output.c_str(),std::ios::app);
	for (int T=T1; T<=T2;T = T+ Tpas){
		for (double A=A1; A<=A2; A=A+Tpas){
			int mark=0;
			if (use_squares){
				for (unsigned int i=0;i<sqs.size();i++){
					if ((T>=std::get<0>(sqs[i]))&&
						(T<=std::get<1>(sqs[i]))&&
						(A>=std::get<2>(sqs[i]))&&
						(A<=std::get<3>(sqs[i]))){
						output<<T <<" "<<A<<" "<<
						 std::get<4>(sqs[i])<< endl;
						mark=1;
						break;
					}
				}
			}
			if ((mark==0)||(!use_squares)){
				cout<<T<<" "<<A<<endl;
				output<<T<<" "<<A<<" "<<
				run.multip_test(T,A,D,rep)<<endl;
			}
		}
	}
	return 0;
}
