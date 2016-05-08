#include "Envir.h"
#include "Ecoli.h"
#include "Run.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

using std::cout;
using std::endl;

int main(int argc, char *argv[]){
	std::string s_output;
	std::string s_squars;
	std::ofstream output;
	double curr_res=0;
	double res;
	int rep=5;
	int T1=1;
	int T2=1500;
	int A1=0;
	int A2=50;
	double D=0.1;
	int a1=0,a2=0,t1=0,t2=0;
	double tmpres=0;
	unsigned long seed = unsigned (std::time(0));
	std::srand (seed);
	cout << "Hello projet C++!" << endl;
	std::string s("config.txt");
	s_output= "output"+ std::string(argv[1]) + "_"+std::string(argv[2])
				+".txt";
	s_squars="squars_"+std::string(argv[3])+"";
	Run run(s, s_squars);
	output.open(s_output.c_str());
	if (!output){
		cout<<"File not opened!"<<endl;
	}
	
	double resTL=run.multip_test(T1,A1,D,rep);
	double resTR=run.multip_test(T2,A1,D,rep);
	double resBL=run.multip_test(T1,A2,D,rep);
	double resBR=run.multip_test(T2,A2,D,rep);
	output<<T1<<" "<<A1<<" "<<resTL<<endl;
	output<<T2<<" "<<A1<<" "<<resTR<<endl;
	output<<T1<<" "<<A2<<" "<<resBL<<endl;
	output<<T2<<" "<<A2<<" "<<resBR<<endl;
	run.toDet_.push_back(std::make_tuple(1,1500,0,50,resTL,resTR,resBL,
	resBR,D,rep));
	output.close();
	while ((run.toDet_.size()!=0)&&(run.toDet_.size()<=300)){
		run.dynamic_analyse(run.toDet_[0],s_output);
		run.toDet_.erase(run.toDet_.begin());
		cout<<"Queue size: "<<run.toDet_.size()<<" "<<run.notDet_.size()<<endl;
	}
	run.squars_.close();
	
	std::ifstream squars;
	squars.open(s_squars,std::ios::in);
	std::vector<std::tuple<int, int, int, int,double>> sqs;
	int area=0;
	while(squars>>t1>>t2>>a1>>a2>>tmpres){
		sqs.push_back(std::make_tuple(t1,t2,a1,a2,tmpres));
		area+=(t2-t1)*(a2-a1);
	}
	cout<<area<<endl;
	//for (double d=0.1; d<=0.1; d+=0.01){
		//for (int T=atoi(argv[1]); T<=atoi(argv[2]);T = T+ 1){
			//for (int A=0; A<=50; A=A+1){
				//int mark=0;
				//for (int i=0;i<sqs.size();i++){
					////cout<<std::get<0>(sqs[i])<<" "<<std::get<1>(sqs[i])<<" "<<endl;
					//if ((T>=std::get<0>(sqs[i]))&&(T<=std::get<1>(sqs[i]))&&
						//(A>=std::get<2>(sqs[i]))&&(A<=std::get<3>(sqs[i]))){
						//output<<T <<" "<<A<<" "<< std::get<4>(sqs[i])<< endl;
						//mark=1;
						//break;
					//}
				//}
				//if (mark==0){
					//output<<T<<" "<<A<<" 3"<<endl;
					////cout<<T<<" "<<A<<" "<<run.multip_test(T,A,d,rep)<<endl;
				//}
			//}
		//}
	//}
	//for (double d=0.1; d<=0.1; d+=0.01){
		//for (int T=atoi(argv[1]); T<=atoi(argv[2]);T = T+ atoi(argv[3])){
			//for (int A=0; A<=50; A=A+1){
				//res=0;
				//for (int n=0;n<rep;n++){
					//curr_res=run.single_test(T,A,d);
					//cout<< T << " " << A << " " << curr_res << endl;
					//res+=curr_res;
				//}
				//output<< T <<" "<<A<<" "<<res<<endl;
			//}
		//}
	//}
	squars.close();
	return 0;
}
