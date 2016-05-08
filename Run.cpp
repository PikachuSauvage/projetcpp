#include "Envir.h"
#include "Run.h"
#include <cmath>
#include <cstring>
#include <iostream>
using std::cout;
using std::endl;
//using std::vector;
//using std::tuple;


Run::Run(std::string s, std::string s_squars){
	std::ifstream config;
	config.open(s.c_str(),std::ios::in);
	squars_.open(s_squars.c_str(),std::ios::out);
	std::string key;
	if (config.is_open()){
		std::getline(config, key, '=');
		config>>W_;
		std::getline(config, key, '=');
		config>>H_;
		std::getline(config, key, '=');
		config>>A_percent_;
		std::getline(config, key, '=');
		config>>RAA_;
		std::getline(config, key, '=');
		config>>RAB_;
		std::getline(config, key, '=');
		config>>RBB_;
		std::getline(config, key, '=');
		config>>RBC_;
		std::getline(config, key, '=');
		config>>Wmin_;
		std::getline(config, key, '=');
		config>>Pdie_;
		std::getline(config, key, '=');
		config>>Pmut_;
		std::getline(config, key, '=');
		config>>odeN_;
		std::getline(config, key, '=');
		config>>odeStep_;
	}
	config.close();
}

Run::~Run(){}

double Run::single_test(int T, int A, double D){
	int refeed_count=0;
	Envir test = Envir(W_,H_,D,A,RAA_,RAB_,RBB_,RBC_,Wmin_);
	test.init(W_,H_,A_percent_);
	//cout<<W_<< " "<< H_<< " "<< A_percent_<< " "<< odeStep_ << endl;
	test.updateMetab(odeStep_, odeN_);
	//test.print();
	for (int i=0; i<10000; i++){
		test.diffuse();
		test.updateFitness();
		test.toSurvive(Pmut_);
		if (test.getStatus()==0)
			break;
		test.plsDie(Pdie_);
		test.updateMetab(odeStep_, odeN_);
		test.reinit();
		refeed_count++;
		if (refeed_count==T){
			refeed_count=0;
			test.refeed(A);
		}
	}
	test.result();
	if (test.getStatus()==3)
		test.print();
	return test.getStatus();
}

void Run::dynamic_analyse(std::tuple<int, int ,int ,int,double,double,
		double,double, double, int> tupin, std::string s_output){
	int T1=std::get<0>(tupin);
	int T2=std::get<1>(tupin);
	int A1=std::get<2>(tupin);
	int A2=std::get<3>(tupin);
	double resTL=std::get<4>(tupin);
	double resTR=std::get<5>(tupin);
	double resBL=std::get<6>(tupin);
	double resBR=std::get<7>(tupin);
	double D=std::get<8>(tupin);
	int rep=std::get<9>(tupin);
	std::ofstream output;
	int Tn=(T1+T2)/2;
	int An=(A1+A2)/2;
	if ((Tn==T1)||(Tn==T2)||(An==A1)||(An==A2))
		return;
	output.clear();
	output.open(s_output.c_str(),std::ios::app);
	if (!output){
		cout<<"File not opened!"<<endl;
	}
	double resTM=multip_test(Tn,A1,D,rep);
	double resMM=multip_test(Tn,An,D,rep);
	double resML=multip_test(T1,An,D,rep);
	double resMR=multip_test(T2,An,D,rep);
	double resBM=multip_test(Tn,A2,D,rep);
	output<<Tn<<" "<<A1<<" "<<resTM<<endl;
	output<<Tn<<" "<<An<<" "<<resMM<<endl;
	output<<T1<<" "<<An<<" "<<resML<<endl;
	output<<T2<<" "<<An<<" "<<resMR<<endl;
	output<<Tn<<" "<<A2<<" "<<resBM<<endl;
	if ((resTL==resTM)&&(resTM==resML)&&(resML==resMM)){
		cout<<"NOT DET:"<<T1<<" "<<Tn<<" "<<A1<<" "<<An<<endl;
		squars_<<T1<<" "<< Tn <<" "<< A1 <<" "<<An<<" "<< resTL << endl;
		notDet_.push_back(std::make_tuple(T1,Tn,A1,An));
	} else {
		if (toDet_.size()>=max_size_){
			//cout<<"Terrain pull up"<<endl;
			return;
		} else {
			cout<<"Pushing:"<<T1<<" "<<Tn<<" "<<A1<<" "<<An<<endl;
			toDet_.push_back(std::make_tuple(T1,Tn,A1,An,resTL, resTM, resML, resMM,D,rep));
		}
	}
	if ((resTM==resTR)&&(resMM==resMR)&&(resTM=resMR)){
		cout<<"Not Det:"<<Tn<<" "<<T2<<" "<<A1<<" "<<An<<endl;
		squars_<<Tn<<" "<<T2<<" "<<A1<<" "<<An<<" "<< resTM << endl;
		notDet_.push_back(std::make_tuple(Tn,T2,A1,An));
	} else {
		if (toDet_.size()>=max_size_){
			return;
		} else {
			cout<<"PUSHING:"<<Tn<<" "<<T2<<" "<<A1<<" "<<An<<endl;
			toDet_.push_back(std::make_tuple(Tn,T2,A1,An,resTM, resTR, resMM, resMR,D,rep));
		}
	}
	if ((resML==resMM)&&(resBL==resBM)&&(resML==resBM)){
		cout<<"not DET:"<<T1<<" "<<Tn<<" "<<An<<" "<<A2<<endl;
		squars_<<T1<<" "<<Tn<<" "<<An<<" "<<A2<< " "<< resML<< endl;
		notDet_.push_back(std::make_tuple(T1,Tn,An,A2));
	} else {
		if (toDet_.size()>=max_size_){
			return;
		} else {
			cout<<"PushinG:"<<T1<<" "<<Tn<<" "<<An<<" "<<A2<<endl;
			toDet_.push_back(std::make_tuple(T1,Tn,An,A2,resML,resMM,resBL,resBM,D,rep));
		}
	}
	if ((resMM==resMR)&&(resBM==resBR)&&(resMM==resBR)){
		cout<<"not det4:"<<Tn<<" "<<T2<<" "<<An<<" "<<A2<<endl;
		squars_<<Tn<<" "<<T2<<" "<<An<<" "<<A2<<" "<<resMM<<endl;
		notDet_.push_back(std::make_tuple(Tn,T2,An,A2));
	} else {
		if (toDet_.size()>=max_size_){
			return;
		} else {
			cout<<"pushing 4:"<<Tn<<" "<<T2<<" "<<An<<" "<<A2<<endl;
			toDet_.push_back(std::make_tuple(Tn,T2,An,A2,resMM,resMR,resBM,resBR,D,rep));
		}
	}
	output.close();
}

double Run::multip_test(int T, int A, double D, int Rep){
	double curr_res=0;
	double res=0;
	cout<<"da wang jiao wo lai xun shan: "<<T<<" "<<A<<endl;
	for (int i=0; i<Rep; i++){
		curr_res=single_test(T,A,D);
		res+=curr_res;
	}
	res=res/Rep;
	cout<<res<<endl;
	return res;
}
