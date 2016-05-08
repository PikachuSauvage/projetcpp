// ===========================================================================
//                                  Includes
// ===========================================================================
#ifndef RUN_H
#define RUN_H
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
class Run{
// =========================================================================
//                              Protected Methods
// ========================================================================= 
	protected:
		int W_ =0;
		int H_ =0; 
		int max_size_=300;
		double A_percent_ = 0;
		double RAA_=0;
		double RAB_=0;
		double RBB_=0;
		double RBC_=0;
		double Wmin_=0;
		double Pdie_=0;
		double Pmut_=0;
		double odeStep_=0;
		double odeN_=1;
	public:
		std::vector<std::tuple<int, int, int, int, int, int, int, int, 
		double, int>> toDet_;
		std::vector<std::tuple<int, int, int, int>> notDet_;
		std::ofstream squars_;
// =========================================================================
//                               Constructors
// =========================================================================
		Run(std::string s, std::string s_squars);
// ===========================================================================
//                                 Destructor
// ===========================================================================
		~Run();
// =========================================================================
//                              Public Methods
// =========================================================================
		double single_test(int T, int A, double D);
		void dynamic_analyse(std::tuple<int,int,int,int,double,double,
					double,double,double,int> tupin, std::string s_output);
		double multip_test(int T, int A, double D, int Rep);

// =========================================================================
//                                  Getters
// =========================================================================

// =========================================================================
//                                  Setters
// =========================================================================

};
#endif // RUN_H
