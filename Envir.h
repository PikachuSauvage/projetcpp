// ===========================================================================
//                                  Includes
// ===========================================================================
#ifndef ENVIR_H
#define ENVIR_H

#include "Ecoli.h"
#include "vector"
class Envir{
// =========================================================================
//                              Protected Methods
// ========================================================================= 
	protected:
		double* Qa_;
		double* Qb_;
		double* Qc_;
		double* F_;
		Ecoli** indiv_;
		unsigned int* reproduced_;
		double D_;
		double RAA_;
		double RAB_;
		double RBB_;
		double RBC_;
		double Wmin_;
		int status_;
		unsigned int W_;
		unsigned int H_;
		void mutation(double prob, int posi);
	public:
// =========================================================================
//                               Constructors
// =========================================================================
		Envir(unsigned int W, unsigned int H, double D, double A_homo,
			  double RAA, double RAB, double RBB, double RBC, double Wmin);
// ===========================================================================
//                                 Destructor
// ===========================================================================
		~Envir();
// =========================================================================
//                              Public Methods
// =========================================================================
		void updateMetab(double p, double n);
		void updateFitness();
		void diffuse();
		void plsDie(double prob);
		void toSurvive(double prob);
		void print();
		void run(int TMAX);
		void init(int W, int H, double percentageGA);
		void reinit();
		void refeed(double A);
		void result();
// =========================================================================
//                                  Getters
// =========================================================================
		double getQa(unsigned int x, unsigned int y);
		double getQb(unsigned int x, unsigned int y);
		double getQc(unsigned int x, unsigned int y);
		int getStatus();
		//unsigned int getRep(unsigned int x, unsigned int y);
		Ecoli* getEcoli(unsigned int x, unsigned int y);
// =========================================================================
//                                  Setters
// =========================================================================
		void setQa(unsigned int x, unsigned int y, double newValue);
		void setQb(unsigned int x, unsigned int y, double newValue);
		void setQc(unsigned int x, unsigned int y, double newValue);
};
#endif // ENVIR_H
