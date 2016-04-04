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
		//unsigned int* reproduced_;
		double D_;
		double RAA_;
		double RAB_;
		double RBB_;
		double RBC_;
		double Wmin_;
		unsigned int W_;
		unsigned int H_;
		void mutation(double prob, int posi, char origin);
	public:
// =========================================================================
//                               Constructors
// =========================================================================
		Envir(unsigned int W, unsigned int H, double D, double A_homo,
			  double RAA, double RAB, double RBB, double RBC, double Wmin);
// ===========================================================================
//                                 Destructor
// ===========================================================================

// =========================================================================
//                              Public Methods
// =========================================================================
		void updateMetab();
		void updateFitness();
		void diffuse();
		void plsDie(double prob);
		void toSurvive(double prob);
		void print();
		void run(int TMAX);
// =========================================================================
//                                  Getters
// =========================================================================
		double getQa(unsigned int x, unsigned int y);
		double getQb(unsigned int x, unsigned int y);
		double getQc(unsigned int x, unsigned int y);
		//unsigned int getRep(unsigned int x, unsigned int y);
		Ecoli* getEcoli(unsigned int x, unsigned int y);
// =========================================================================
//                                  Setters
// =========================================================================
		void setQa(unsigned int x, unsigned int y, double newValue);
		void setQb(unsigned int x, unsigned int y, double newValue);
		void setQc(unsigned int x, unsigned int y, double newValue);
		void init(int W, int H, double percentageGA);
};
#endif // ENVIR_H
