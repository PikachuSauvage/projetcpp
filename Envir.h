// ===========================================================================
//                                  Includes
// ===========================================================================
#ifndef ENVIR_H
#define ENVIR_H

#include "Ecoli.h"
class Envir{
// =========================================================================
//                              Protected Methods
// ========================================================================= 
	protected:
		double* Qa_;
		double* Qb_;
		double* Qc_;
		Ecoli** indiv_;
		unsigned int* reproduced_;
		double D_;
		unsigned int W_;
		unsigned int H_;
	public:
// =========================================================================
//                               Constructors
// =========================================================================
		Envir(unsigned int W, unsigned int H, double D, double A_MAX);
// ===========================================================================
//                                 Destructor
// ===========================================================================

// =========================================================================
//                              Public Methods
// =========================================================================

// =========================================================================
//                                  Getters
// =========================================================================
		double getQa(unsigned int x, unsigned int y);
		double getQb(unsigned int x, unsigned int y);
		double getQc(unsigned int x, unsigned int y);
		unsigned int getRep(unsigned int x, unsigned int y);
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
