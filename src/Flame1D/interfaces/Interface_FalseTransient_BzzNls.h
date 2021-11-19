/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ Suite.                               |
|                                                                         |
|   Copyright(C) 2015, 2014, 2013  Alberto Cuoci                          |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

class OpenSMOKE_Flame1D_MyFalseTransientSystem_Premixed_FLAMESPEED : public BzzMyNonLinearSystemSparseObject
{
public:

	void assign(OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D *flame)
	{
		ptFlame = flame;
	}
	
	void SetInitialConditions(const BzzVector& yInitial)
	{
		yInitial_ = yInitial;
	}

	void SetDifferentialAlgebraic(const BzzVectorInt& indices_differential_algebraic)
	{
		indices_differential_algebraic_ = indices_differential_algebraic;
	}
	
	void SetTimeStep(const double deltat)
	{
		deltat_ = deltat;
	}

	virtual void GetResiduals(BzzVector &y, BzzVector &f)
	{
		double* pty = y.GetHandle();
		double* ptf = f.GetHandle();

		ptFlame->Equations(0., pty, ptf);

		for (int i = 1; i <= y.Size(); i++)
		if (indices_differential_algebraic_[i] == 1)
			f[i] = y[i] - yInitial_[i] - f[i] * deltat_;
	}

	virtual void ObjectBzzPrint(void)
	{
	}

public:

	double deltat() const { return deltat_; }
	const BzzVector& InitialConditions() const { return yInitial_; }
	BzzVector& InitialConditions() { return yInitial_; }

private:

	OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D *ptFlame;
	double deltat_;
	BzzVector yInitial_;
	BzzVectorInt indices_differential_algebraic_;
};

void FalseTransientPrint(BzzVector &y, double t)
{
	//flame_premixed->Print(t, y.GetHandle());
}

#include "math/native-nls-solvers/interfaces-false-transient/Band_BzzNlsFalseTransient.h"

