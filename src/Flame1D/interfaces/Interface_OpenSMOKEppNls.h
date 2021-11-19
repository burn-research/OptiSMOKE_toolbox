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

#include "math/native-nls-solvers/NonLinearSystemSolver"

class OpenSMOKE_Flame1D_MyNlsSystem_OpenSMOKEpp_Premixed_FLAMESPEED
{
public:

	OpenSMOKE_Flame1D_MyNlsSystem_OpenSMOKEpp_Premixed_FLAMESPEED() 
	{
	};

	void assign(OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D *flame)
	{
		ptFlame = flame;
	}

private:

	OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D *ptFlame;

protected:

	unsigned int ne_;

	void MemoryAllocation()
	{
	}

	virtual void Equations(const Eigen::VectorXd &x, Eigen::VectorXd &f)
	{
		ptFlame->Equations(0., x.data(), f.data());
	}

	void Jacobian(const Eigen::VectorXd &y, Eigen::MatrixXd &J)
	{
	};

	void Print(const int call, const double t, const double phiW, const Eigen::VectorXd &x, const Eigen::VectorXd &f)
	{
		//flame_premixed->Print(x.data(), phiW);
	}
};

#include "math/native-nls-solvers/interfaces/Band_OpenSMOKEppNls.h"
#include "math/native-nls-solvers/interfaces/Sparse_OpenSMOKEppNls.h"

