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

class OpenSMOKE_Flame1D_MyDaeSystem_Premixed_FLAMESPEED : public BzzDaeSystemObject
{
public:

	OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D *ptFlame;

	void assign(OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D *flame)
	{
		ptFlame = flame;
	}

	virtual void GetSystemFunctions(BzzVector &x, double t, BzzVector &f)
	{
		double* ptx = x.GetHandle();
		double* ptf = f.GetHandle();

		ptFlame->Equations(t, ptx, ptf);
	}
	
	virtual void ObjectBzzPrint(void)
	{
	}
};

void DaePrint(BzzVector &y, double t)
{
	//flame_premixed->Print(t, y.GetHandle());
}

#include "math/native-dae-solvers/interfaces/TridiagonalBlock_BzzDae.h"
