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

typedef struct
{
	double  deltat;
	double* yInitial;
} *FalseTransient_UserData;

static int kinsol_equations_false_transient(N_Vector u, N_Vector f, void *user_data)
{
	realtype *pt_y = NV_DATA_S(u);
	realtype *pt_res = NV_DATA_S(f);
	FalseTransient_UserData data = (FalseTransient_UserData)user_data;

	flame_premixed->Equations(0., pt_y, pt_res);

	for (int i = 0; i < flame_premixed->NumberOfEquations(); i++)
	if (flame_premixed->id_equations()[i] == true)
		pt_res[i] = pt_y[i] - data->yInitial[i] - pt_res[i] * data->deltat;

	return 0;
}

#include "math/native-nls-solvers/interfaces-false-transient/Band_KinSolFalseTransient.h"

