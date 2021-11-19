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

#include <boost/timer/timer.hpp>

typedef struct
{
	N_Vector J;
	N_Vector invJ;
} *IDAUserData;


int ida_equations(realtype t, N_Vector y, N_Vector yp, N_Vector res, void *user_data)
{
	realtype *pt_y = NV_DATA_S(y);
	realtype *pt_res = NV_DATA_S(res);
	realtype *pt_yp = NV_DATA_S(yp);

	flame_premixed->Equations(t, pt_y, pt_res);
	flame_premixed->CorrectDifferentialEquations(pt_yp, pt_res);

	return 0;
}

int ida_initial_derivatives(realtype t, N_Vector y, N_Vector yp, void *user_data)
{
	realtype *pt_y = NV_DATA_S(y);
	realtype *pt_yp = NV_DATA_S(yp);

	flame_premixed->Equations(t, pt_y, pt_yp);
	flame_premixed->CorrectAlgebraicEquations(pt_yp);
	
	return 0;
}

int ida_preconditioner_setup(realtype t, N_Vector y, N_Vector yp, N_Vector rr, realtype c_j, void *user_data)
{
	IDAUserData data;
	data = (IDAUserData)user_data;

	realtype *pt_y = NV_DATA_S(y);
	realtype *pt_J = NV_DATA_S(data->J);
	
	flame_premixed->DiagonalJacobian(t, pt_y, pt_J);
	flame_premixed->DiagonalJacobianForIDA(c_j, pt_J);

	realtype *pt_invJ = NV_DATA_S(data->invJ);
	for (int i = 0; i < NV_LENGTH_S(y); i++)
		pt_invJ[i] = 1. / pt_J[i];

	return(0);
}

int ida_preconditioner_solution(realtype t, N_Vector y, N_Vector yp, N_Vector rr, N_Vector rvec, N_Vector zvec, realtype c_j, realtype delta, void *user_data)
{
	IDAUserData data;
	data = (IDAUserData)user_data;

	N_VProd(data->invJ, rvec, zvec);

	return(0);
}

int ida_print_solution(realtype t, N_Vector y)
{
	double* ydata = N_VGetArrayPointer(y);
	flame_premixed->Print(t, ydata);
	return 0;
}

#include "math/native-dae-solvers/interfaces/Band_Ida.h"
