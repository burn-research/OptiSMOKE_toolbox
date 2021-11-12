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
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunmatrix/sunmatrix_band.h>  /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver       */

#include <ida/ida.h>
#include <ida/ida_direct.h>            /* access to IDADls interface           */

#include <kinsol/kinsol.h>
#include <kinsol/kinsol_direct.h>      /* access to KINDls interface      */

#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>

static int check_flag(void *flagvalue, char *funcname, int opt)
{
	// 0. Check if SUNDIALS function returned NULL pointer - no memory allocated 
	if (opt == 0 && flagvalue == NULL) 
	{
		fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
		return(1);
	}
	// 1. Check if flag < 0
	else if (opt == 1) 
	{
		
		int *errflag = (int *)flagvalue;
		if (*errflag < 0) 
		{
			fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
			return(1);
		}
	}
	// 2. Check if function returned NULL pointer - no memory allocated
	else if (opt == 2 && flagvalue == NULL) 
	{
		fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
		return(1);
	}

	return(0);
}

realtype N_SumAbs(N_Vector x)
{
	realtype *xd;
	xd = NULL;

	long int N = NV_LENGTH_S(x);
	xd = NV_DATA_S(x);

	realtype sum = 0.;
	for (long int i = 0; i < N; i++)
		sum += std::fabs(xd[i]);

	return(sum);
}

realtype N_Norm2(N_Vector x)
{
	realtype *xd;
	xd = NULL;

	long int N = NV_LENGTH_S(x);
	xd = NV_DATA_S(x);

	realtype sum = 0.;
	for (long int i = 0; i < N; i++)
		sum += xd[i] * xd[i];

	return(std::sqrt(sum));
}