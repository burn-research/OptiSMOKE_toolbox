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

static void DaspkEquations(double *t, double *y, double *dy, double *cj, double *delta, int *ires, double *rpar, int *ipar)
{
	flame_premixed->Equations(*t, y, delta);
	flame_premixed->CorrectDifferentialEquations(dy, delta);
}

static void DaspkAlgebraicDifferentialVector(double* v)
{
	// Returns 1. if differential, 0. if algebraic
	flame_premixed->AlgebraicDifferentialVector(v);
}

void DaspkInitialDerivatives(double t, double *y, double *yp)
{
	flame_premixed->Equations(t, y, yp);
	flame_premixed->CorrectAlgebraicEquations(yp);
}

static void DaspkAnalyticalJacobian(double *x, double *y, double *dy, double *pd, double *cj, double *rpar, int *ipar)
{
}

static void DaspkKrylovSolver(int *n, double *x, double *y, double *dy, double *savr, double *wk, double *cj, double *wght, double *wp, int *iwp, double *b, double *eplin, int *ier, double *rpar, int *ipar)
{
}

static void DaspkPrintSolution(double *x, double *y)
{
	flame_premixed->Print(*x, y);
}

#include "math/native-dae-solvers/interfaces/Band_Daspk.h"
