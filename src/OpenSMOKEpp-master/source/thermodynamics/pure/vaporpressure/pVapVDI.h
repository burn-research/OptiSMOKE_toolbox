/*-----------------------------------------------------------------------*\
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
|	License                                                           |
|                                                                         |
|   Copyright(C) 2016, 2015 Alberto Cuoci                                 |
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

/**
 \class PVapVDI
 \brief Class to determine the components' vapor pressure \f$p_{vap,i}\f$ by use of polynomial and coefficients given in VDI Waermeatlas
 
 This class is used to determine the liquid components' vapor pressure of each component using a polynomial and coefficients given in VDI Waermeatlas.
 */
#ifndef _PVAPVDI_H_
#define _PVAPVDI_H_


#include "vaporPressureModel.h"
#include <cmath>
#include <iostream>

using namespace std;

class PVapVDI : public VaporPressureModel {
private:
    double _pc;
public:
    /** \brief Constructor with coefficients given in configuration file, critical temperature and critical pressure for evaluating polynomial for liquid phase components
     @param coefficients coefficients for polynomial
     @param Tc critical temperature of related component
     @param pc critical pressure of related component
     
     */
    PVapVDI(double* coefficients, double Tc, double pc);
    /** \brief Function to determine vapor pressure \f$p_{vap,i}\f$ at given temperature
     @param T temperature
     
     \f[
     \ln\frac{p_{vap}}{p_c}=\frac{1}{T_r}\left(a_i(1-T_r)+b_i(1-T_r)^{1.5}+c_i(1-T_r)^3+d_i(1-T_r)^6\right)
     \f]
     */
    virtual double pVap(double T);
    virtual double dpVapdT(double T);
};

#include "pVapVDI.hpp"

#endif
