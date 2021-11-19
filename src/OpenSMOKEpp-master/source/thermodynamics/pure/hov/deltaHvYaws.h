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
 \class DeltaHvYaws
 \brief Class to determine the components' heat of vaporization \f$\Delta h_{v,i}\f$ by use of polynomial and coefficients given in Yaws: Chemical Properties Handbook
 
 This class is used to determine the components' heat of vaporization using a polynomial and coefficients given in Yaws: Chemical Properties Handbook.
 */

#ifndef _DELTAHVYAWS_H_
#define _DELTAHVYAWS_H_

#include "heatOfVaporizationModel.h"
#include <iostream>

using namespace std;

class DeltaHvYaws : public HeatOfVaporizationModel {
private:
    double _M;
public:
    /** \brief Constructor with coefficients given in configuration file, critical temperature and molar weight for evaluating polynomial for each component
     @param coefficients coefficients for polynomial
     @param Tc critical temperature of related species
     @param M molar weight of related species
     
     */
    DeltaHvYaws(double* coefficients, double Tc, double M);
    /** \brief Function to determine heat of vaporization \f$\Delta h_{v,i}\f$ at given temperature
     @param T temperature
     
     \f[
     \Delta h_{v,i}=\Delta h_{v,i}=a_i\left(1-\frac{T}{T_c}\right)^b_i\frac{10^3}{M_i}
     \f]
     */
    virtual double deltaHv(double T);
};

#include "deltaHvYaws.hpp"
#endif
