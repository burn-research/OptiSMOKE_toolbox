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
 \class EtaLPrausnitz
 \brief Class to determine the liquid components' dynamic viscosity \f$\eta_{i}^l\f$ by use of polynomial and coefficients given in Poling, Prausnitz, O'Connell: The Properties of Gases and Liquids
 
 This class is used to determine the liquid components' dynamic viscosity of each component using a polynomial and coefficients given in Poling, Prausnitz, O'Connell: The Properties of Gases and Liquids.
 */
#ifndef _ETALPRAUSNITZ_H_
#define _ETALPRAUSNITZ_H_

#include <cmath>
#include "viscosityModel.h"

#include <iostream>


using namespace std;

class EtaLPrausnitz : public ViscosityModel{
private:
    double _Tc;
public:
    /** \brief Constructor with coefficients given in configuration file and critical temperature for evaluating polynomial for liquid phase components
     @param coefficients coefficients for polynomial
     @param Tc critical temperature of related component
     
     */
    EtaLPrausnitz(double* coefficients, double Tc);
    /** \brief Function to determine dynamic viscosity \f$\eta_{i}^l\f$ at given temperature
     @param T temperature
     
     \f[
     \log\eta_{i}^l=a_i+b_i\left(1-\frac{T}{c_i}\right)^{2/7}
     \f]
     */
    virtual double eta(double T);
};

#include "etaLPrausnitz.hpp"
#endif
