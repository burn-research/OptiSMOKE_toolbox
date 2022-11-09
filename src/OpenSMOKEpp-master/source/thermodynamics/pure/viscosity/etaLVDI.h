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
 \class EtaLVDI
 \brief Class to determine the liquid components' dynamic viscosity \f$\eta_{i}^l\f$ by use of polynomial and coefficients given in VDI Waermeatlas
 
 This class is used to determine the liquid components' dynamic viscosity of each component using a polynomial and coefficients given in VDI Waermeatlas.
 */

#ifndef _ETALVDI_H_
#define _ETALVDI_H_

#include <cmath>
#include "viscosityModel.h"

#include <iostream>


using namespace std;

class EtaLVDI : public ViscosityModel{
    
public:
    /** \brief Constructor with coefficients given in configuration file for evaluating polynomial for liquid phase components
     @param coefficients coefficients for polynomial
     
     */
    EtaLVDI(double* coefficients);
    /** \brief Function to determine dynamic viscosity \f$\eta_{i}^l\f$ at given temperature
     @param T temperature
     
     \f[
     \ln\eta_i^l=a_i+b_i/T+c_iT+d_iT^2+e_iT^3
     \f]
     */
    virtual double eta(double T);
};

#include "etaLVDI.hpp"
#endif