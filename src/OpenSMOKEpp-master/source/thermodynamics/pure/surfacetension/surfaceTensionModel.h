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
 \class SurfaceTensionModel
 \brief Class to determine the components' surface tension \f$\sigma_{i}\f$ by use of polynomial and coefficients given in VDI Waermeatlas
 
 This class is used to determine the components' surface tension of each component using a polynomial and coefficients given in VDI Waermeatlas.
 */
#ifndef _SURFACETENSIONMODEL_H_
#define _SURFACETENSIONMODEL_H_

#include <iostream>
#include <cmath>

using namespace std;

class SurfaceTensionModel {
    
private:
    double* _coefficients;
    double _Tc;
public:
    /** \brief Constructor with coefficients given in configuration file and critical temperature for evaluating polynomial for each component
     @param coefficients coefficients for polynomial
     @param Tc critical temperature of related component
     
     */
    SurfaceTensionModel(double* coefficients, double Tc);
    /** \brief Function to determine surface tension \f$\sigma_i\f$ at given temperature
     @param T temperature
     
     \f[
     \sigma_i=10^{-3}\cdot a_i\left(1-\frac{T}{T_c}\right)^{b_i}
     \f]
     */
    double sigma(double T);
};

#include "surfaceTensionModel.hpp"

#endif
