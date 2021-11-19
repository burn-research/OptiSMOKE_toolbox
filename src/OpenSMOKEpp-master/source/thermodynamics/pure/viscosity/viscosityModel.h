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
 \class ViscosityModel
 \brief Abstract base class to determine dynamic viscosity \f$\eta_i\f$
 
 This class is used as base class containing the coefficients of each species to determine the dynamic viscosity of each phase component.
 */
#ifndef _VISCOSITYMODEL_H_
#define _VISCOSITYMODEL_H_

#include <iostream>


using namespace std;

class ViscosityModel {
    
private:
    double* _coefficients;
    friend class EtaLVDI;
    friend class EtaLYaws;
    friend class EtaLPrausnitz;
    friend class EtaGVDI;
    friend class EtaGYaws;
public:
    /** \brief Constructor with coefficients given in configuration file for evaluating polynomials in derived classes
     @param coefficients coefficients for polynomial
     
     */
    ViscosityModel(double* coefficients);
    ~ViscosityModel();
    /** \brief Pure virtual function declared to determine dynamic viscosity \f$\eta_{i}\f$ at given temperature (has to be declared in derived classes)
     @param T temperature
     */
    virtual double eta(double T) = 0;
};

#include "viscosityModel.hpp"

#endif
