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
 \class RhoLPengRobinson
 \brief Class to determine the liquid components' density \f$\rho_{i}^l\f$ by use of Peng Robinson equation of state
 
 This class is used to determine the liquid components' density of each component using a polynomial and coefficients given in Yaws: Chemical Properties Handbook.
 */
#ifndef _RHOLPENGROBINSON_H_
#define _RHOLPENGROBINSON_H_

#include "densityLModel.h"
#include <iostream>
#include "../eos/pengrobinson.h"


using namespace std;

class RhoLPengRobinson : public DensityModel {
    
public:

    RhoLPengRobinson(const double Tc,
                const double Pc, const double omega,
                const double MW);
    ~RhoLPengRobinson();

    virtual double rho(double T, double P);
    
    pengrobinson* pr;
    
private:
  double Tc_;
  double Pc_;
  double omega_;
  double MW_;
};

#include "rhoLPengRobinson.hpp"

#endif
