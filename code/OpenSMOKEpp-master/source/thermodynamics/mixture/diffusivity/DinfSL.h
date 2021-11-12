/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alessandro Stagni <alessandro.stagni@polimi.it>               |
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
|   Copyright(C) 2016, 2015 Alessandro Stagni & Alberto Cuoci             |
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
 \class DinfSL
 \brief Class to calculate infinite-dilution diffusion coefficients via Siddiqi-Lucas rule
 
 Class to calculate infinite-dilution diffusion coefficients via Siddiqi-Lucas rule
 */

#ifndef OPENSMOKEPP_DINFSL_H
#define OPENSMOKEPP_DINFSL_H

#include "DinfModel.h"
#include "thermodynamics/mixture/mixtureL/mixtureL.h"

class DinfSL : public DinfModel
{
	public:

		DinfSL(const unsigned int NS, const std::vector<double>& MW,
		       const std::vector<double>& Tnbp, std::vector<DensityModel*>& rhoL);
  
  
		void UpdateInfiniteDiffusivities(const std::vector<double>& etaL, const double T, const double P);
  
	private:

	const std::vector<double>& MW_;
	const std::vector<double>& Tnbp_;
	std::vector<DensityModel*>& rhoL_;
};

#include "DinfSL.hpp"

#endif	// OPENSMOKEPP_DINFSL_H
