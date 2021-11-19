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
 \class FugacityGPengRobinson_mix
 \brief Class to calculate gas fugacity with Peng-Robinson EoS
 
 Class to calculate gas fugacity with Peng-Robinson EoS
 */

#ifndef OPENSMOKEPP_FUGACITYG_PENGROBINSON_MIX_H
#define OPENSMOKEPP_FUGACITYG_PENGROBINSON_MIX_H

class FugacityGPengRobinson_mix : public FugacityModel_mix
{
	public:

		FugacityGPengRobinson_mix(	const unsigned int NS, 
									const std::vector<double>& Tc, const std::vector<double>& Pc,
									const std::vector<double>& omega, const std::vector<double>& MW);

		~FugacityGPengRobinson_mix();

		std::vector<double> fugacity(const double T, const double P, const std::vector<double>& x);
		std::vector<double> fugacity_coefficient_pure(const double T, const double P);
		std::vector<double> fugacity_coefficient(const double T, const double P, const std::vector<double>& x);

	private:

		const std::vector<double>& Tc_;
		const std::vector<double>& Pc_;
		const std::vector<double>& omega_;
		const std::vector<double>& MW_;

		pengrobinson_mix* pr;
};


#include "fugacityGPengRobinson_mix.hpp"

#endif	// OPENSMOKEPP_FUGACITYG_PENGROBINSON_MIX_H
