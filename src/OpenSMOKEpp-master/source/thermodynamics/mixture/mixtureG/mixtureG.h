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
 \class mixtureG
 \brief Class to manage gas-phase mixture properties
 
 Class to manage gas-phase mixture properties
 */

#ifndef OPENSMOKEPP_MIXTUREG_H
#define OPENSMOKEPP_MIXTUREG_H

#include "Grammar_GasMixture.h"
#include "thermodynamics/maps/speciesMap.h"
#include "thermodynamics/mixture/eos/eos_mix.h"
#include "thermodynamics/mixture/fugacity/fugacity_mix.h"

class mixtureG
{
	public:

		mixtureG(const std::vector<string>& species, speciesMap& map);
		mixtureG(const std::vector<string>& species, speciesMap& map, OpenSMOKE::OpenSMOKE_Dictionary& dictionary);
		~mixtureG();

		const std::vector<double> fugacity(const double T, const double P, const std::vector<double>& x) const;

	private:

		void CheckBasicSpecies();
		void IndependentProperties();
		void SetupProperties();
		void MemoryAllocation();
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

	private:

		const std::vector<string>& species_;
		speciesMap& map_;

		unsigned int NS_;
		std::vector<int> mapindex_;
  
		//Independent properties
		std::vector<double> Tc_;
		std::vector<double> Pc_;
		std::vector<double> MW_;
		std::vector<double> Rhoc_;
		std::vector<double> Tnbp_;
		std::vector<double> dipolMoment_;
		std::vector<double> omega_;
		std::vector<double> specific_volume_stp_;

		//Options rule
		unsigned int fugacity_rule_;

		//Model objects
		FugacityModel_mix* fugacityM_;
};

#include "mixtureG.hpp"

#endif	// OPENSMOKEPP_MIXTUREG_H
