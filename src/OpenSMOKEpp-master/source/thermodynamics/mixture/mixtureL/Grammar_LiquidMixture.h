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

#ifndef OPENSMOKEPP_GRAMMAR_LIQUIDMIXTURE_H
#define OPENSMOKEPP_GRAMMAR_LIQUIDMIXTURE_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "math/OpenSMOKEFunctions.h"

class Grammar_LiquidMixture : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord	(	OpenSMOKE::OpenSMOKE_DictionaryKeyWord
						(	"@Density", OpenSMOKE::SINGLE_STRING,
							"Density rule: Peng-Robinson", false	));

			AddKeyWord	(	OpenSMOKE::OpenSMOKE_DictionaryKeyWord
						(	"@HeatCapacity", OpenSMOKE::SINGLE_STRING,
							"Heat Capacity rule: Ideal", false	));

			AddKeyWord	(	OpenSMOKE::OpenSMOKE_DictionaryKeyWord
						(	"@Viscosity", OpenSMOKE::SINGLE_STRING,
							"Viscosity rule: GrunbergNissanIdeal", false	));

			AddKeyWord	(	OpenSMOKE::OpenSMOKE_DictionaryKeyWord
						(	"@Conductivity", OpenSMOKE::SINGLE_STRING,
							"Conductivity rule: Vredeveld", false	));

			AddKeyWord	(	OpenSMOKE::OpenSMOKE_DictionaryKeyWord
						(	"@Dinf", OpenSMOKE::SINGLE_STRING,
							"Correlation for infinite dilution diffusion coefficients: Siddiqi-Lucas", false	));

			AddKeyWord	(	OpenSMOKE::OpenSMOKE_DictionaryKeyWord
						(	"@DiffusionCorrection", OpenSMOKE::SINGLE_DOUBLE,
							"Correction for infinite dilution diffusion coefficients", false	));

			AddKeyWord	(	OpenSMOKE::OpenSMOKE_DictionaryKeyWord
						(	"@ConductivityCorrection", OpenSMOKE::SINGLE_DOUBLE,
							"Correction for liquid conductivity", false	));

			AddKeyWord	(	OpenSMOKE::OpenSMOKE_DictionaryKeyWord
						(	"@Activity", OpenSMOKE::SINGLE_STRING,
							"Method for activity calculation: Ideal || Unifac", false	));

			AddKeyWord	(	OpenSMOKE::OpenSMOKE_DictionaryKeyWord
						(	"@Fugacity", OpenSMOKE::SINGLE_STRING,
							"Fugacity rule: Raoult | PengRobinsonIndirect | PengRobinsonDirect", false	));
		}
};

#endif	// OPENSMOKEPP_GRAMMAR_LIQUIDMIXTURE_H
