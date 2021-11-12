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
|	License                                                               |
|                                                                         |
|   Copyright(C) 2018 Alberto Cuoci and Mattia Bissoli                    |
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

#ifndef OpenSMOKE_Grammar_IgnitionDelayTimes_H
#define	OpenSMOKE_Grammar_IgnitionDelayTimes_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{

	class Grammar_IgnitionDelayTimes : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Temperature", 
																OpenSMOKE::SINGLE_BOOL, 
																"Ignition delay time defined on the basis of maximum value and maximum slope of temperature (default: true)", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Pressure", 
																OpenSMOKE::SINGLE_BOOL, 
																"Ignition delay time defined on the basis of maximum value and maximum slope of pressure (default: true)", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Species", 
																OpenSMOKE::VECTOR_STRING, 
																"Ignition delay time defined on the basis of maximum value and maximum slope of selected species (default: none)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SpeciesThreshold", 
																OpenSMOKE::SINGLE_DOUBLE, 
																"Minimum amount (mole fraction) of a given species for evaluating the ignition delay time  (default: 1e-12)", 
																false) );

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MinimumTime",
																OpenSMOKE::SINGLE_MEASURE,
																"Minimum time for evaluating the ignition delay time (default: 1e-12 s)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MinimumTimeInterval",
																OpenSMOKE::SINGLE_MEASURE,
																"Minimum time interval for evaluating the ignition delay time based on the slope values of T, P or species (default: 1e-12 s)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SpeciesSlope",
																OpenSMOKE::SINGLE_BOOL,
																"If true, the maximum slope criterion is also adopted for species (default: true)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureIncrease",
																OpenSMOKE::SINGLE_MEASURE,
																"Ignition delay time defined on the temperature increase with respect to the initial value (default: 0, i.e. false)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PressureIncrease",
																OpenSMOKE::SINGLE_MEASURE,
																"Ignition delay time defined on the pressure increase with respect to the initial value (default: 0, i.e. false)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RapidCompressionMachine",
																OpenSMOKE::SINGLE_BOOL,
																"If true, the analysis is carried out in a post processing phase to select the local ignition delay times (default: false)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FilterWidth",
																OpenSMOKE::SINGLE_MEASURE,
																"Filter width for post processing derivative profile in case of @RapidCompressionMachine (default: 1 ms)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureDerivativeThreshold",
																OpenSMOKE::SINGLE_MEASURE,
																"Threshold for temperature derivative in case of @RapidCompressionMachine (default: 1e3 K/s)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RegularizationTimeInterval",
																OpenSMOKE::SINGLE_MEASURE,
																"Regularization time interval in case of @RapidCompressionMachine (default: 1 ms)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Verbose",
																OpenSMOKE::SINGLE_BOOL,
																"If true, additional data are written in the output files (default: false)",
																false));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TargetMoleFractions", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"The ignition delay is the time at which the mole fraction of target species reaches a specified value (i.e. OH 1e-5 HO2 2e-5)", 
																false));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TargetConcentrations", 
																OpenSMOKE::VECTOR_STRING, 
																"The ignition delay is the time at which the concentration of target species reaches a specified value (i.e. OH 1e-5 HO2 2e-5 mol/cm3)", 
																false));


			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TargetRelativeMoleFractions",
																OpenSMOKE::VECTOR_STRING_DOUBLE,
																"The ignition delay is the time at which the mole fraction of target species reaches a specified value relative to the maximum value (i.e. OH 0.90 HO2 0.80)",
																false));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SpeciesMaxIntercept", 
																OpenSMOKE::VECTOR_STRING, 
																"Extrapolation to the initial baseline mole fraction, from the maximum slope.", 
																false));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SpeciesMinIntercept", 
																OpenSMOKE::VECTOR_STRING, 
																"Extrapolation to the initial baseline mole fraction, from the minimum slope.", 
																false));
		}
	};
}

#endif	/* OpenSMOKE_Grammar_IgnitionDelayTimes_H */

