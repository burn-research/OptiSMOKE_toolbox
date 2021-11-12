/*-----------------------------------------------------------------------*\
|     ____            _  ______ __  __  ____  _  ________                 |
|    / __ \       _  (_)/  ___ |  \/  |/ __ \| |/ /  ____|                |
|   | |  | |_ __ | |_ _|  (___ | \  / | |  | | ' /| |__    _     _        |
|   | |  | | '_ \|  _| |\___  \| |\/| | |  | |  < |  __| _| |_ _| |_      |
|   | |__| | |_) | |_| |____)  | |  | | |__| | . \| |___|_   _|_   _|     |
|    \____/| .__/\___|_|______/|_|  |_|\____/|_|\_\______||_|   |_|       |
|          | |                                                            |
|          |_|                                                            |
|                                                                         |
|            Authors: Magnus Fürst <magnus.furst@ulb.ac.be>               |
|                     Andrea Bertolino <andrea.bertolino@ulb.be>          |
|-------------------------------------------------------------------------|
|   License                                                               |
|                                                                         |
|   This file is part of OptiSMOKE.                                       |
|   Copyright (C) 2020 by Magnus Fürst and Andrea Bertolino               |
|                                                                         |
|   OptiSMOKE is free software: you can redistribute it and/or modify     |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OptiSMOKE is distributed in the hope that it will be useful,          |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OptiSMOKE. If not, see <http://www.gnu.org/licenses/>.     |
|                                                                         |
\*-----------------------------------------------------------------------*/


#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	class Grammar_BatchReactor : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsFolder", 
																OpenSMOKE::SINGLE_PATH, 
																"Name of the folder containing the kinetic scheme (XML Version)", 
																true,
																"@KineticsPreProcessor",
																"none",
																"none") );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Name of the dictionary containing the list of kinetic files to be interpreted", 
																true,
																"@KineticsFolder",
																"none",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type", 
																OpenSMOKE::SINGLE_STRING, 
																"Batch reactor type: Isothermal-ConstantVolume | NonIsothermal-ConstantVolume | Isothermal-ConstantPressure | NonIsothermal-ConstantPressure | NonIsothermal-UserDefinedVolume", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InitialStatus", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Name of the dictionary defining the initial gas composition, temperature and pressure", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EndTime", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Ending time for transient simulation (i.e. 0.1 s)", 
																true) );	

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@StartTime", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Start time for transient simulation (default 0)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Volume", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Volume of the reactor (eg 30 cm3)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SensitivityAnalysis", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary containing additional options for solving the sensitivity analysis", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Options", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary containing additional options for solving the batch reactor", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@GlobalThermalExchangeCoefficient", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Global thermal exchange coefficient U: Q = UA(T-Tenv)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EnvironmentTemperature", 
																OpenSMOKE::SINGLE_MEASURE, 
																"EnvironmentTemperature Tenv: Q = UA(T-Tenv)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ExchangeArea", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Exchange area A: Q = UA(T-Tenv)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OdeParameters", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary containing the numerical parameters for solving the stiff ODE system", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@VolumeProfile", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary defining the volume profile",
																false ,
																"none",
																"none",
																"@PressureCoefficient @PressureProfile") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PressureProfile", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary defining the pressure profile",
																false ,
																"none",
																"none",
																"@PressureCoefficient @VolumeProfile") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PressureCoefficient", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Coefficient for pressure increase (eg 1e5 Pa/s)",
																false ,
																"none",
																"none",
																"@VolumeProfile @PressureProfile") );

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ParametricAnalysis",
															   OpenSMOKE::SINGLE_DICTIONARY,
															   "Dictionary containing additional options for performing a parametric analysis",
															   false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OnTheFlyROPA",
																OpenSMOKE::SINGLE_DICTIONARY,
																"Dictionary specifying the details for carrying out the ROPA (on the fly)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OnTheFlyCEMA",
																OpenSMOKE::SINGLE_DICTIONARY,
																"Dictionary specifying the details for carrying out the CEMA (on the fly)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OnTheFlyPostProcessing",
																OpenSMOKE::SINGLE_DICTIONARY,
																"Dictionary specifying the details for carrying out the post-processing analyses (on the fly)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PolimiSoot",
																OpenSMOKE::SINGLE_DICTIONARY,
																"Name of the dictionary defining the rules for analyzing soot calculated using the Polimi mechanism",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@IgnitionDelayTimes",
															   OpenSMOKE::SINGLE_DICTIONARY,
															   "Dictionary containing additional options for estimating the ignition delay times",
															   false));
		}
	};
}

