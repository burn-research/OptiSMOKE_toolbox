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
	class Grammar_PerfectlyStirredReactor : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
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
																"Perfectly Stirred Reactor type: Isothermal-ConstantPressure | NonIsothermal-ConstantPressure", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InitialStatus", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Name of the dictionary defining the initial gas composition, temperature and pressure", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletStatus", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Name of the dictionary defining the inlet gas composition, temperature and pressure", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ResidenceTime", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Residence time (i.e. 0.1 s)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Volume", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Volume of the reactor (eg 30 cm3)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MassFlowRate", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Inlet mass flow rate (eg 10 g/s)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EndTime", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Ending time for transient simulation (i.e. 0.1 s)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SensitivityAnalysis", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary containing additional options for solving the sensitivity analysis", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Options", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary containing additional options for solving the perfectly stirred reactor", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@GlobalThermalExchangeCoefficient", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Global thermal exchange coefficient U: Q = UA(T-Tenv)", 
																false,
																"none",
																"@ExchangeArea",
																"none") );

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

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ParametricAnalysis",
																OpenSMOKE::SINGLE_DICTIONARY,
																"Dictionary containing additional options for performing a parametric analysis",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OnTheFlyROPA",
																OpenSMOKE::SINGLE_DICTIONARY,
																"Dictionary specifying the details for carrying out the ROPA (on the fly)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OnTheFlyPostProcessing",
																OpenSMOKE::SINGLE_DICTIONARY,
																"Dictionary specifying the details for carrying out the post-processing analyses (on the fly)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PolimiSoot",
																OpenSMOKE::SINGLE_DICTIONARY,
																"Name of the dictionary defining the rules for analyzing soot calculated using the Polimi mechanism",
																false));
		}
	};
}

