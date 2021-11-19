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
	class Grammar_PlugFlowReactor : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
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
																"Plug Flow reactor type: Isothermal | NonIsothermal", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletStatus", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Name of the dictionary defining the inlet gas composition, temperature and pressure", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ConstantPressure", 
																OpenSMOKE::SINGLE_BOOL, 
																"Constant pressure simulation vs Constant velocity simulation", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ResidenceTime", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Total residence time to simulate(i.e. 0.1 s)", 
																true,
																"@Length",
																"none",
																"none") );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Length", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Total length to simulate(i.e. 0.1 s)", 
																true,
																"@ResidenceTime",
																"none",
																"none") );
			
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Velocity", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Inlet velocity of the mixture (e.g. 1 m/s)", 
																true,
																"@MassFlowRate @MoleFlowRate @VolumetricFlowRate",
																"none",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MassFlowRate", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Inlet mass flow rate of the mixture (e.g. 1 kg/s)", 
																true,
																"@Velocity @MoleFlowRate @VolumetricFlowRate",
																"@Diameter",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MoleFlowRate", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Inlet mole flow rate of the mixture (e.g. 1 kmol/s)", 
																true,
																"@Velocity @MassFlowRate @VolumetricFlowRate",
																"@Diameter",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@VolumetricFlowRate", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Inlet volumetric flow rate of the mixture (e.g. 1 cm3/s)", 
																true,
																"@Velocity @MassFlowRate @MoleFlowRate",
																"@Diameter",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Diameter", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Diameter (e.g. 3 cm)", 
																false ) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureProfile", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary defining the temperature profile", 
																false ) );
			
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SensitivityAnalysis", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary containing additional options for solving the sensitivity analysis", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Options", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary containing additional options for solving the plug flow reactor", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@GlobalThermalExchangeCoefficient", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Global thermal exchange coefficient U: Q = UA(T-Tenv)", 
																false,
																"none",
																"@CrossSectionOverPerimeter",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EnvironmentTemperature", 
																OpenSMOKE::SINGLE_MEASURE, 
																"EnvironmentTemperature Tenv: Q = UA(T-Tenv)", 
																false,
																"none",
																"@CrossSectionOverPerimeter",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CrossSectionOverPerimeter", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Ratio between the cross section and the perimeter (for a circular section it is equal to D/4)", 
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

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@IgnitionDelayTimes",
															   OpenSMOKE::SINGLE_DICTIONARY,
															   "Dictionary containing additional options for estimating the ignition delay times",
															   false));
		}
	};
}

