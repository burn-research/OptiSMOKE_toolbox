/*----------------------------------------------------------------------*\
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
|   This file is part of OpenSMOKE++ Suite.                               |
|                                                                         |
|   Copyright(C) 2015, 2014, 2013  Alberto Cuoci                          |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_Grammar_LaminarFlame_H
#define OpenSMOKE_Grammar_LaminarFlame_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	class Grammar_LaminarFlame : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsFolder",
				OpenSMOKE::SINGLE_PATH,
				"Name of the folder containing the kinetic scheme (XML Version)",
				true,
				"@KineticsPreProcessor",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Backup",
				OpenSMOKE::SINGLE_PATH,
				"Name of backup file (XML Version)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DontUseBackupGrid",
				OpenSMOKE::SINGLE_BOOL,
				"If true, the user defined grid is used, instead of the grid corresponding to the backup file (default:false)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary containing the list of kinetic files to be interpreted",
				true,
				"@KineticsFolder",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type",
				OpenSMOKE::SINGLE_STRING,
				"Type of simulation: BurnerStabilized | FlameSpeed",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Grid",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary defining the mesh",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletStream",
				OpenSMOKE::VECTOR_STRING,
				"Name of the dictionary/dictionaries defining the inlet gas composition, temperature and pressure",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OutletStream",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary defining the outlet gas composition, temperature and pressure (this is only for initialization)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletVelocity",
				OpenSMOKE::SINGLE_MEASURE,
				"Inlet velocity (first guess in case of flame speed calculations)",
				true,
				"@InletMassFlux",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletMassFlux",
				OpenSMOKE::SINGLE_MEASURE,
				"Inlet mass flux (first guess in case of flame speed calculations)",
				true,
				"@InletVelocity",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SensitivityAnalysis",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing additional options for solving the sensitivity analysis",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DaeParameters",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing the numerical parameters for solving the stiff DAE system",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NlsParameters",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing the numerical parameters for solving the NL system",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FalseTransientParameters",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing the numerical parameters for solving the pseudo-transient phase",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Output",
				OpenSMOKE::SINGLE_PATH,
				"Name of the folder where to write the output files",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@UseNlsSolver",
				OpenSMOKE::SINGLE_BOOL,
				"Use the NLS solver to solve the steady state problems, after the DAE solver (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@UseDaeSolver",
				OpenSMOKE::SINGLE_BOOL,
				"Use the Dae solver (instead of NLS) to solve the steady state problems (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PolimiSoot",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary defining the rules for analyzing soot calculated using the Polimi mechanism",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@HMOM",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary defining the rules for applying the Hybrid Method of Moments (HMOM)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OnTheFlyPostProcessing",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary specifying the details for carrying out the post-processing analyses (on the fly)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FixedTemperatureProfile",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of dictionary describing the temperature profile",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FixedSpecificMassFlowRateProfile",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of dictionary describing the specific (i.e. per unit area) mass flow rate profile",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Soret",
				OpenSMOKE::SINGLE_BOOL,
				"Add Soret effect (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SpeciesBundling",
				OpenSMOKE::SINGLE_DOUBLE,
				"Estimation of mass diffusion coefficients through the species bundling according the the specified maximum relative error (example: @SpeciesBundling 0.1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Radiation",
				OpenSMOKE::SINGLE_BOOL,
				"Radiative heat transfer (default: none)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EnvironmentTemperature",
				OpenSMOKE::SINGLE_MEASURE,
				"Environment temperature (default: 298.15 K)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DerivativeGasMassFractions",
				OpenSMOKE::SINGLE_STRING,
				"Derivative gas mass fractions (default: backward)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DerivativeGasTemperature",
				OpenSMOKE::SINGLE_STRING,
				"Derivative gas temperature (default: backward)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@LewisNumbers",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of dictionary containing the list of Lewis numbers of species",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FixedOutletTemperature",
				OpenSMOKE::SINGLE_MEASURE,
				"Fixed outlet temperature (to be used only if @FixedSpecificMassFlowRateProfile is on)",
				false,
				"none",
				"@FixedSpecificMassFlowRateProfile",
				"none"));


			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@WallHeatExchangeCoefficient",
				OpenSMOKE::SINGLE_MEASURE,
				"Wall heat exchange coefficient (default: 0 W/m2/K)",
				false,
				"none",
				"none",
				"@WallHeatNusseltNumber"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@WallHeatNusseltNumber",
				OpenSMOKE::SINGLE_DOUBLE,
				"Wall heat Nusselt number (default: 0)",
				false,
				"none",
				"none",
				"@WallHeatExchangeCoefficient"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InternalDiameter",
				OpenSMOKE::SINGLE_MEASURE,
				"Internal diameter (default: 1 cm)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@WallTemperatureProfile",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of dictionary describing the wall temperature profile",
				false));
		}
	};
}

#endif /* OpenSMOKE_Grammar_LaminarFlame_H */