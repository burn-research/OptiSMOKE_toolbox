#ifndef GRAMMAR_PLUGFLOWREACTOR_H
#define GRAMMAR_PLUGFLOWREACTOR_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OptiSMOKE
{
	class Grammar_PlugFlowReactor : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
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

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary containing the list of kinetic files to be interpreted",
				true,
				"@KineticsFolder",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type",
				OpenSMOKE::SINGLE_STRING,
				"Plug Flow reactor type: Isothermal | NonIsothermal",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletStatus",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Name of the dictionary defining the inlet gas composition, temperature and pressure",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ConstantPressure",
				OpenSMOKE::SINGLE_BOOL,
				"Constant pressure simulation vs Constant velocity simulation",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ResidenceTime",
				OpenSMOKE::SINGLE_MEASURE,
				"Total residence time to simulate(i.e. 0.1 s)",
				true,
				"@Length",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Length",
				OpenSMOKE::SINGLE_MEASURE,
				"Total length to simulate(i.e. 0.1 s)",
				true,
				"@ResidenceTime",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Velocity",
				OpenSMOKE::SINGLE_MEASURE,
				"Inlet velocity of the mixture (e.g. 1 m/s)",
				true,
				"@MassFlowRate @MoleFlowRate @VolumetricFlowRate",
				"none",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MassFlowRate",
				OpenSMOKE::SINGLE_MEASURE,
				"Inlet mass flow rate of the mixture (e.g. 1 kg/s)",
				true,
				"@Velocity @MoleFlowRate @VolumetricFlowRate",
				"@Diameter",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MoleFlowRate",
				OpenSMOKE::SINGLE_MEASURE,
				"Inlet mole flow rate of the mixture (e.g. 1 kmol/s)",
				true,
				"@Velocity @MassFlowRate @VolumetricFlowRate",
				"@Diameter",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@VolumetricFlowRate",
				OpenSMOKE::SINGLE_MEASURE,
				"Inlet volumetric flow rate of the mixture (e.g. 1 cm3/s)",
				true,
				"@Velocity @MassFlowRate @MoleFlowRate",
				"@Diameter",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Diameter",
				OpenSMOKE::SINGLE_MEASURE,
				"Diameter (e.g. 3 cm)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureProfile",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary defining the temperature profile",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SensitivityAnalysis",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing additional options for solving the sensitivity analysis",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Options",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing additional options for solving the plug flow reactor",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@GlobalThermalExchangeCoefficient",
				OpenSMOKE::SINGLE_MEASURE,
				"Global thermal exchange coefficient U: Q = UA(T-Tenv)",
				false,
				"none",
				"@CrossSectionOverPerimeter",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EnvironmentTemperature",
				OpenSMOKE::SINGLE_MEASURE,
				"EnvironmentTemperature Tenv: Q = UA(T-Tenv)",
				false,
				"none",
				"@CrossSectionOverPerimeter",
				"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CrossSectionOverPerimeter",
				OpenSMOKE::SINGLE_MEASURE,
				"Ratio between the cross section and the perimeter (for a circular section it is equal to D/4)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OdeParameters",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing the numerical parameters for solving the stiff ODE system",
				false));

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
				true));
		}
	};
} // namespace OptiSMOKE

#endif // GRAMMAR_PLUGFLOWREACTOR_H