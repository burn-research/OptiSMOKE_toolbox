#ifndef GRAMMAR_BATCHREACTOR_H
#define GRAMMAR_BATCHREACTOR_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OptiSMOKE
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
																"@PressureCoefficient @PressureProfile @TemperatureProfile") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PressureProfile", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary defining the pressure profile",
																false ,
																"none",
																"none",
																"@PressureCoefficient @VolumeProfile") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureProfile", 
																OpenSMOKE::SINGLE_DICTIONARY, 
																"Dictionary defining the temperature profile",
																false ,
																"none",
																"none",
																"@VolumeProfile") );

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

#endif // GRAMMAR_BATCHREACTOR_H
