#ifndef OPTISMOKE_PLUGFLOWREACTOR_H
#define OPTISMOKE_PLUGFLOWREACTOR_H

// Utilities
#include "idealreactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"
#include "utilities/cema/OnTheFlyCEMA.h"
#include "utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h"
#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"

// Standard plug flow reactors
#include "grammar/grammar_plugflowreactor.h"
#include "idealreactors/plugflow/PlugFlowReactor"

namespace OptiSMOKE
{
	class PlugFlowReactor
	{
	public:

		void Setup(	const std::string input_file_name,
					OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
					OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML);

		void Solve();

		double GetMolefraction(std::string species_name);

	private:

		OpenSMOKE::PlugFlowReactor_Type type_;
		OpenSMOKE::PlugFlowReactor_Isothermal* plugflow_isothermal_;
		OpenSMOKE::PlugFlowReactor_NonIsothermal* plugflow_non_isothermal_;


		// Read thermodynamics and kinetics maps
		OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;

		OpenSMOKE::PolimiSoot_Analyzer* polimi_soot_;
		OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing_;
		OpenSMOKE::OnTheFlyROPA* onTheFlyROPA_;
		OpenSMOKE::PlugFlowReactor_Options*	plugflow_options_;
		OpenSMOKE::ODE_Parameters* ode_parameters_;
		OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options_;
		OpenSMOKE::IgnitionDelayTimes_Analyzer* idt_;
		OpenSMOKE::PlugFlowReactor_Profile* profile_;

		double end_value_;
		bool temperature_profile;

		void CleanMemory();
	};
}

#include "PlugFlowReactor.hpp"

#endif // OPTISMOKE_PLUGFLOWREACTOR_H