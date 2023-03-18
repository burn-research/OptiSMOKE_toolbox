#ifndef OPTISMOKE_BATCHREACTOR_H
#define OPTISMOKE_BATCHREACTOR_H

// Utilities
#include "idealreactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"
#include "utilities/cema/OnTheFlyCEMA.h"
#include "utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h"
#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"

// Standard batch reactors
#include "idealreactors/batch/BatchReactor"
#include "grammar/grammar_batchreactor.h"

namespace OptiSMOKE
{
	class BatchReactor
	{
	public:
		void Setup(const std::string input_file_name,
			OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
			OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML);

		void Solve();

		double GetIgnitionDelayTime(std::string criterion);

	private:
		OpenSMOKE::BatchReactor_Type type_;
		OpenSMOKE::BatchReactor_NonIsothermal_ConstantPressure* batch_nonisothermal_constantp_;
		OpenSMOKE::BatchReactor_Isothermal_ConstantPressure* batch_isothermal_constantp_;
		OpenSMOKE::BatchReactor_NonIsothermal_ConstantVolume* batch_nonisothermal_constantv_;
		OpenSMOKE::BatchReactor_Isothermal_ConstantVolume* batch_isothermal_constantv_;
		OpenSMOKE::BatchReactor_NonIsothermal_UserDefinedVolume* batch_nonisothermal_userdefinedv_;


		// Read thermodynamics and kinetics maps
		OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;

		OpenSMOKE::PolimiSoot_Analyzer* polimi_soot_;
		OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing_;
		OpenSMOKE::OnTheFlyROPA* onTheFlyROPA_;
		OpenSMOKE::BatchReactor_Options* batch_options_;
		OpenSMOKE::ODE_Parameters* ode_parameters_;
		OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options_;
		OpenSMOKE::IgnitionDelayTimes_Analyzer* idt;
		OpenSMOKE::OnTheFlyCEMA* onTheFlyCEMA;

		double tStart_;
		double tEnd_;
	};
} // namespace OptiSMOKE

#include "BatchReactor.hpp"
#endif // OPTISMOKE_BATCHREACTOR_H