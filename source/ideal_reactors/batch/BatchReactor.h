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

		OpenSMOKE::BatchReactor_Options* batch_options;
		OpenSMOKE::ODE_Parameters* ode_parameters;
		OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options;
		OpenSMOKE::OnTheFlyROPA* onTheFlyROPA;
		OpenSMOKE::OnTheFlyCEMA* onTheFlyCEMA;
		OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing;
		OpenSMOKE::IgnitionDelayTimes_Analyzer* idt;
		OpenSMOKE::BatchReactor_VolumeProfile* batchreactor_volumeprofile;
		OpenSMOKE::PolimiSoot_Analyzer* polimi_soot;

		double tStart_;
		double tEnd_;

		void CleanMemory();

		bool volume_profile_;

		bool temperature_profile_;
	};
} // namespace OptiSMOKE

#include "BatchReactor.hpp"
#endif // OPTISMOKE_BATCHREACTOR_H
