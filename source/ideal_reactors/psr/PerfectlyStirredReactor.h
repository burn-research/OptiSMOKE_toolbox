#ifndef OPTISMOKE_PERFECTLYSTIRREDREACTOR_H
#define OPTISMOKE_PERFECTLYSTIRREDREACTOR_H

// Utilities
#include "idealreactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"
#include "utilities/cema/OnTheFlyCEMA.h"
#include "utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h"
#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"

// Oscillating perfectly stirred reactor
#include "perfectlystirredreactor/oscillations/OscillatingPerfectlyStirredReactor.h"

// Standard perfectly stirred
#include "grammar/grammar_perfectlystirredreactor.h"
#include "idealreactors/psr/PerfectlyStirredReactor"

namespace OptiSMOKE
{
	class PerfectlyStirredReactor
	{
	public:

		void Setup(	const std::string input_file_name,
					OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
					OpenSMOKE::KineticsMap_CHEMKIN*	kineticsMapXML);

		void Solve();

		std::vector<double> GetMolefractionsOut(std::vector<std::string> targets_names);

	private:

		OpenSMOKE::PerfectlyStirredReactor_Type type_;
		OpenSMOKE::PerfectlyStirredReactor_Isothermal_ConstantPressure* psr_isothermal_;
		OpenSMOKE::PerfectlyStirredReactor_NonIsothermal_ConstantPressure* psr_non_isothermal_;

		// Read thermodynamics and kinetics maps
		OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;

		OpenSMOKE::PolimiSoot_Analyzer*	polimi_soot_;
		OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing_;
		OpenSMOKE::OnTheFlyROPA* onTheFlyROPA_;
		OpenSMOKE::PerfectlyStirredReactor_Options*	psr_options_;
		OpenSMOKE::ODE_Parameters* ode_parameters_;
		OpenSMOKE::SensitivityAnalysis_Options*	sensitivity_options_;
		OpenSMOKE::OscillatingPerfectlyStirredReactor* oscillating_psr_;

		double P_Pa_Initial;
		double T_Initial;
		double T_Inlet;
		double P_Pa_Inlet;
		double residence_time;
		double mass_flow_rate;
		double volume;
		double exchange_area;
		double global_thermal_exchange_coefficient;
		double T_environment;
		double tEnd;
		OpenSMOKE::OpenSMOKEVectorDouble omega_Initial;
		OpenSMOKE::OpenSMOKEVectorDouble omega_Inlet;
		void CleanMemory();
	};
} // namespace OptiSMOKE

#include "PerfectlyStirredReactor.hpp"

#endif // OPTISMOKE_PERFECTLYSTIRREDREACTOR_H
