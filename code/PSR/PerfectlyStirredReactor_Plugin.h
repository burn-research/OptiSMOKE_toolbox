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
|            Author: Magnus Fürst <magnus.furst@ulb.ac.be>                |
|-------------------------------------------------------------------------|
|   License                                                               |
|                                                                         |
|   This file is part of OptiSMOKE.                                       |
|   Copyright (C) 2019 by Magnus Fürst                                    |
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

#ifndef OpenSMOKE_PerfectlyStirredReactor_Plugin_H
#define OpenSMOKE_PerfectlyStirredReactor_Plugin_H


namespace OpenSMOKE
{
	class PerfectlyStirredReactor_Plugin
	{
	public:

		void Setup(	const std::string input_file_name,
					OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML,
					OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML);

		void Update_and_Solve_PSR(OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML, OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_updated);
		double Solve_Species(std::string Species);
		std::vector<double> Solve_Species(std::vector<std::string> Species);
		std::vector<double> Solve_Multipl_Species(std::vector<std::string> Species_vec);

		void clean_up();




	private:

		OpenSMOKE::PerfectlyStirredReactor_Type type_;
		OpenSMOKE::PerfectlyStirredReactor_Isothermal_ConstantPressure* psr_isothermal_;
		OpenSMOKE::PerfectlyStirredReactor_NonIsothermal_ConstantPressure* psr_non_isothermal_;


		// Read thermodynamics and kinetics maps
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML_;

		OpenSMOKE::PolimiSoot_Analyzer*			polimi_soot_;
		OpenSMOKE::OnTheFlyPostProcessing*		on_the_fly_post_processing_;
		OpenSMOKE::OnTheFlyROPA*			onTheFlyROPA_;
		OpenSMOKE::PerfectlyStirredReactor_Options*		psr_options_;
		OpenSMOKE::ODE_Parameters*			ode_parameters_;
		OpenSMOKE::SensitivityAnalysis_Options*		sensitivity_options_;
		//OpenSMOKE::VirtualChemistry*			virtual_chemistry_;

		double P_Pa_Initial;
		double T_Initial;
		OpenSMOKE::OpenSMOKEVectorDouble omega_Initial;
		double T_Inlet;
		double P_Pa_Inlet;
		OpenSMOKE::OpenSMOKEVectorDouble omega_Inlet;
		double residence_time;
		double mass_flow_rate;
		double volume;
		double exchange_area;
		double global_thermal_exchange_coefficient;
		double T_environment;
		double tEnd;

		double Mole_frac_temp;
		std::vector<double> Mole_frac_temp_vec;
		double T_Final;
		double P_Pa_Final;
		double MW_Final;
	};
}

#include "PerfectlyStirredReactor_Plugin.hpp"

#endif // OpenSMOKE_PerfectlyStirredReactor_Plugin_H
