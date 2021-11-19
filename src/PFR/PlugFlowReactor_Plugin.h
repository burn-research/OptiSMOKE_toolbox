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

#ifndef OpenSMOKE_PlugFlowReactor_Plugin_H
#define OpenSMOKE_PlugFlowReactor_Plugin_H

namespace OpenSMOKE
{
	class PlugFlowReactor_Plugin
	{
	public:

		void Setup(	const std::string input_file_name,
					OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML,
					OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML);


		double Solve_tau(std::string tau_calc_type_temp);
		void Update_and_Solve_PFR(OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML, OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML);
		double Solve_Species(std::string Species);
		std::vector<std::vector<double>> Solve_Multipl_Species_time_profile(std::vector<std::vector<std::vector<double>>> &Exp_data_temp, std::vector<std::string> Species_vec);
		std::vector<double> Solve_Multipl_Species_outlet(std::vector<std::string> Species_vec);

		double Interpolate(std::vector<double> Time_vec_temp, std::vector<OpenSMOKE::OpenSMOKEVectorDouble> Species_matrix_temp, double Abscissa_temp, std::string Species_name);
		//std::vector<double> Solve_Species(std::vector<std::string> Species);
		void clean_up();

		double time_fifty_percent_conv_exp;
		double time_fifty_percent_conv_sim;
		double time_shift;

	private:

		OpenSMOKE::PlugFlowReactor_Type 		type_;
		OpenSMOKE::PlugFlowReactor_Isothermal* 		plugflow_isothermal_;
		OpenSMOKE::PlugFlowReactor_NonIsothermal* 	plugflow_non_isothermal_;


		// Read thermodynamics and kinetics maps
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*		thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML_;

//		OpenSMOKE::OptimizationRules*			optimization_;
		OpenSMOKE::PolimiSoot_Analyzer*			polimi_soot_;
		OpenSMOKE::OnTheFlyPostProcessing*		on_the_fly_post_processing_;
		OpenSMOKE::OnTheFlyROPA*			onTheFlyROPA_;
		OpenSMOKE::PlugFlowReactor_Options*		plugflow_options_;
		OpenSMOKE::ODE_Parameters*			ode_parameters_;
		OpenSMOKE::SensitivityAnalysis_Options*		sensitivity_options_;
//		OpenSMOKE::VirtualChemistry*			virtual_chemistry_;
		OpenSMOKE::IgnitionDelayTimes_Analyzer*		idt;

		double 						end_value_;
		double 						T;
		double 						velocity;

		bool 						constant_pressure;
		double 						P_Pa;
		OpenSMOKE::OpenSMOKEVectorDouble 		omega;
		bool 						time_independent_variable;
		double 						cross_section_over_perimeter;
		double 						global_thermal_exchange_coefficient;
		double 						T_environment;

		bool 						temperature_profile;
		OpenSMOKE::PlugFlowReactor_Profile*		profile;

		double 						tau_ign_temp;
		int 						pos_OH_max;
		double 						Mole_frac_temp;
		double 						T_Final;
		double 						P_Pa_Final; 
		double 						MW_Final;
		std::vector<double> 				Mole_frac_temp_vec;
	};
}

#include "PlugFlowReactor_Plugin.hpp"

#endif // OpenSMOKE_PlugFlowReactor_Plugin_H
