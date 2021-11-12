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

#ifndef BatchReactor_Plugin_H
#define BatchReactor_Plugin_H


namespace OpenSMOKE
{
	
	class BatchReactor_Plugin
	{
	public:

		void Setup(const std::string input_file_name,
			OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML,
			OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML);

		void Update_and_Solve_Batch(OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML, OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML);
		// Function for solving tau according to taul_calc_type
		double Solve_tau(std::string tau_calc_type_temp);
		std::vector<std::vector<double>> Solve_Multiple_Species(std::vector<std::vector<std::vector<double>>> &Exp_data_temp, std::vector<std::string> Species_vec);
		double Interpolate(std::vector<double> Time_vec_temp, std::vector<OpenSMOKE::OpenSMOKEVectorDouble> Species_matrix_temp, double Abscissa_temp, std::string Species_name);

		std::vector<double> Solve_Species(std::vector<double> Abscissa, std::string Species);

		void clean_up();
		
	private:

		OpenSMOKE::BatchReactor_Type type_;
		OpenSMOKE::BatchReactor_NonIsothermal_ConstantPressure* batch_nonisothermal_constantp_;
		OpenSMOKE::BatchReactor_Isothermal_ConstantPressure* batch_isothermal_constantp_;
		OpenSMOKE::BatchReactor_NonIsothermal_ConstantVolume* batch_nonisothermal_constantv_;
		OpenSMOKE::BatchReactor_Isothermal_ConstantVolume* batch_isothermal_constantv_;
		OpenSMOKE::BatchReactor_NonIsothermal_UserDefinedVolume* batch_nonisothermal_userdefinedv_;


		// Read thermodynamics and kinetics maps
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML_;

		OpenSMOKE::PolimiSoot_Analyzer*			polimi_soot_;
		OpenSMOKE::OnTheFlyPostProcessing*		on_the_fly_post_processing_;
		OpenSMOKE::OnTheFlyROPA*				onTheFlyROPA_;
		OpenSMOKE::BatchReactor_Options*		batch_options_;
		OpenSMOKE::ODE_Parameters*				ode_parameters_;
		OpenSMOKE::SensitivityAnalysis_Options*	sensitivity_options_;
		OpenSMOKE::IgnitionDelayTimes_Analyzer*	idt;
		OpenSMOKE::OnTheFlyCEMA*				onTheFlyCEMA_;

		bool volume_profile;
		double tStart_;
		double tEnd_;
		double T, P_Pa;
		OpenSMOKE::OpenSMOKEVectorDouble omega;
		double volume;					// default value [1 m3]
		double exchange_area;
		double global_thermal_exchange_coefficient;
		double T_environment;
		OpenSMOKE::BatchReactor_VolumeProfile* batchreactor_volumeprofile;
		double Mole_frac_temp;
		double T_Final;
		double P_Pa_Final; 
		double MW_Final;

	};
}

#include "BatchReactor_Plugin.hpp"

#endif // OpenSMOKE_BatchReactorExperiment_H
