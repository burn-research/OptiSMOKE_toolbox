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

#ifndef Flame1D_Plugin_H
#define Flame1D_Plugin_H

namespace OpenSMOKE
{
	//OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D* flame_premixed;
	
	class Flame1D_Plugin
	{
	public:

		std::vector<double> Setup_and_Solve(const std::string input_file_name,
			OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML,
			OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML,
			OpenSMOKE::TransportPropertiesMap_CHEMKIN* 		transportMapXML,
			std::string Species_of_interest,
			std::vector<double> Abscissa_x_m);

		double Interpolate(Eigen::VectorXd grid_coordinates, std::vector<double> species_profile, double Abscissa_temp);
		
		void clean_up();

	private:


		// Read thermodynamics and kinetics maps
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML_;
		OpenSMOKE::TransportPropertiesMap_CHEMKIN*	transportMapXML_;
	
		OpenSMOKE::Grid1D* grid;
		OpenSMOKE::PolimiSoot_Analyzer* polimi_soot;
		OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing;
		OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options;
		DaeSMOKE::DaeSolver_Parameters* dae_parameters;
		NlsSMOKE::NonLinearSolver_Parameters* nls_parameters;
		NlsSMOKE::FalseTransientSolver_Parameters* false_transient_parameters;
		OpenSMOKE::HMOM* hmom;
		std::shared_ptr<OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D> flame_premixed;
		bool polimi_soot_used;
		bool sensitivity_options_used;
		bool hmom_used;
		

	};
}

//#include "Flame1D_Plugin.hpp"

#endif // 
