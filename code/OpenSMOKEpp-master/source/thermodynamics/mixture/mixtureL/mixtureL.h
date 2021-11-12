/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alessandro Stagni <alessandro.stagni@polimi.it>               |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2016, 2015 Alessandro Stagni & Alberto Cuoci             |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

/**
 \class mixtureL
 \brief Class to manage liquid-phase mixture properties
 
 Class to manage liquid-phase mixture properties
 */

#ifndef OPENSMOKEPP_MIXTURE_L_H
#define OPENSMOKEPP_MIXTURE_L_H

#include "thermodynamics/maps/speciesMap.h"
#include "Grammar_LiquidMixture.h"
#include "thermodynamics/mixture/eos/eos_mix.h"
#include "thermodynamics/mixture/density/densityL_mix.h"
#include "thermodynamics/mixture/heatcapacity/heatcapacity_mix.h"
#include "thermodynamics/mixture/viscosity/viscosity_mix.h"
#include "thermodynamics/mixture/conductivity/conductivity_mix.h"
#include "thermodynamics/mixture/diffusivity/diffusivityL.h"
#include "thermodynamics/mixture/activity/activity.h"
#include "thermodynamics/mixture/fugacity/fugacity_mix.h"

class mixtureL
{
	public:

		mixtureL(const std::vector<std::string>& species, speciesMap& map);
		mixtureL(const std::vector<std::string>& species, speciesMap& map, OpenSMOKE::OpenSMOKE_Dictionary& dictionary);
		~mixtureL();

		double rhoL_mix(const double T, const double P, const std::vector<double>& x) const;
		
		double d_rhoL_over_dT_mix(const double T, const double P, const std::vector<double>& x) const;
		
		double cpL_mix(const double T, const std::vector<double>& x) const;
		
		double etaL_mix(const double T, const std::vector<double>& x) const;
		
		double lambdaL_mix(const double T, const std::vector<double>& x) const;
		
		const std::vector<double> gamma(const double T, const double P, const std::vector<double>& x) const;
		
		const std::vector<double> fugacity(const double T, const double P, const std::vector<double>& x) const;
		
		const Eigen::MatrixXd Dinf(const double T, const double P) const;
		
		const Eigen::MatrixXd D_sm(const double T, const double P, const std::vector<double>& x) const;
		
		const Eigen::MatrixXd B(const double T, const double P, const std::vector<double>& x) const;
		
		const Eigen::MatrixXd thermodynamicmatrix_sm(const double T, const double P, const std::vector<double>& x) const;
		
		const Eigen::VectorXd molar_diffusion_fluxes(	const double T, const double P,
														const std::vector<double>& x, const double ctot, 
														const std::vector<double>& dx) const;
  
		void SetStefanMaxwellLastSpecies(const std::string species);

	private:

		void IndependentProperties();
		
		void SetupProperties();
		
		void MemoryAllocation();
		
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

	private:

		const std::vector<std::string>& species_;
		speciesMap& map_;

		unsigned int NS_;
		std::vector<int> mapindex_;
  
		//Independent properties
		std::vector<double> Tc_;
		std::vector<double> Pc_;
		std::vector<double> MW_;
		std::vector<double> Rhoc_;
		std::vector<double> Tnbp_;
		std::vector<double> dipolMoment_;
		std::vector<double> omega_;
		std::vector<double> specific_volume_stp_;
  
		//Unifac properties
		std::vector< std::vector<double> > unifacgroup_R_;
		std::vector< std::vector<double> > unifacgroup_Q_;
		std::vector< std::vector<int> > unifacgroup_subgroup_;
		std::vector< std::vector<int> > unifacgroup_maingroup_;
		std::vector< std::vector<int> > unifacgroup_number_;
  
		const vector<double> rhoL(double T, double P) const;

		//Options rule
		unsigned int density_rule_;
		unsigned int capacity_rule_;
		unsigned int viscosity_rule_;
		unsigned int conductivity_rule_;
		unsigned int dinf_rule_;
		unsigned int activity_rule_;
		unsigned int fugacity_rule_;
  
		//Single-species model objects
		std::vector<DensityModel*> dM_;
		std::vector<VaporPressureModel*> vpM_;
  
		//Mixture model objects
		densityLModel_mix* rhoLM_;
		CapacityModel_mix* cpLM_;
		ViscosityModel_mix* etaLM_;
		ConductivityModel_mix* lambdaLM_;
		DinfModel* DinfM_;
		stefanmaxwell* smM_;
		ActivityModel* activityM_;
		FugacityModel_mix* fugacityM_;
  
		double diffusion_correction_;
		double conductivity_correction_;
};

#include "mixtureL.hpp"

#endif	// OPENSMOKEPP_MIXTURE_L_H
