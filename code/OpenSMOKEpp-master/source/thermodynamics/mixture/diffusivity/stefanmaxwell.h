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
 \class stefanmaxwell
 \brief Class managing liquid diffusion via Stefan-Maxwell liquid multicomponent model
 
 Class managing liquid diffusion via Stefan-Maxwell liquid multicomponent model
 */

#ifndef OPENSMOKEPP_STEFANMAXWELL_H
#define OPENSMOKEPP_STEFANMAXWELL_H

#include "DinfModel.h"
#include "thermodynamics/mixture/activity/activity.h"

class stefanmaxwell
{
	public:

		stefanmaxwell(DinfModel& dinf_model, ActivityModel& activity_model);
  
		void UpdateBinaryDiffusivities(const std::vector<double>& x);
		void UpdateB(const double T, const double P, const std::vector<double>& x);
		void UpdateThermodynamicMatrix(const double T, const double P, const std::vector<double>& x);
  
		void Update(const double T, const double P, const std::vector<double>& x);
  
		void SetLastIndex(const unsigned int index) { last_index_ = index; }
  
		inline const Eigen::MatrixXd& D_sm() const { return D_sm_; }
		inline const Eigen::MatrixXd& B() const { return B_; }
		inline const Eigen::MatrixXd& tm() const { return tm_; }
  
		Eigen::VectorXd molar_diffusion_fluxes(const double ctot, const std::vector<double>& dx);
		Eigen::VectorXd molar_diffusion_fluxes(const double T, const double P, const std::vector<double>& x, const double ctot, const std::vector<double>& dx);
  
	private:  

		unsigned int NS_;
		Eigen::MatrixXd D_sm_;
		Eigen::MatrixXd B_;
		Eigen::MatrixXd tm_;
  
		//Fluxes evaluation
		Eigen::MatrixXd A_;
		Eigen::VectorXd b_;
  
		Eigen::MatrixXd* Dinf_;
		DinfModel& dinf_model_;
		ActivityModel& activity_model_;
  
		unsigned int last_index_;
};

#include "stefanmaxwell.hpp"

#endif	// OPENSMOKEPP_STEFANMAXWELL_H
