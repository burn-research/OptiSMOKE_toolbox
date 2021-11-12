/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
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
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
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

#ifndef OpenSMOKE_BatchReactor_NonIsothermal_ConstantVolume_H
#define	OpenSMOKE_BatchReactor_NonIsothermal_ConstantVolume_H

// Parent class
#include "BatchReactor.h"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"

// Options
#include "BatchReactor_Options.h"
#include "math/external-ode-solvers/ODE_Parameters.h"


namespace OpenSMOKE
{
	//!  A class for simulating batch reactors with constant volume in adiabatic conditions
	/*!
		 The purpose of this class is to simulate a batch reactor with constant volume, in adiabatic conditions
		 The conservation equations of species and energy are solved in terms of mass fractions and
		 temperature respectively. The internal energy (massive) is conserved
	*/

	class BatchReactor_NonIsothermal_ConstantVolume : public virtual BatchReactor
	{

	public:


		std::vector<double> time_vector;
		std::vector<OpenSMOKE::OpenSMOKEVectorDouble> species_matrix;

		/**
		*@brief Default constructor
		*@param thermodynamicsMap map containing the thermodynamic data
		*@param kineticsMap map containing the kinetic mechanism
		*@param ode_parameters parameters governing the solution of the stiff ODE system
		*@param batch_options options governing the output
		*@param on_the_fly_ropa rate of production analysis (on the fly)
		*@param on_the_fly_cema chemical explosive mode analysis (on the fly)
		*@param on_the_fly_post_processing post-processing analysis (on the fly)
		*@param idt_analyzer ignition delay time analyzer (on the fly)
		*@param polimi_soot_analyzer soot analyzer based on the POLIMI model (on the fly)
		*@param V0 initial volume [m3]
		*@param T0 initial temperature [K]
		*@param P0 initial pressure [Pa]
		*@param omega0 initial mass fractions of species
		*@param global_thermal_exchange_coefficient global thermal exchange coefficient [W/m2/K]
        *@param exchange_area thermal exchange area available [m2]
        *@param T_environment environment temperature [K]
		*/
		BatchReactor_NonIsothermal_ConstantVolume(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
													OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
													OpenSMOKE::ODE_Parameters& ode_parameters,
													OpenSMOKE::BatchReactor_Options& batch_options,
													OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa,
													OpenSMOKE::OnTheFlyCEMA& on_the_fly_cema,
													OpenSMOKE::OnTheFlyPostProcessing& on_the_fly_post_processing,
													OpenSMOKE::IgnitionDelayTimes_Analyzer& idts_analyzer,
													OpenSMOKE::PolimiSoot_Analyzer& polimi_soot_analyzer,
													const double V0, const double T0, const double P0, 
													const OpenSMOKE::OpenSMOKEVectorDouble& omega0,
                									const double global_thermal_exchange_coefficient,
                                                    const double exchange_area,
                                                    const double T_environment);

		/**
		*@brief Solves the batch reactor
                *@param t0 the starting time of integration [s]
		*@param tf the final time of integration [s]
		*/
		virtual void Solve(const double t0, const double tf);

		/**
		*@brief Ordinary differential Equations corresponding to the reactor under investigation
		*@param t current time [s]
		*@param y current solution
		*@param dy current unsteady terms
		*/
		virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

		/**
		*@brief Sparse Analystical Jacobian
		*@param t current time [s]
		*@param y current solution
		*@param J sparse analytical Jacobian
		*/
		virtual void SparseAnalyticalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::SparseMatrix<double> &J);

		/**
		*@brief Dense Analytical Jacobian
		*@param t current time [s]
		*@param y current solution
		*@param J sparse analytical Jacobian
		*/
		virtual void DenseAnalyticalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::MatrixXd &J);

		/**
		*@brief Writes the output (called at the end of each time step)
		*@param t current time [s]
		*@param y current solution
		*/
		virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	protected:

		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param t current time [s]
		*@param y current solution
		*/
		void SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);
		void ChemicalExplosiveModeAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	private:

		OpenSMOKE::OpenSMOKEVectorDouble U_;		//!< current internal energy [J/kmol]
		OpenSMOKE::OpenSMOKEVectorDouble u_;		//!< current internal energy [J/kg]

		Eigen::SparseMatrix<double> Jomega_;
	};
}

#include "BatchReactor_NonIsothermal_ConstantVolume.hpp"

#endif	/* OpenSMOKE_BatchReactor_NonIsothermal_ConstantVolume_H */

