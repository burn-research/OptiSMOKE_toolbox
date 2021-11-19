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

#ifndef OpenSMOKE_PerfectlyStirredReactor_Isothermal_ConstantPressure_H
#define	OpenSMOKE_PerfectlyStirredReactor_Isothermal_ConstantPressure_H

// Parent class
#include "PerfectlyStirredReactor.h"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"

// Options
#include "PerfectlyStirredReactor_Options.h"
#include "math/external-ode-solvers/ODE_Parameters.h"


namespace OpenSMOKE
{
	//!  A class for simulating perfectly stirred reactors with constant pressure in isothermal conditions
	/*!
		 The purpose of this class is to simulate a perfectly stirred reactor with constant pressure, in isothermal conditions
		 The conservation equations of species are solved in terms of mass fractions.
	*/

	class PerfectlyStirredReactor_Isothermal_ConstantPressure : public PerfectlyStirredReactor
	{

	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap map containing the thermodynamic data
		*@param kineticsMap map containing the kinetic mechanism
		*@param ode_parameters parameters governing the solution of the stiff ODE system
		*@param psr_options options governing the output
		*@param on_the_fly_ropa rate of production analysis (on the fly)
		*@param on_the_fly_post_processing post-processing analysis (on the fly)
		*@param T0 initial temperature [K]
		*@param P0 initial pressure [Pa]
		*@param omega0 initial mass fractions of species
		*@param TI initial temperature [K]
		*@param PI initial pressure [Pa]
		*@param omegaI initial mass fractions of species
		*@param constraint_a constraint: mass-flow-rate, volume, residence-time
		*@param value_a constraint value: mass-flow-rate [kg/s], volume [m3], residence-time [s]
		*@param constraint_b constraint: mass-flow-rate, volume, residence-time
		*@param valueba constraint value: mass-flow-rate [kg/s], volume [m3], residence-time [s]
		*/
		PerfectlyStirredReactor_Isothermal_ConstantPressure(	
								OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								OpenSMOKE::ODE_Parameters& ode_parameters,
								OpenSMOKE::PerfectlyStirredReactor_Options& psr_options,
								OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa,
								OpenSMOKE::OnTheFlyPostProcessing& on_the_fly_post_processing,
								OpenSMOKE::PolimiSoot_Analyzer& polimi_soot_analyzer,
								const double T0, const double P0, const OpenSMOKE::OpenSMOKEVectorDouble& omega0,
								const double TInlet, const double PInlet, const OpenSMOKE::OpenSMOKEVectorDouble& omegaI,
								const double residence_time, const double volume, const double mass_flow_rate);

		/**
		*@brief Solves the perfectly stirred reactor
		*@param tf the final time of integration [s]
		*/
		virtual void Solve(const double tf);

		/**
		*@brief Ordinary differential Equations corresponding to the reactor under investigation
		*@param t current time [s]
		*@param y current solution
		*@param dy current unsteady terms
		*/
		virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

		/**
		*@brief Writes the output (called at the end of each time step)
		*@param t current time [s]
		*@param y current solution
		*/
		virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

		/**
		*@brief Print the final status on output stream
		*@param fOutput output stream
		*/
		void PrintFinalStatus(std::ostream& fOutput);

		/**
		*@brief Print the final status on output stream
		*@param fOutput output stream
		*/
		void PrintParametricFinalStatus(std::ostream& fOutput);

		void PrintSolution();

	protected:

		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param t current time [s]
		*@param y current solution
		*/
		void SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);
	};
}

#include "PerfectlyStirredReactor_Isothermal_ConstantPressure.hpp"

#endif	/* OpenSMOKE_PerfectlyStirredReactor_Isothermal_ConstantPressure_H */

