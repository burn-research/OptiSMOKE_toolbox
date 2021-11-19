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

#ifndef OpenSMOKE_ShockTubeReactor_CHEMKIN_Model_H
#define	OpenSMOKE_ShockTubeReactor_CHEMKIN_Model_H

// Parent class
#include "ShockTubeReactor.h"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"

// Options
#include "ShockTubeReactor_Options.h"
#include "math/external-ode-solvers/ODE_Parameters.h"


namespace OpenSMOKE
{
	//!  A class for simulating batch reactors with constant pressure in adiabatic conditions
	/*!
		 The purpose of this class is to simulate a batch reactor with constant pressure, in adiabatic conditions
		 The conservation equations of species and energy are solved in terms of mass fractions and
		 temperature respectively. The internal energy (massive) is conserved
	*/

	class ShockTubeReactor_CHEMKIN_Model : public ShockTubeReactor
	{

	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap map containing the thermodynamic data
		*@param kineticsMap map containing the kinetic mechanism
		*@param ode_parameters parameters governing the solution of the stiff ODE system
		*@param shocktube_options options governing the output
		*@param on_the_fly_ropa rate of production analysis (on the fly)
		*@param on_the_fly_post_processing post-processing analysis (on the fly)
		*@param v0 initial velocity [m/s]
		*@param T0 initial temperature [K]
		*@param P0 initial pressure [Pa]
		*@param omega0 initial mass fractions of species
		*@param d0 initial diameter of the cross sectional area [m]
		*@param lm boundary layer characteristic length [m]
		*@param rho1 density of the mixture before the shock [kg/m3/s]
		*/
		ShockTubeReactor_CHEMKIN_Model(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
										OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
										OpenSMOKE::ODE_Parameters& ode_parameters,
										OpenSMOKE::ShockTubeReactor_Options& shocktube_options,
										OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa,
										OpenSMOKE::OnTheFlyPostProcessing& on_the_fly_post_processing,
										OpenSMOKE::PolimiSoot_Analyzer& polimi_soot_analyzer,
										const double v0, const double T0, const double P0, 
										const OpenSMOKE::OpenSMOKEVectorDouble& omega0, const double rho1,
										const double d0, const double lm );

		/**
		*@brief Solves the batch reactor
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

	protected:

		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param t current time [s]
		*@param y current solution
		*/
		void SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);
	};
}

#include "ShockTubeReactor_CHEMKIN_Model.hpp"

#endif	/* OpenSMOKE_ShockTubeReactor_CHEMKIN_Model_H */

