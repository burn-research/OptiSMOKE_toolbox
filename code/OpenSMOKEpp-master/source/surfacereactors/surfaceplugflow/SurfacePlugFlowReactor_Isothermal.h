/*-----------------------------------------------------------------------*\
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
|   Copyright(C) 2016, 2015 Alberto Cuoci                                 |
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

#ifndef OpenSMOKE_SurfaceSurfacePlugFlowReactor_Isothermal_H
#define	OpenSMOKE_SurfacePlugFlowReactor_Isothermal_H

// Parent class
#include "SurfacePlugFlowReactor.h"
#include "SurfacePlugFlowReactor_Profile.h"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"

// Options
#include "SurfacePlugFlowReactor_Options.h"
#include "math/external-ode-solvers/ODE_Parameters.h"


namespace OpenSMOKE
{
	//!  A class for simulating plug flow reactors with constant pressure in isothermal conditions
	/*!
		 The purpose of this class is to simulate a plug flow reactor with constant pressure, in isothermal conditions
		 The conservation equations of species are solved in terms of mass fractions.
	*/

	class SurfacePlugFlowReactor_Isothermal : public SurfacePlugFlowReactor
	{

	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap map containing the thermodynamic data
		*@param kineticsMap map containing the kinetic mechanism
		*@param ode_parameters parameters governing the solution of the stiff ODE system (consistent inlet conditions)
		*@param dae_parameters parameters governing the solution of the stiff DAE system
		*@param plugflow_options options governing the output
		*@param v0 initial velocity [m/s]
		*@param T0 initial temperature [K]
		*@param P0 initial pressure [Pa]
		*@param omega0 initial mass fractions of species
		*/
		SurfacePlugFlowReactor_Isothermal(		OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
												OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
												OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap,
												OpenSMOKE::KineticsMap_Surface_CHEMKIN& kineticsSurfaceMap,
												OpenSMOKE::ODE_Parameters& ode_parameters,
												OpenSMOKE::DAE_Parameters& dae_parameters,
												OpenSMOKE::SurfacePlugFlowReactor_Options& plugflow_options,
												const bool time_independent, 
												const bool constant_pressure, 
												const double v0, const double T0, const double P0, 
												const OpenSMOKE::OpenSMOKEVectorDouble& omega0,
												const OpenSMOKE::OpenSMOKEVectorDouble& Z0,
												const OpenSMOKE::OpenSMOKEVectorDouble& Gamma0,
												const std::vector<bool>& site_non_conservation,
												const double cross_section_over_perimeter );

		/**
		*@brief Moves the calculated inlet conditions as initial conditions for DAE solution
		*/
		virtual void SetInletConditions();

		/**
		*@brief Solves the plugflow reactor
		*@param tf the final time of integration [s]
		*/
		virtual void Solve(const double tf);

		/**
		*@brief Ordinary Differential Equations corresponding to the inlet conditions of reactor under investigation
		*@param t current time [s]
		*@param y current solution
		*@param dy current unsteady terms
		*/
		virtual int OdeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

		/**
		*@brief Algebraic Differential Equations corresponding to the reactor under investigation
		*@param t current time [s]
		*@param y current solution
		*@param dy current unsteady terms
		*/
		virtual int DaeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

		/**
		*@brief Writes the output during the solution of inlet conditions
		*@param t current time [s]
		*@param y current solution
		*/
		virtual int OdePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

		/**
		*@brief Writes the output (called at the end of each time step)
		*@param t current time [s]
		*@param y current solution
		*/
		virtual int DaePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

		/**
		*@brief Set a temperature profile
		*@param profile temperature profile
		*/		
		void SetTemperatureProfile(OpenSMOKE::SurfacePlugFlowReactor_Profile&  profile);

		/**
		*@brief Print the current status in the plug flow reactor
		*@param fOutput output stream
		*@param fOutput output stream for bulk status
		*/
		void PrintCurrentStatus(std::ostream& fOutput, std::ostream& fBulkOutput);

		/**
		*@brief Print the final status on output stream
		*@param fOutput output stream
		*/
		void PrintParametricFinalStatus(std::ostream& fOutput);

	protected:

		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param t current time [s]
		*@param y current solution
		*/
		virtual void SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);
                
        bool fixed_temperature_profile_;
		OpenSMOKE::SurfacePlugFlowReactor_Profile*  temperature_profile_;

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&			thermodynamicsMap_;				//!< thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&					kineticsMap_;					//!< kinetic map
		OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN&	thermodynamicsSurfaceMap_;		//!< surface thermodynamic map
		OpenSMOKE::KineticsMap_Surface_CHEMKIN&			kineticsSurfaceMap_;			//!< surface kinetic map
		OpenSMOKE::ODE_Parameters&						ode_parameters_;				//!< ode parameters
		OpenSMOKE::DAE_Parameters&						dae_parameters_;				//!< dae parameters
		OpenSMOKE::SurfacePlugFlowReactor_Options&		plugflow_options_;				//!< options

	};
}

#include "SurfacePlugFlowReactor_Isothermal.hpp"

#endif	/* OpenSMOKE_SurfacePlugFlowReactor_Isothermal_H */

