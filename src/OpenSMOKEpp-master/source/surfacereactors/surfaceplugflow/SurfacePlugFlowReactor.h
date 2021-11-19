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

#ifndef OpenSMOKE_SurfacePlugFlowReactor_H
#define	OpenSMOKE_SurfacePlugFlowReactor_H

#include "math/OpenSMOKEVector.h"
#include "SurfacePlugFlowReactor_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysis_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"
#include "boost/filesystem.hpp"
#include "math/external-ode-solvers/ODE_Parameters.h"
#include "math/external-dae-solvers/DAE_Parameters.h"

namespace OpenSMOKE
{
	enum SurfacePlugFlowReactor_Type {	SURFACEPLUGFLOW_REACTOR_ISOTHERMAL,
										SURFACEPLUGFLOW_REACTOR_NONISOTHERMAL };

	//!  A class for simulating plug flow reactors
	/*!
		 The purpose of this class is to simulate plug flow reactors. It provides a common interface to different
		 plug flow reactors
	*/

	class SurfacePlugFlowReactor
	{

	public:

		/**
		*@brief Solves the equations governing the inlet conditions to the plug flow reactor
		*@param tInf the final time of integration (must be sufficiently long to reach steady-state) [s]
		*/
		virtual void SolveInletConditions(const double tInf, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& plugflow_options, OpenSMOKE::ODE_Parameters& ode_parameters);

		/**
		*@brief Solves the plug flow reactor
		*@param tf the final time of integration [s]
		*/
		virtual void Solve(const double tf) = 0;

		/**
		*@brief Moves the calculated inlet conditions as initial conditions for DAE solution
		*/
		virtual void SetInletConditions() = 0;

		/**
		*@brief Ordinary Differential Equations corresponding to the inlet conditions of reactor under investigation
		*@param t current time [s]
		*@param y current solution
		*@param dy current unsteady terms
		*/
		virtual int OdeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy) = 0;

		/**
		*@brief Algebraic Differential Equations corresponding to the reactor under investigation
		*@param t current time [s]
		*@param y current solution
		*@param dy current unsteady terms
		*/
		virtual int DaeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy) = 0;

		/**
		*@brief Writes the output during the solution of inlet conditions
		*@param t current time [s]
		*@param y current solution
		*/
		virtual int OdePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y) = 0;

		/**
		*@brief Writes the output (called at the end of each time step)
		*@param t current time [s]
		*@param y current solution
		*/
		virtual int DaePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y) = 0;

		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param t current time [s]
		*@param y current solution
		*/
		virtual void SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y) = 0;

		/**
		*@brief Return the total number of ODE equations
		*@return the total number of ODE equations
		*/
		unsigned int NumberOfOdeEquations() const { return NE_ODE_; };

		/**
		*@brief Return the total number DAE of equations
		*@return the total number of DAE equations
		*/
		unsigned int NumberOfDaeEquations() const { return NE_DAE_; };
		
		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param sensitivityMap map of sensitivity analysis
		*@param thermodynamicsMap map of thermodynamic data
		*@param sensitivity_options options for sensitivity analysis
		*/		
		void EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap& sensitivityMap, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& batch_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options);
	
		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param sensitivityMap map of sensitivity analysis
		*@param thermodynamicsMap map of thermodynamic data
		*@param sensitivity_options options for sensitivity analysis
		*/		
		void EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap& sensitivityMap, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& batch_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options, bool iXmlMaps);
        
		/**
		*@brief Writes on a file the final summary of the plugflow reactor
		*@param summary_file file where the summary will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param plugflow_options options for solving the plug flow reactor
		*/		
		void DaeFinalSummary(const boost::filesystem::path summary_file, const double tf, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap, SurfacePlugFlowReactor_Options& plugflow_options);

		/**
		*@brief Returns the final composition of the plugflow reactor, together with temperature and pressure
		*@param T temperature [K]
		*@param P pressure [Pa]
		*@param G specific mass flow rate [kg/s/m2]
		*@param omega mass fractions [-]
		*/	
		void DaeGetFinalStatus(double& T, double& P, double& G, OpenSMOKE::OpenSMOKEVectorDouble& omega, OpenSMOKE::OpenSMOKEVectorDouble& Z, OpenSMOKE::OpenSMOKEVectorDouble& Gamma, OpenSMOKE::OpenSMOKEVectorDouble& thicknessBulk);

		/**
		*@brief Writes on the video the final status of the plugflow reactor
		*@param tf final time of integration
		*@param thermodynamicsMap map of thermodynamic data
		*@param plugflow_options options for solving the plug flow reactor
		*/			
		void DaeFinalStatus(const double tf, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap, SurfacePlugFlowReactor_Options& plugflow_options);

		/**
		*@brief Writes on a file the calculated inlet status of the plugflow reactor
		*@param summary_file file where the summary will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param plugflow_options options for solving the plug flow reactor
		*/
		void OdeFinalSummary(const boost::filesystem::path summary_file, const double tf, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap, SurfacePlugFlowReactor_Options& plugflow_options);

		/**
		*@brief Prepares the ASCII file
		*@param fOutput ASCII file where the output will be written
		*@param output_file_ascii name of ASCII file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param plugflow_options options for solving the plug flow reactor
		*/
		void PrepareASCIIFile(std::ofstream& fOutput, boost::filesystem::path output_file_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& plugflow_options);

		/**
		*@brief Prepares the ASCII file for bulk species
		*@param output_file_bulk_ascii ASCII file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param batch_options options for solving the batch reactor
		*/
		void PrepareBulkASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_bulk_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& plugflow_options);

		/**
		*@brief Prepares the ASCII file
		*@param fOutput ASCII file where the output will be written
		*@param output_file_ascii name of ASCII file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param plugflow_options options for solving the plug flow reactor
		*/
		void PrepareParametricASCIIFile(std::ofstream& fOutput, boost::filesystem::path output_file_ascii, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& plugflow_options);

	protected:

		/**
		*@brief Prepares the ASCII file
		*@param output_file_ascii ASCII file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param plugflow_options options for solving the plug flow reactor
		*/
		void PrepareASCIIFile(const boost::filesystem::path output_file_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& plugflow_options);

		/**
		*@brief Prepares the ASCII file
		*@param output_file_ascii ASCII file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param plugflow_options options for solving the plug flow reactor
		*/
		void PrepareBulkASCIIFile(const boost::filesystem::path output_file_bulk_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& plugflow_options);

		/**
		*@brief Prepares the XML file
		*@param output_file_xml XML file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*/
		void PrepareXMLFile(const boost::filesystem::path output_file_xml, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap);
		
		/**
		*@brief Closes the XML file
		*/		
		void CloseXMLFile();

		/**
		*@brief Calculates the numerical Jacobian
		*@param t current time
		*@param y current solution
		*@param J Jacobian matrix
		*/
		void NumericalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEMatrixDouble& J);	
		
		/**
		*@brief Allocates the memory for the vectors and the matrices contained in the plugflow reactor
		*/			
		void MemoryAllocation();

		/**
		*@brief Prepares the XML files for the sensitivity analysis
		*@param thermodynamicsMap map of thermodynamic data
		*@param plugflow_options options for solving the plugflow reactor
		*/
		void PrepareSensitivityXMLFiles(OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& plugflow_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options);
		
		/**
		*@brief Closes the XML files for the sensitivity analysis
		*/			
		void CloseSensitivityXMLFiles();

		/**
		*@brief Open all the files (ASCII and XML)
		*@param thermodynamicsMap map of thermodynamic data
		*@param plugflow_options options for solving the plugflow reactor
		*/
		void OpenAllFiles(OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, OpenSMOKE::SurfacePlugFlowReactor_Options& plugflow_options);
		
		/**
		*@brief Closes all the files (ASCII and XML)
		*@param plugflow_options options for solving the plugflow reactor
		*/		
		void CloseAllFiles(OpenSMOKE::SurfacePlugFlowReactor_Options& plugflow_options);

		/**
		*@brief Solves the Ordinary Differential Equations using one of the provided open/source ODE solvers
		*@param tf final time of integrations
		*@param ode_parameters parameters governing the ODE integration
		*/	
		void OdeSolveOpenSourceSolvers(const double tf, OpenSMOKE::ODE_Parameters& ode_parameters);
		
		/**
		*@brief Solves the Differential-Algebraic Equation System using one of the provided open/source DAE solvers
		*@param tf final time of integrations
		*@param algebraic_equations true is the equation is algebraic
		*@param dae_parameters parameters governing the DAE integration
		*/			
		void DaeSolveOpenSourceSolvers(const double tf, const bool* algebraic_equations, OpenSMOKE::DAE_Parameters& dae_parameters);

	
protected:

		SurfacePlugFlowReactor_Type type_;			//!< type of plugflow reactor

		double rho0_;								//!< initial density [kg/m3]
		double T0_;									//!< initial temperature [K]
		double P0_;									//!< initial pressure [Pa]
		double v0_;									//!< initial velocity [m/s]
		double G0_;									//!< initial mass flow rate per unit of surface [kg/m2/s]
		double thicknessBulk0_;						//!< initial thickness of bulk phase [kg] (always equal to 0, to be considered in a relative way)
		double MW0_;								//!< initial molecular weight [kg/kmol]
		double H0_;									//!< initial enthalpy [J/kmol]
		double H0_Surface_over_Volume_;				//!< initial enthalpy of surface per unit of volume [J/m]
		double U0_;									//!< initial internal energy [J]
		double U0_Surface_over_Volume_;				//!< initial internal energy of surface per unit of volume [J/m]
		OpenSMOKE::OpenSMOKEVectorDouble omega0_;	//!< initial composition [mass fractions]
		OpenSMOKE::OpenSMOKEVectorDouble Z0_;		//!< initial composition [surface fractions]
		OpenSMOKE::OpenSMOKEVectorDouble Gamma0_;	//!< initial surface densities [kmol/m2]
		OpenSMOKE::OpenSMOKEVectorDouble x0_;		//!< initial composition [mole fractions]

		bool time_independent_;
		bool constant_pressure_;

		double global_thermal_exchange_coefficient_;	//!< global thermal exchange coefficient [W/m2/K]
		double cross_section_over_perimeter_;           //!< ratio between cross section and perimeter [m]
		double T_environment_;							//!< environment temperature                
        
		double artificial_speedup_coefficient_;			//!< artificial speedup coefficient for algebraic equations
        
		std::vector<bool> site_non_conservation_;		//!< non conservation of sites is allowed (if true)

		unsigned int NC_;		//!< number of species [-]
		unsigned int NR_;		//!< number of reactions [-]
		unsigned int NE_ODE_;	//!< number of equations for the ODE system (inlet conditions) [-]
		unsigned int NE_DAE_;	//!< number of equations for the DAE system (whole plug flow reactor) [-]

		unsigned int SURF_NC_;	//!< number of surface species [-]
		unsigned int SURF_NR_;	//!< number of surface reactions [-]
		unsigned int SURF_NP_;	//!< number of surface phases [-]

		unsigned int BULK_NC_;	//!< number of surface species [-]
		unsigned int BULK_NP_;	//!< number of surface phases [-]

		OpenSMOKE::OpenSMOKEVectorDouble omega_;			//!< current mass fractions
		OpenSMOKE::OpenSMOKEVectorDouble Z_;				//!< current surface fractions
		OpenSMOKE::OpenSMOKEVectorDouble Gamma_;			//!< current surface densities [kmol/m2]
		OpenSMOKE::OpenSMOKEVectorDouble GammaFromEqn_;		//!< surface densities from equation [kmol/m2]
		OpenSMOKE::OpenSMOKEVectorDouble a_;				//!< current bulk species activities [-] (TODO: they are always assumed to be equal to 1)
		OpenSMOKE::OpenSMOKEVectorDouble thicknessBulk_;	//!< current cumulative bulk thickness [m]
		OpenSMOKE::OpenSMOKEVectorDouble x_;				//!< current mole fractions
		OpenSMOKE::OpenSMOKEVectorDouble c_;				//!< current concentrations [kmol/m3]
		OpenSMOKE::OpenSMOKEVectorDouble RfromGas_;			//!< current formation rates (from homogeneous reactions)  [kmol/m3/s]
		OpenSMOKE::OpenSMOKEVectorDouble RfromSurface_;		//!< current formation rates (from heterogeneous reations) [kmol/m3/s]
		OpenSMOKE::OpenSMOKEVectorDouble Rsurface_;			//!< current surface species formation rates [kmol/m2/s]
		OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases_;	//!< current surface site phases formation rates [kmol/m2/s]
		OpenSMOKE::OpenSMOKEVectorDouble Rbulk_;			//!< current bulk species formation rates [kmol/m2/s]

		double tau_;								//!< current residence time [s]
		double csi_;								//!< current axial coordinate [m]
		double rho_;								//!< current density [kg/m3]
		double T_;									//!< current temperature [K]
		double P_;									//!< current pressure [Pa]
		double v_;									//!< current velocity [m/s]
		double G_;									//!< mass flow rate per unit of surface [kg/m2/s]
		double cTot_;								//!< current concentration [kmol/m3]
		double MW_;									//!< current molecular weight [kg/kmol]
		double QRGas_;								//!< current heat release [W/m3]
		double QRSurface_;							//!< current heat release [W/m3]
		double CvMixMass_;							//!< current specific heat (constant volume) [J/kg/K]
		double CpMixMass_;							//!< current specific heat (constant pressure) [J/kg/K]

		unsigned int ode_iteration_;				//!< iteration index (ODE)
		unsigned int dae_iteration_;				//!< iteration index (DAE)
		unsigned int counter_file_video_;			//!< iteration counter for writing on video
		unsigned int counter_file_ASCII_;			//!< iteration counter for writing on ASCII file
		unsigned int counter_file_XML_;				//!< iteration counter for writing on XML file
		unsigned int counter_sensitivity_XML_;		//!< iteration counter for writing on XML sensitivity files

		OpenSMOKE::OpenSMOKEVectorDouble yOde0_;	//!< vector containing the initial values for all the ODE variables
		OpenSMOKE::OpenSMOKEVectorDouble yOdef_;	//!< vector containing the final values for all the ODE variables

		OpenSMOKE::OpenSMOKEVectorDouble yDae0_;	//!< vector containing the initial values for all the DAE variables
		OpenSMOKE::OpenSMOKEVectorDouble yDaef_;	//!< vector containing the final values for all the DAE variables
		OpenSMOKE::OpenSMOKEVectorDouble dyDae0_;	//!< vector containing the initial values for all the DAE variables
		OpenSMOKE::OpenSMOKEVectorDouble dyDaef_;	//!< vector containing the final values for all the DAE variables

		std::ofstream fXML_;						//!< XML file where the output is written
		std::ofstream fASCII_;						//!< ASCII file where the output is written
		std::ofstream fBulkASCII_;					//!< ASCII file where the output for bulk species is written
		std::ofstream  fSensitivityParentXML_;		//!< XML file where the sensitivity analysis is written (parent file)
		std::ofstream* fSensitivityChildXML_;		//!< XML files where the sensitivity analysis is written (children files)

		std::vector<unsigned int> indices_of_output_species_;			//!< indices of species for which the profiles are written on the ASCII file
		std::vector<unsigned int> indices_of_sensitivity_species_;		//!< indices of species for which the sensitivity analysis is written on file
		std::vector<unsigned int> widths_of_output_species_;			//!< width   of species for which the profiles are written on the ASCII file
		std::vector<unsigned int> widths_of_output_surface_species_;	//!< width   of surface species for which the profiles are written on the ASCII file

		OpenSMOKE::SensitivityMap* sensitivityMap_;				//!< sensitivity map
		OpenSMOKE::OpenSMOKEVectorDouble scaling_Jp;					//!< vector containing datafor performing the sensitivity analysis
		OpenSMOKE::OpenSMOKEMatrixDouble J;								//!< Jacobian matrix (required for sensitivity analysis)

		bool iXmlMaps_;													//!< Activate Xml maps writing (required for sensitivity analysis)
	};
}

#include "SurfacePlugFlowReactor.hpp"

#endif	/* OpenSMOKE_SurfacePlugFlowReactor_H */

