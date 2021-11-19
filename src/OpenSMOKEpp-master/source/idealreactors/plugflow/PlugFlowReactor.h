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

#ifndef OpenSMOKE_PlugFlowReactor_H
#define	OpenSMOKE_PlugFlowReactor_H

#include "math/OpenSMOKEVector.h"
#include "PlugFlowReactor_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysis_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"
#include "boost/filesystem.hpp"
#include "math/external-ode-solvers/ODE_Parameters.h"

namespace OpenSMOKE
{
	enum PlugFlowReactor_Type { PLUGFLOW_REACTOR_ISOTHERMAL,
								PLUGFLOW_REACTOR_NONISOTHERMAL };

	//!  A class for simulating plug flow reactors
	/*!
		 The purpose of this class is to simulate plug flow reactors. It provides a common interface to different
		 plug flow reactors
	*/

	class PlugFlowReactor
	{

	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap map containing the thermodynamic data
		*@param kineticsMap map containing the kinetic mechanism
		*@param ode_parameters parameters governing the solution of the stiff ODE system
		*@param plugflow_options options governing the output
		*@param on_the_fly_ropa rate of production analysis (on the fly)
		*@param on_the_fly_post_processing post-processing analysis (on the fly)
		*@param idts_analyzer ignition delay time analyzer (on the fly)
		*@param polimi_soot_analyzer soot analyzer (on the fly)
		*/
		PlugFlowReactor(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
							OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
							OpenSMOKE::ODE_Parameters& ode_parameters,
							OpenSMOKE::PlugFlowReactor_Options& plugflow_options,
							OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa,
							OpenSMOKE::OnTheFlyPostProcessing& on_the_fly_post_processing,
							OpenSMOKE::IgnitionDelayTimes_Analyzer& idts_analyzer,
							OpenSMOKE::PolimiSoot_Analyzer& polimi_soot_analyzer);

		/**
		*@brief Solves the plug flow reactor
		*@param tf the final time of integration [s]
		*/
		virtual void Solve(const double tf) = 0;

		/**
		*@brief Ordinary differential Equations corresponding to the reactor under investigation
		*@param t current time [s]
		*@param y current solution
		*@param dy current unsteady terms
		*/
		virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy) = 0;

		/**
		*@brief Writes the output (called at the end of each time step)
		*@param t current time [s]
		*@param y current solution
		*/
		virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y) = 0;

		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param t current time [s]
		*@param y current solution
		*/
		virtual void SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y) = 0;

		/**
		*@brief Return the total number of equations
		*@return the total number of equations
		*/
		unsigned int NumberOfEquations() const { return NE_; };

		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param sensitivityMap map of sensitivity analysis
		*@param sensitivity_options options for sensitivity analysis
		*/		
		void EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap& sensitivityMap, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options);
	
		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param sensitivityMap map of sensitivity analysis
		*@param sensitivity_options options for sensitivity analysis
		*/		
		void EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap& sensitivityMap, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options, bool iXmlMaps);
        
		/**
		*@brief Writes on a file the final summary of the plugflow reactor
		*@param summary_file file where the summary will be written
		*/		
		void FinalSummary(const boost::filesystem::path summary_file, const double tf);

		/**
		*@brief Returns the final composition of the plugflow reactor, together with temperature and pressure
		*@param T temperature [K]
		*@param P pressure [Pa]
		*@param omega mass fractions [-]
		*/	
		void GetFinalStatus(double& T, double& P, OpenSMOKE::OpenSMOKEVectorDouble& omega);

		/**
		*@brief Writes on the video the final status of the plugflow reactor
		*@param tf final time of integration
		*/			
		void FinalStatus(const double tf);

		/**
		*@brief Prepares the ASCII file
		*@param fOutput ASCII file where the output will be written
		*@param output_file_ascii name of ASCII file where the output will be written
		*/
		void PrepareASCIIFile(std::ofstream& fOutput, boost::filesystem::path output_file_ascii);

		/**
		*@brief Prepares the ASCII file
		*@param fOutput ASCII file where the output will be written
		*@param output_file_ascii name of ASCII file where the output will be written
		*/
		void PrepareParametricASCIIFile(std::ofstream& fOutput, boost::filesystem::path output_file_ascii);

		/**
		*@brief Print Polimi Soot Analysis on file
		*@param t current time [s]
		*/
		void PrintPolimiSoot(const double t);

		/**
		*@brief Prepares the ASCII file for ignition delay times
		*@param fOutput ASCII file where the output will be written
		*@param output_file_ascii name of ASCII file where the output will be written
		*/
		void PrepareParametricIDTASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii);

		/**
		*@brief Print the ignition delay times on output stream
		*@param fOutput output stream
		*@param t time
		*/
		void PrintParametricIDT(std::ostream& fOutput, const double t);

	protected:

		/**
		*@brief Prepares the ASCII file
		*@param output_file_ascii ASCII file where the output will be written
		*/
		void PrepareASCIIFile(const boost::filesystem::path output_file_ascii);

		/**
		*@brief Prepares the ASCII file for ROPA (on the fly)
		*@param output_file_ropa name of ASCII file where the ROPA will be written
		*/
		void PrepareROPAFile(const boost::filesystem::path output_file_ropa);

		/**
		*@brief Prepares the ASCII file for soot analysis (Polimi)
		*@param output_file_polimi_soot name of ASCII file where the Polimi Soot Analysis will be written
		*@param output_file_polimi_soot_distribution name of ASCII file where the Polimi Soot Distribution will be written
		*/
		void PreparePolimiSootFiles(const boost::filesystem::path output_file_polimi_soot, const boost::filesystem::path output_file_polimi_soot_distribution);

		/**
		*@brief Prepares the XML file
		*@param output_file_xml XML file where the output will be written
		*/
		void PrepareXMLFile(const boost::filesystem::path output_file_xml);
		
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
		*/
		void PrepareSensitivityXMLFiles(OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options);
		
		/**
		*@brief Closes the XML files for the sensitivity analysis
		*/			
		void CloseSensitivityXMLFiles();

		/**
		*@brief Open all the files (ASCII and XML)
		*/
		void OpenAllFiles();
		
		/**
		*@brief Closes all the files (ASCII and XML)
		*/		
		void CloseAllFiles();

		/**
		*@brief Solves the Ordinary Differential Equations using one of the provided open/source ODE solvers
		*@param tf final time of integrations
		*/	
		void SolveOpenSourceSolvers(const double tf);

		/**
		*@brief Performs the Polimi Soot Analysis
		*@param t current time [s]
		*/
		void PolimiSootAnalysis(const double t);


protected:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;				//!< thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&			kineticsMap_;					//!< kinetic map
		OpenSMOKE::ODE_Parameters&				ode_parameters_;				//!< ode parameters
		OpenSMOKE::PlugFlowReactor_Options&		plugflow_options_;				//!< options
		OpenSMOKE::OnTheFlyROPA&				on_the_fly_ropa_;				//!< on the fly ROPA
		OpenSMOKE::OnTheFlyPostProcessing&		on_the_fly_post_processing_;	//!< on the fly post-processing
		OpenSMOKE::PolimiSoot_Analyzer&			polimi_soot_analyzer_;			//!< on the fly soot analyzer (based on POLIMI)
		OpenSMOKE::IgnitionDelayTimes_Analyzer&			idts_analyzer_;					//!< on the fly signition delay time analyzer

protected:

		PlugFlowReactor_Type type_;					//!< type of plugflow reactor

		double rho0_;								//!< initial density [kg/m3]
		double T0_;									//!< initial temperature [K]
		double P0_;									//!< initial pressure [Pa]
		double v0_;									//!< initial velocity [m/s]
		double MW0_;								//!< initial molecular weight [kg/kmol]
		double H0_;									//!< initial enthalpy [J/kmol]
		double U0_;									//!< initial internal energy [J/kmol]
		OpenSMOKE::OpenSMOKEVectorDouble omega0_;	//!< initial composition [mass fractions]
		OpenSMOKE::OpenSMOKEVectorDouble x0_;		//!< initial composition [mole fractions]

		bool time_independent_;
		bool constant_pressure_;
		double specificmassflowrate_;				//!< mass [kg]

		double global_thermal_exchange_coefficient_;	//!< global thermal exchange coefficient [W/m2/K]
		double cross_section_over_perimeter_;           //!< ratio between cross section and perimeter [m]
		double T_environment_;				//!< environment temperature                
                
                
		unsigned int NC_;							//!< number of species [-]
		unsigned int NR_;							//!< number of reactions [-]
		unsigned int NE_;							//!< number of equations [-]

		OpenSMOKE::OpenSMOKEVectorDouble omega_;	//!< current mass fractions
		OpenSMOKE::OpenSMOKEVectorDouble x_;		//!< current mole fractions
		OpenSMOKE::OpenSMOKEVectorDouble c_;		//!< current concentrations [kmol/m3]
		OpenSMOKE::OpenSMOKEVectorDouble R_;		//!< current formation rates [kmol/m3/s]
		double tau_;								//!< current residence time [s]
		double csi_;								//!< current axial coordinate [m]
		double rho_;								//!< current density [kg/m3]
		double T_;									//!< current temperature [K]
		double P_;									//!< current pressure [Pa]
		double v_;									//!< current velocity [m/s]
		double cTot_;								//!< current concentration [kmol/m3]
		double MW_;									//!< current molecular weight [kg/kmol]
		double QR_;									//!< current heat release [W/m3]
		double CvMixMass_;							//!< current specific heat (constant volume) [J/kg/K]
		double CpMixMass_;							//!< current specific heat (constant pressure) [J/kg/K]

		unsigned int iteration_;					//!< iteration index
		unsigned int counter_file_video_;			//!< iteration counter for writing on video
		unsigned int counter_file_ASCII_;			//!< iteration counter for writing on ASCII file
		unsigned int counter_file_XML_;				//!< iteration counter for writing on XML file
		unsigned int counter_sensitivity_XML_;		//!< iteration counter for writing on XML sensitivity files

		OpenSMOKE::OpenSMOKEVectorDouble y0_;		//!< vector containing the initial values for all the variables
		OpenSMOKE::OpenSMOKEVectorDouble yf_;		//!< vector containing the final values for all the variables

		std::ofstream fXML_;						//!< XML file where the output is written
		std::ofstream fASCII_;						//!< ASCII file where the output is written
		std::ofstream fROPA_;						//!< ASCII file where the ROPA (on the fly) is written
		std::ofstream fPolimiSoot_;					//!< ASCII file where the Polimi Soot Analysis is written
		std::ofstream fPolimiSootDistribution_;		//!< ASCII file where the Polimi Soot Distribution is written
		std::ofstream  fSensitivityParentXML_;		//!< XML file where the sensitivity analysis is written (parent file)
		std::ofstream* fSensitivityChildXML_;		//!< XML files where the sensitivity analysis is written (children files)

		std::vector<unsigned int> indices_of_output_species_;		//!< indices of species for which the profiles are written on the ASCII file
		std::vector<unsigned int> indices_of_sensitivity_species_;	//!< indices of species for which the sensitivity analysis is written on file
		std::vector<unsigned int> widths_of_output_species_;        //!< width   of species for which the profiles are written on the ASCII file
		
		OpenSMOKE::SensitivityMap* sensitivityMap_;			//!< sensitivity map
		OpenSMOKE::OpenSMOKEVectorDouble scaling_Jp;		//!< vector containing datafor performing the sensitivity analysis
		OpenSMOKE::OpenSMOKEMatrixDouble J;					//!< Jacobian matrix (required for sensitivity analysis)

		std::vector<std::string>	reaction_names_;	//!< Reaction names
		bool iXmlMaps_;									//!< Activate Xml maps writing (required for sensitivity analysis)

		// AB
public:
		std::vector<double> time_vector;

		std::vector<OpenSMOKE::OpenSMOKEVectorDouble> species_matrix;
		//AB
	};
}

#include "PlugFlowReactor.hpp"

#endif	/* OpenSMOKE_PlugFlowReactor_H */

