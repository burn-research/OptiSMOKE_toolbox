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

#ifndef OpenSMOKE_SurfacePerfectlyStirredReactor_H
#define	OpenSMOKE_SurfacePerfectlyStirredReactor_H

#include "math/OpenSMOKEVector.h"
#include "SurfacePerfectlyStirredReactor_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysis_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"
#include "boost/filesystem.hpp"
#include "math/external-ode-solvers/ODE_Parameters.h"

namespace OpenSMOKE
{
	enum SurfacePerfectlyStirredReactor_Type {	SURFACEPERFECTLYSTIRRED_REACTOR_NONISOTHERMAL_CONSTANTP,
												SURFACEPERFECTLYSTIRRED_REACTOR_ISOTHERMAL_CONSTANTP };

	enum SurfacePerfectlyStirredReactor_ConstraintType { MASSFLOWRATE_RESIDENCETIME, MASSFLOWRATE_VOLUME, VOLUME_RESIDENCETIME};

	//!  A class for simulating perfectly stirred reactors
	/*!
		 The purpose of this class is to simulate perfectly stirred reactors. It provides a common interface to different
		 perfectly stirred reactors
	*/

	class SurfacePerfectlyStirredReactor
	{

	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap map containing the thermodynamic data
		*@param kineticsMap map containing the kinetic mechanism
		*@param thermodynamicsSurfaceMap map containing the surface thermodynamic data
		*@param kineticsSurfaceMap map containing the surface kinetic mechanism
		*@param ode_parameters parameters governing the solution of the stiff ODE system
		*@param psr_options options governing the output
		*@param on_the_fly_ropa rate of production analysis (on the fly)
		*/
		SurfacePerfectlyStirredReactor(
			OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
			OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
			OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap,
			OpenSMOKE::KineticsMap_Surface_CHEMKIN& kineticsSurfaceMap,
			OpenSMOKE::ODE_Parameters& ode_parameters,
			OpenSMOKE::SurfacePerfectlyStirredReactor_Options& psr_options,
			OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa);

		/**
		*@brief Solves the perfectly stirred reactor
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
		*@brief Writes on the video the final status of the reactor
		*@param tf final time of integration
		*/			
		void FinalStatus(const double tf);

		/**
		*@brief Writes on a file the final summary of the reactor
		*@param summary_file file where the summary will be written
		*/		
		void FinalSummary(const boost::filesystem::path summary_file, const double tf);

		/**
		*@brief Returns the final composition of the reactor, together with temperature and pressure
		*@param T temperature [K]
		*@param P pressure [Pa]
		*@param omega mass fractions [-]
		*/	
		void GetFinalStatus(double& T, double& P, OpenSMOKE::OpenSMOKEVectorDouble& omega, OpenSMOKE::OpenSMOKEVectorDouble& Z, OpenSMOKE::OpenSMOKEVectorDouble& Gamma, double& mass, OpenSMOKE::OpenSMOKEVectorDouble& massBulk);

		/**
		*@brief Prepares the ASCII file
		*@param fOutput ASCII file where the output will be written
		*@param output_file_ascii name of ASCII file where the output will be written
		*/
		void PrepareASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii);

		/**
		*@brief Prepares the ASCII file for bulk species
		*@param fOutput ASCII file where the output will be written
		*@param output_file_ascii name of ASCII file where the output will be written
		*/
		void PrepareBulkASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii);

	protected:
	
		/**
		*@brief Prepares the ASCII file for ROPA (on the fly)
		*@param output_file_ropa name of ASCII file where the ROPA will be written
		*/
		void PrepareROPAFile(const boost::filesystem::path output_file_ropa);

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
		*@brief Allocates the memory for the vectors and the matrices contained in the perfectly stirred reactor
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

protected:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&			thermodynamicsMap_;			//!< thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN&					kineticsMap_;				//!< kinetic map
		OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN&	thermodynamicsSurfaceMap_;	//!< surface thermodynamic map
		OpenSMOKE::KineticsMap_Surface_CHEMKIN&			kineticsSurfaceMap_;		//!< surface kinetic map
		OpenSMOKE::ODE_Parameters&								ode_parameters_;			//!< ode parameters
		OpenSMOKE::SurfacePerfectlyStirredReactor_Options&		psr_options_;				//!< options
		OpenSMOKE::OnTheFlyROPA&								on_the_fly_ropa_;			//!< on the fly ROPA
	
protected:

		SurfacePerfectlyStirredReactor_Type type_;							//!< type of psr
		SurfacePerfectlyStirredReactor_ConstraintType constraint_type_;	//!< constraint type

		double T0_;										//!< initial temperature [K]
		double P0_;										//!< initial pressure [Pa]
		double MW0_;									//!< initial molecular weight [kg/kmol]
		double H0_;										//!< initial enthalpy [J/kmol]
		double U0_;										//!< initial internal energy [J/kmol]
		double rho0_;									//!< initial density [kg/m3]
		double mass0_;									//!< initial mass [kg]
		OpenSMOKE::OpenSMOKEVectorDouble omega0_;		//!< initial composition [mass fractions]
		OpenSMOKE::OpenSMOKEVectorDouble x0_;			//!< initial composition [mole fractions]
		OpenSMOKE::OpenSMOKEVectorDouble Z0_;			//!< initial surface fractions
		OpenSMOKE::OpenSMOKEVectorDouble Gamma0_;		//!< initial surface densities [kmol/m2]
		OpenSMOKE::OpenSMOKEVectorDouble massBulk0_;	//!< initial mass of bulk species [kg]

		double TInlet_;									//!< inlet temperature [K]
		double PInlet_;									//!< inlet pressure [Pa]
		double MWInlet_;								//!< inlet molecular weight [kg/kmol]
		double HInlet_;								//!< inlet enthalpy [J/kmol]
		double UInlet_;									//!< inlet internal energy [J/kmol]
		double rhoInlet_;								//!< inlet density [kg/m3]
		OpenSMOKE::OpenSMOKEVectorDouble omegaInlet_;	//!< inlet composition [mass fractions]
		OpenSMOKE::OpenSMOKEVectorDouble xInlet_;		//!< inlet composition [mole fractions]

		double rho_;				//!< density [kg/m3]
		double mass_;				//!< mass [kg]

		double V0_;								//!< initial volume [m3]
		double V_;								//!< current volume [m3]
		double mass_flow_rate_in_0_;			//!< initial inlet mass flow rate [kg/s]
		double mass_flow_rate_in_;				//!< current inlet mass flow rate [kg/s]
		double mass_flow_rate_out_;				//!< current outlet mass flow rate [kg/s]
		double mass_flow_rate_loss_surface_;	//!< current mass flow rate loss from the surface [kg/s]
		double tau0_;							//!< initial residence time [s]
		double tau_;							//!< current residence time [s]

		double global_thermal_exchange_coefficient_;	//!< global thermal exchange coefficient [W/m2/K]
		double exchange_area_;							//!< exchange area [m2]
		double T_environment_;							//!< environment temperature

		double internal_area_;							//!< internal area [m2]

		std::vector<bool> site_non_conservation_;		//!< non-conservation of sites is allowed (true)

		unsigned int NC_;		//!< number of species [-]
		unsigned int NR_;		//!< number of reactions [-]
		unsigned int NE_;		//!< number of equations [-]

		unsigned int SURF_NC_;	//!< number of surface species [-]
		unsigned int SURF_NR_;	//!< number of surface reactions [-]
		unsigned int SURF_NP_;	//!< number of surface phases [-]

		unsigned int BULK_NC_;	//!< number of surface species [-]
		unsigned int BULK_NP_;	//!< number of surface phases [-]

		OpenSMOKE::OpenSMOKEVectorDouble omega_;			//!< current mass fractions
		OpenSMOKE::OpenSMOKEVectorDouble x_;				//!< current mole fractions
		OpenSMOKE::OpenSMOKEVectorDouble c_;				//!< current concentrations [kmol/m3]

		OpenSMOKE::OpenSMOKEVectorDouble Z_;				//!< current surface fractions
		OpenSMOKE::OpenSMOKEVectorDouble Gamma_;			//!< current surface densities [kmol/m2]
		OpenSMOKE::OpenSMOKEVectorDouble GammaFromEqn_;		//!< surface densities accounting for site non conservation [kmol/m2]
		OpenSMOKE::OpenSMOKEVectorDouble a_;				//!< current bulk species activities [-] (TODO: they are always assumed to be equal to 1)

		OpenSMOKE::OpenSMOKEVectorDouble RfromGas_;			//!< current formation rates (from homogeneous reactions)  [kmol/m3/s]
		OpenSMOKE::OpenSMOKEVectorDouble RfromSurface_;		//!< current formation rates (from heterogeneous reations) [kmol/m3/s]

		OpenSMOKE::OpenSMOKEVectorDouble Rsurface_;			//!< current surface species formation rates [kmol/m2/s]
		OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases_;	//!< current surface site phases formation rates [kmol/m2/s]

		OpenSMOKE::OpenSMOKEVectorDouble Rbulk_;			//!< current bulk species formation rates [kmol/m2/s]
		OpenSMOKE::OpenSMOKEVectorDouble massBulk_;			//!< current mass of bulk species [kg]

		double T_;									//!< current temperature [K]
		double P_;									//!< current pressure [Pa]
		double cTot_;								//!< current concentration [kmol/m3]
		double MW_;									//!< current molecular weight [kg/kmol]
		double QRGas_;								//!< current heat release from gas [W/m3]
		double QRSurface_;							//!< current heat release from surface [W/m2]
		double CvMixMass_;							//!< current specific heat (constant volume) [J/kg/K]
		double CpMixMass_;							//!< current specific heat (constant pressure) [J/kg/K]
		double H_;									//!< current enthalpy [J/kmol]
		double U_;									//!< current internal energy [J/kmol]

		unsigned int iteration_;					//!< iteration index
		unsigned int counter_file_video_;			//!< iteration counter for writing on video
		unsigned int counter_file_ASCII_;			//!< iteration counter for writing on ASCII file
		unsigned int counter_file_XML_;				//!< iteration counter for writing on XML file
		unsigned int counter_sensitivity_XML_;		//!< iteration counter for writing on XML sensitivity files

		OpenSMOKE::OpenSMOKEVectorDouble y0_;		//!< vector containing the initial values for all the variables
		OpenSMOKE::OpenSMOKEVectorDouble yf_;		//!< vector containing the final values for all the variables

		std::ofstream fXML_;						//!< XML file where the output is written
		std::ofstream fASCIITime_;					//!< ASCII file where the output is written during time integration
		std::ofstream fBulkASCIITime_;				//!< ASCII file where the output is written during time integration (bulk species)
		std::ofstream fASCIIFinal_;					//!< ASCII file where the final output is written
		std::ofstream fBulkASCIIFinal_;				//!< ASCII file where the final output is written (bulk species)
		std::ofstream fROPA_;						//!< ASCII file where the ROPA (on the fly) is written
		std::ofstream  fSensitivityParentXML_;		//!< XML file where the sensitivity analysis is written (parent file)
		std::ofstream* fSensitivityChildXML_;		//!< XML files where the sensitivity analysis is written (children files)

		std::vector<unsigned int> indices_of_output_species_;			//!< indices of species for which the profiles are written on the ASCII file
		std::vector<unsigned int> indices_of_sensitivity_species_;		//!< indices of species for which the sensitivity analysis is written on file
		std::vector<unsigned int> widths_of_output_species_;			//!< width of species for which the profiles are written on the ASCII file
		std::vector<unsigned int> widths_of_output_surface_species_;	//!< width of surface species for which the profiles are written on the ASCII file

		OpenSMOKE::SensitivityMap* sensitivityMap_;			//!< sensitivity map
		OpenSMOKE::OpenSMOKEVectorDouble scaling_Jp;		//!< vector containing datafor performing the sensitivity analysis
		OpenSMOKE::OpenSMOKEMatrixDouble J;					//!< Jacobian matrix (required for sensitivity analysis)
	};
}

#include "SurfacePerfectlyStirredReactor.hpp"

#endif	/* OpenSMOKE_SurfacePerfectlyStirredReactor_H */

