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

#ifndef OpenSMOKE_SurfaceBatchReactor_H
#define	OpenSMOKE_SurfaceBatchReactor_H

#include "math/OpenSMOKEVector.h"
#include "SurfaceBatchReactor_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysis_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"
#include "boost/filesystem.hpp"
#include "math/external-ode-solvers/ODE_Parameters.h"

namespace OpenSMOKE
{
	enum SurfaceBatchReactor_Type { SURFACEBATCH_REACTOR_NONISOTHERMAL_CONSTANTP, SURFACEBATCH_REACTOR_NONISOTHERMAL_CONSTANTV,
									SURFACEBATCH_REACTOR_ISOTHERMAL_CONSTANTP, SURFACEBATCH_REACTOR_ISOTHERMAL_CONSTANTV };

	//!  A class for simulating batch reactors
	/*!
		 The purpose of this class is to simulate batch reactors. It provides a common interface to different
		 batch reactors
	*/

	class SurfaceBatchReactor
	{

	public:

		/**
		*@brief Solves the batch reactor
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
		*@param thermodynamicsMap map of thermodynamic data
		*@param sensitivity_options options for sensitivity analysis
		*/		
		void EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap& sensitivityMap, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfaceBatchReactor_Options& batch_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options);
	
		/**
		*@brief Returns the current value of temperature [K]
		*/	
		double T() const { return T_; }

		/**
		*@brief Returns the current value of pressure [Pa]
		*/	
		double P_Pa() const { return P_; }

		/**
		*@brief Returns the current value of volume [m3]
		*/	
		double V() const { return V_; }

		/**
		*@brief Returns the current values of mass fractions [-]
		*/	
		const OpenSMOKE::OpenSMOKEVectorDouble& omega() const { return omega_; }
	
	protected:

		/**
		*@brief Prepares the ASCII file
		*@param output_file_ascii ASCII file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param batch_options options for solving the batch reactor
		*/
		void PrepareASCIIFile(const boost::filesystem::path output_file_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfaceBatchReactor_Options& batch_options);
	
		/**
		*@brief Prepares the ASCII file for bulk species
		*@param output_file_bulk_ascii ASCII file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param batch_options options for solving the batch reactor
		*/
		void PrepareBulkASCIIFile(const boost::filesystem::path output_file_bulk_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfaceBatchReactor_Options& batch_options);

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
		*@brief Allocates the memory for the vectors and the matrices contained in the batch reactor
		*/			
		void MemoryAllocation();

		/**
		*@brief Prepares the XML files for the sensitivity analysis
		*@param thermodynamicsMap map of thermodynamic data
		*@param batch_options options for solving the batch reactor
		*/
		void PrepareSensitivityXMLFiles(OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfaceBatchReactor_Options& batch_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options);
		
		/**
		*@brief Closes the XML files for the sensitivity analysis
		*/			
		void CloseSensitivityXMLFiles();

		/**
		*@brief Open all the files (ASCII and XML)
		*@param thermodynamicsMap map of thermodynamic data
		*@param batch_options options for solving the batch reactor
		*/
		void OpenAllFiles(OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, OpenSMOKE::SurfaceBatchReactor_Options& batch_options);
		
		/**
		*@brief Closes all the files (ASCII and XML)
		*@param batch_options options for solving the batch reactor
		*/		
		void CloseAllFiles(OpenSMOKE::SurfaceBatchReactor_Options& batch_options);

		/**
		*@brief Solves the Ordinary Differential Equations using one of the provided open/source ODE solvers
		*@param tf final time of integrations
		*@param ode_parameters parameters governing the ODE integration
		*/	
		void SolveOpenSourceSolvers(const double tf, OpenSMOKE::ODE_Parameters& ode_parameters);
		
		/**
		*@brief Writes on the video the final status of the batch reactor
		*@param tf final time of integration
		*@param thermodynamicsMap map of thermodynamic data
		*@param thermodynamicsSurfaceMap map of surface thermodynamic data
		*/			
		void FinalStatus(const double tf, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap);

		/**
		*@brief Writes on a file the final summary of the batch reactor
		*@param summary_file file where the summary will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param thermodynamicsSurfaceMap map of surface thermodynamic data
		*/		
		void FinalSummary(const boost::filesystem::path summary_file, const double tf, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap);

		/**
		*@brief Returns the final composition of the batch reactor, together with temperature and pressure
		*@param T temperature [K]
		*@param P pressure [Pa]
		*@param omega mass fractions [-]
		*@param surface fractions [-]
		*@param site densities fractions [kmol/m2]
		*/	
		void GetFinalStatus(double& T, double& P, OpenSMOKE::OpenSMOKEVectorDouble& omega, OpenSMOKE::OpenSMOKEVectorDouble& Z, OpenSMOKE::OpenSMOKEVectorDouble& Gamma, double& mass, OpenSMOKE::OpenSMOKEVectorDouble& massBulk);

	
protected:

		SurfaceBatchReactor_Type type_;				//!< type of batch reactor

		double T0_;									//!< initial temperature [K]
		double P0_;									//!< initial pressure [Pa]
		double mass0_;								//!< initial mass [kg]
		double massBulk0_;							//!< initial mass of bulk phase [kg] (always equal to 0, to be considered in a relative way)
		double V0_;									//!< initial volume [m3]
		double A0_;									//!< initial internal surface area [m2]
		double MW0_;								//!< initial molecular weight [kg/kmol]
		double H0_;									//!< initial enthalpy [J/kmol]
		double H0_Surface_;							//!< initial surface enthalpy [J]
		double U0_;									//!< initial internal energy [J/kmol]
		double U0_Surface_;							//!< initial surface internal energy [J]
		OpenSMOKE::OpenSMOKEVectorDouble omega0_;	//!< initial composition [mass fractions]
		OpenSMOKE::OpenSMOKEVectorDouble Z0_;		//!< initial composition [surface fractions]
		OpenSMOKE::OpenSMOKEVectorDouble Gamma0_;	//!< initial surface densities [kmol/m2]
		OpenSMOKE::OpenSMOKEVectorDouble x0_;		//!< initial composition [mole fractions]

		double rho_;	//!< density [kg/m3]
		double mass_;	//!< mass [kg]

		std::vector<bool> site_non_conservation_;	//!< non conservation of sites is allowed (if true)

		unsigned int NC_;		//!< number of species [-]
		unsigned int NR_;		//!< number of reactions [-]

		unsigned int SURF_NC_;	//!< number of surface species [-]
		unsigned int SURF_NR_;	//!< number of surface reactions [-]
		unsigned int SURF_NP_;	//!< number of surface phases [-]

		unsigned int BULK_NC_;	//!< number of surface species [-]
		unsigned int BULK_NP_;	//!< number of surface phases [-]

		unsigned int NE_;		//!< number of equations [-]

		OpenSMOKE::OpenSMOKEVectorDouble omega_;			//!< current mass fractions
		OpenSMOKE::OpenSMOKEVectorDouble Z_;				//!< current surface fractions
		OpenSMOKE::OpenSMOKEVectorDouble Gamma_;			//!< current surface densities [kmol/m2]
		OpenSMOKE::OpenSMOKEVectorDouble GammaFromEqn_;		//!< surface densities from equation [kmol/m2]
		OpenSMOKE::OpenSMOKEVectorDouble a_;				//!< current bulk species activities [-] (TODO: they are always assumed to be equal to 1)
		OpenSMOKE::OpenSMOKEVectorDouble massBulk_;			//!< current mass of bulk species [kg]
		OpenSMOKE::OpenSMOKEVectorDouble x_;				//!< current mole fractions
		OpenSMOKE::OpenSMOKEVectorDouble c_;				//!< current concentrations [kmol/m3]
		OpenSMOKE::OpenSMOKEVectorDouble RfromGas_;			//!< current formation rates (from homogeneous reactions)  [kmol/m3/s]
		OpenSMOKE::OpenSMOKEVectorDouble RfromSurface_;		//!< current formation rates (from heterogeneous reations) [kmol/m2/s]

		OpenSMOKE::OpenSMOKEVectorDouble Rsurface_;			//!< current surface species formation rates [kmol/m2/s]
		OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases_;	//!< current surface site phases formation rates [kmol/m2/s]

		OpenSMOKE::OpenSMOKEVectorDouble Rbulk_;			//!< current bulk species formation rates [kmol/m2/s]

		double T_;									//!< current temperature [K]
		double P_;									//!< current pressure [Pa]
		double V_;									//!< current volume [m3]
		double A_;									//!< current internal surface area [m2]
		double cTot_;								//!< current concentration [kmol/m3]
		double MW_;									//!< current molecular weight [kg/kmol]
		double QRGas_;								//!< current heat release from gas phase[W/m3]
		double QRSurface_;							//!< current heat release from surface [W/m2]
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
		std::ofstream fBulkASCII_;					//!< ASCII file where the output for bulk species is written
		std::ofstream  fSensitivityParentXML_;		//!< XML file where the sensitivity analysis is written (parent file)
		std::ofstream* fSensitivityChildXML_;		//!< XML files where the sensitivity analysis is written (children files)

		std::vector<unsigned int> indices_of_output_species_;			//!< indices of species for which the profiles are written on the ASCII file
		std::vector<unsigned int> indices_of_sensitivity_species_;		//!< indices of species for which the sensitivity analysis is written on file
		std::vector<unsigned int> widths_of_output_species_;			//!< width   of species for which the profiles are written on the ASCII file
		std::vector<unsigned int> widths_of_output_surface_species_;	//!< width   of surface species for which the profiles are written on the ASCII file
		
		OpenSMOKE::SensitivityMap* sensitivityMap_;			//!< sensitivity map
		OpenSMOKE::OpenSMOKEVectorDouble scaling_Jp;		//!< vector containing datafor performing the sensitivity analysis
		OpenSMOKE::OpenSMOKEMatrixDouble J;					//!< Jacobian matrix (required for sensitivity analysis)
	};
}

#include "SurfaceBatchReactor.hpp"

#endif	/* OpenSMOKE_SurfaceBatchReactor_H */

