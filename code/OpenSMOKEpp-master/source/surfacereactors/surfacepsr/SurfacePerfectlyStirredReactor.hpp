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

#include "math/OpenSMOKEVector.h"
#include "SurfacePerfectlyStirredReactor_OdeInterfaces.h"
#include "SurfacePerfectlyStirredReactor_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"
#include "math/external-ode-solvers/ODE_Parameters.h"

namespace OpenSMOKE
{
	SurfacePerfectlyStirredReactor::SurfacePerfectlyStirredReactor(
		OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
		OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
		OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap,
		OpenSMOKE::KineticsMap_Surface_CHEMKIN& kineticsSurfaceMap,
		OpenSMOKE::ODE_Parameters& ode_parameters,
		OpenSMOKE::SurfacePerfectlyStirredReactor_Options& psr_options,
		OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa) :

		thermodynamicsMap_(thermodynamicsMap),
		kineticsMap_(kineticsMap),
		thermodynamicsSurfaceMap_(thermodynamicsSurfaceMap),
		kineticsSurfaceMap_(kineticsSurfaceMap),
		ode_parameters_(ode_parameters),
		psr_options_(psr_options),
		on_the_fly_ropa_(on_the_fly_ropa)

	{
	}

	void SurfacePerfectlyStirredReactor::EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap& sensitivityMap, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options)
	{
		sensitivityMap_ = &sensitivityMap;

		PrepareSensitivityXMLFiles(sensitivity_options);
		
		ChangeDimensions(NE_, &scaling_Jp, true);
		ChangeDimensions(NE_, NE_, &J, true);
	}

	void SurfacePerfectlyStirredReactor::MemoryAllocation()
	{
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &RfromGas_, true);
		ChangeDimensions(NC_, &RfromSurface_, true);

		ChangeDimensions(NE_, &y0_, true);
		ChangeDimensions(NE_, &yf_, true);

		ChangeDimensions(NC_, &x0_, true);
		ChangeDimensions(NC_, &xInlet_, true);

		ChangeDimensions(SURF_NC_, &Z_, true);
		ChangeDimensions(SURF_NP_, &Gamma_, true);
		ChangeDimensions(SURF_NC_, &Rsurface_, true);
		ChangeDimensions(SURF_NP_, &RsurfacePhases_, true);

		ChangeDimensions(BULK_NC_, &massBulk_, true);
		ChangeDimensions(BULK_NC_, &a_, true);
		ChangeDimensions(BULK_NC_, &Rbulk_, true);

		// TODO: in the current implementation activities are always equal to 1
		a_ = 1.;
	}

	void SurfacePerfectlyStirredReactor::NumericalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEMatrixDouble& J)
	{
		// Calculated as suggested by Buzzi (private communication)

		const double ZERO_DER = sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);			
		const double TOLR = 100. * OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT;
		const double TOLA = 1.e-10;
	
		OpenSMOKE::OpenSMOKEVectorDouble y_plus = y;
		OpenSMOKE::OpenSMOKEVectorDouble dy_original(y.Size());
		OpenSMOKE::OpenSMOKEVectorDouble dy_plus(y.Size());

		Equations(t, y, dy_original);

		// Derivatives with respect to y[kd]
		for(int kd=1;kd<=y.Size();kd++)
		{
			double hf = 1.e0;
			double error_weight = 1./(TOLA+TOLR*std::fabs(y[kd]));
			double hJ = ETA2 * std::fabs(std::max(y[kd], 1./error_weight));
			double hJf = hf/error_weight;
			hJ = std::max(hJ, hJf);
			hJ = std::max(hJ, ZERO_DER);

			// This is what is done by Buzzi
			double dy = std::min(hJ, 1.e-3 + 1e-3*std::fabs(y[kd]));
			double udy = 1. / dy;
			y_plus[kd] += dy;
			Equations(t, y_plus, dy_plus);

			for(int j=1;j<=y.Size();j++)
				J[j][kd] = (dy_plus[j]-dy_original[j]) * udy;

			y_plus[kd] = y[kd];
		}
	}

	void SurfacePerfectlyStirredReactor::PrepareROPAFile(const boost::filesystem::path output_file_ropa)
	{
		// Open ROPA file
		fROPA_.open(output_file_ropa.c_str(), std::ios::out);

		// Write head in file
		on_the_fly_ropa_.WriteHead(fROPA_, "Perfectly Stirred Reactor");
	}

	void SurfacePerfectlyStirredReactor::PrepareXMLFile(const boost::filesystem::path output_file_xml)
	{
		fXML_.open(output_file_xml.c_str(), std::ios::out);
		OpenSMOKE::SetXMLFile(fXML_);
		fXML_ << "<Type> SurfacePerfectlyStirredReactor </Type>" << std::endl;

		unsigned int counter = 2;
		fXML_ << "<additional>" << std::endl;
		fXML_ << 6 << std::endl;
		fXML_ << "time [s] " << counter++ << std::endl;
		fXML_ << "temperature [K] " << counter++ << std::endl;
		fXML_ << "pressure [Pa] " << counter++ << std::endl;
		fXML_ << "mol-weight [kg/kmol] " << counter++ << std::endl;
		fXML_ << "density [kg/m3] " << counter++ << std::endl;
		fXML_ << "heat-release [W/m3] " << counter++ << std::endl;
		fXML_ << "</additional>" << std::endl;

		fXML_ << "<t-p-mw>" << std::endl;
		fXML_ << 1 << " " << 2 << " " << 3 << std::endl;
		fXML_ << "</t-p-mw>" << std::endl;

		fXML_ << "<mass-fractions>" << std::endl;
		fXML_ << thermodynamicsMap_.NumberOfSpecies() << std::endl;
		for (unsigned int j=0;j<NC_;j++)
			fXML_ << thermodynamicsMap_.NamesOfSpecies()[j] << " " << thermodynamicsMap_.MW(j) << " " << counter++ << std::endl;
		fXML_ << "</mass-fractions>" << std::endl;
		fXML_ << "<profiles>" << std::endl;
	}

	void SurfacePerfectlyStirredReactor::CloseXMLFile()
	{
		fXML_ << "</profiles>" << std::endl;
		fXML_ << "<profiles-size> " << std::endl;
		fXML_ << counter_file_XML_ << " " << 1 + (NC_+1) << std::endl;
		fXML_ << "</profiles-size> " << std::endl;
		fXML_ << "</opensmoke>" << std::endl;
	}

	void SurfacePerfectlyStirredReactor::PrepareSensitivityXMLFiles(OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options)
	{
		indices_of_sensitivity_species_.resize(sensitivity_options.list_of_species().size());
		for(unsigned int i=0;i<indices_of_sensitivity_species_.size();i++)
			indices_of_sensitivity_species_[i] = thermodynamicsMap_.IndexOfSpecies(sensitivity_options.list_of_species()[i]);

		const boost::filesystem::path parent_file = psr_options_.output_path() / "Sensitivities.xml";
		fSensitivityParentXML_.open(parent_file.c_str(), std::ios::out);
		SetXMLFile(fSensitivityParentXML_);
		fSensitivityParentXML_ << "<variables>" << std::endl;
		
		if (type_ == SURFACEPERFECTLYSTIRRED_REACTOR_NONISOTHERMAL_CONSTANTP)
			fSensitivityParentXML_ << indices_of_sensitivity_species_.size()+1 << std::endl;
		else if (type_ == SURFACEPERFECTLYSTIRRED_REACTOR_ISOTHERMAL_CONSTANTP)
			fSensitivityParentXML_ << indices_of_sensitivity_species_.size() << std::endl;

		for (unsigned int j=0;j<indices_of_sensitivity_species_.size();j++)
			fSensitivityParentXML_ << thermodynamicsMap_.NamesOfSpecies()[indices_of_sensitivity_species_[j]-1] << " " << j << " " << indices_of_sensitivity_species_[j] << std::endl;
		
		if (type_ == SURFACEPERFECTLYSTIRRED_REACTOR_NONISOTHERMAL_CONSTANTP )
			fSensitivityParentXML_ << "temperature" << " " << indices_of_sensitivity_species_.size() << " " << NC_+1 << std::endl;
		
		fSensitivityParentXML_ << "</variables>" << std::endl;
		fSensitivityParentXML_ << "<n-parameters> " << std::endl;
		fSensitivityParentXML_ << sensitivityMap_->number_of_parameters() << std::endl;
		fSensitivityParentXML_ << "</n-parameters> " << std::endl;

		fSensitivityChildXML_ = new std::ofstream[indices_of_sensitivity_species_.size()+1];
		for (unsigned int j=0;j<indices_of_sensitivity_species_.size();j++)
		{
			const std::string name = "Sensitivities." + thermodynamicsMap_.NamesOfSpecies()[indices_of_sensitivity_species_[j]-1] +".xml";
			const boost::filesystem::path child_file = psr_options_.output_path() / name;
			fSensitivityChildXML_[j].open(child_file.c_str(), std::ios::out);
		}
		{
			const boost::filesystem::path child_file = psr_options_.output_path() / "Sensitivities.temperature.xml";
			fSensitivityChildXML_[indices_of_sensitivity_species_.size()].open(child_file.c_str(), std::ios::out);
		}

		for (unsigned int j=0;j<indices_of_sensitivity_species_.size()+1;j++)
		{
			SetXMLFile(fSensitivityChildXML_[j]);
			fSensitivityChildXML_[j] << std::setprecision(5);
			fSensitivityChildXML_[j] << "<coefficients>" << std::endl;
		}		
	}

	void SurfacePerfectlyStirredReactor::CloseSensitivityXMLFiles()
	{
		fSensitivityParentXML_ << "<points> " << std::endl;
		fSensitivityParentXML_ << counter_sensitivity_XML_ << std::endl;
		fSensitivityParentXML_ << "</points> " << std::endl;
		fSensitivityParentXML_ << "<constant-parameters> " << std::endl;
		for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
			fSensitivityParentXML_ << sensitivityMap_->parameters()[j] << std::endl;
		fSensitivityParentXML_ << "</constant-parameters> " << std::endl;
		fSensitivityParentXML_ << "</opensmoke>" << std::endl;

		for (unsigned int j=0;j<indices_of_sensitivity_species_.size()+1;j++)
		{
			fSensitivityChildXML_[j] << "</coefficients>" << std::endl;
			fSensitivityChildXML_[j] << "</opensmoke>" << std::endl;
		}
	}

	void SurfacePerfectlyStirredReactor::PrepareASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii)
	{
		indices_of_output_species_.resize(psr_options_.output_species().size());
		for (unsigned int i = 0; i<psr_options_.output_species().size(); i++)
			indices_of_output_species_[i] = thermodynamicsMap_.IndexOfSpecies(psr_options_.output_species()[i]);

		if (indices_of_output_species_.size() != 0)
		{
			widths_of_output_species_.resize(psr_options_.output_species().size());
			for (unsigned int i = 0; i<psr_options_.output_species().size(); i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(psr_options_.output_species()[i], NC_);
		}
		else
		{
			widths_of_output_species_.resize(NC_);
			for (unsigned int i = 0; i<NC_; i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap_.NamesOfSpecies()[i], NC_);
		}

		{
			widths_of_output_surface_species_.resize(SURF_NC_);
			for (unsigned int i = 0; i<SURF_NC_; i++)
				widths_of_output_surface_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsSurfaceMap_.NamesOfSpecies()[NC_ + i], NC_);
		}

		fOutput.open(output_file_ascii.c_str(), std::ios::out);

		unsigned int counter = 1;
		fOutput.setf(std::ios::scientific);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "t[s]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T0[K]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P0[Pa]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "V0[m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P[Pa]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "V[m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "rho[kg/m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "MW[kg/kmol]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Qgas[W/m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Qsurf[W/m2]", counter);

		if (BULK_NC_ != 0)
		{
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "bulk[kg]", counter);

			for (unsigned int i = 0; i < BULK_NC_; i++)
			{
				const std::string name = thermodynamicsSurfaceMap_.vector_names_bulk_species()[i].substr(0, 11) + "[kg]";
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, name, counter);
			}
		}

		for (unsigned int i = 0; i< thermodynamicsSurfaceMap_.number_of_site_phases(0); i++)
		{
			const std::string name = thermodynamicsSurfaceMap_.matrix_names_site_phases()[0][i].substr(0, 12) + "[kmol/m2]";
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, name, counter);
		}

		for (unsigned int i = 0; i<SURF_NC_; i++)
			OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_surface_species_[i], fOutput, thermodynamicsSurfaceMap_.NamesOfSpecies()[NC_ + i] + "_z", counter);


		if (indices_of_output_species_.size() != 0)
		{
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap_.NamesOfSpecies()[indices_of_output_species_[i] - 1] + "_x", counter);
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap_.NamesOfSpecies()[indices_of_output_species_[i] - 1] + "_w", counter);
		}
		else
		{
			for (unsigned int i = 0; i<NC_; i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap_.NamesOfSpecies()[i] + "_x", counter);
			for (unsigned int i = 0; i<NC_; i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap_.NamesOfSpecies()[i] + "_w", counter);
		}

		fOutput << std::endl;
	}

	void SurfacePerfectlyStirredReactor::PrepareBulkASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_bulk_ascii)
	{
		fOutput.open(output_file_bulk_ascii.c_str(), std::ios::out);

		unsigned int counter = 1;
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "t[s]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P[Pa]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "V[m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "rho[kg/m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "MW[kg/kmol]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Qgas[W/m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Qsurf[W/m2]", counter);

		OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, "total[kg]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, "total[kg/m2/s]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, "total[kg/s]", counter);

		for (unsigned int i = 0; i < BULK_NC_; i++)
		{
			const std::string name_species = thermodynamicsSurfaceMap_.vector_names_bulk_species()[i].substr(0, 9);

			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[kg]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[kg/m2/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[kg/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[m]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[micron]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[m/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[mm/hr]", counter);
		}

		fOutput << std::endl;
	}

	void SurfacePerfectlyStirredReactor::SolveOpenSourceSolvers(const double tf)
	{
		#if OPENSMOKE_USE_DVODE == 1
		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_DVODE)
		{
			// ODE system
			ODESystem_DVODE_SurfacePerfectlyStirredReactor *odeSystemObject;
			odeSystemObject = ODESystem_DVODE_SurfacePerfectlyStirredReactor::GetInstance();
			odeSystemObject->SetSurfacePerfectlyStirredReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DVODE<ODESystem_DVODE_SurfacePerfectlyStirredReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters_.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters_.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters_.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters_.SetCPUTime(tEnd-tStart);
			ode_parameters_.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters_.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters_.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters_.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters_.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters_.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters_.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters_.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters_.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif

		#if OPENSMOKE_USE_ODEPACK == 1
		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_DLSODE)
		{
			// ODE system
			ODESystem_DLSODE_SurfacePerfectlyStirredReactor *odeSystemObject;
			odeSystemObject = ODESystem_DLSODE_SurfacePerfectlyStirredReactor::GetInstance();
			odeSystemObject->SetSurfacePerfectlyStirredReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DLSODE<ODESystem_DLSODE_SurfacePerfectlyStirredReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters_.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters_.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters_.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters_.SetCPUTime(tEnd-tStart);
			ode_parameters_.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters_.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters_.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters_.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters_.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters_.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters_.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters_.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters_.SetLastStepUsed(o.GetLastStepUsed());
		}

		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_DLSODA)
		{
			// ODE system
			ODESystem_DLSODA_SurfacePerfectlyStirredReactor *odeSystemObject;
			odeSystemObject = ODESystem_DLSODA_SurfacePerfectlyStirredReactor::GetInstance();
			odeSystemObject->SetSurfacePerfectlyStirredReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DLSODA<ODESystem_DLSODA_SurfacePerfectlyStirredReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters_.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters_.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters_.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters_.SetCPUTime(tEnd-tStart);
			ode_parameters_.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters_.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters_.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters_.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters_.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters_.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters_.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters_.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters_.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif

		#if OPENSMOKE_USE_DASPK == 1
		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_DASPK)
		{
			// ODE system
			ODESystem_DASPK_SurfacePerfectlyStirredReactor *odeSystemObject;
			odeSystemObject = ODESystem_DASPK_SurfacePerfectlyStirredReactor::GetInstance(NE_);
			odeSystemObject->SetSurfacePerfectlyStirredReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DASPK<ODESystem_DASPK_SurfacePerfectlyStirredReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters_.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters_.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters_.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters_.SetCPUTime(tEnd-tStart);
			ode_parameters_.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters_.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters_.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters_.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters_.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters_.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters_.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters_.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters_.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif

		#if OPENSMOKE_USE_RADAU == 1
		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_RADAU5)
		{
			// ODE system
			ODESystem_RADAU5_SurfacePerfectlyStirredReactor *odeSystemObject;
			odeSystemObject = ODESystem_RADAU5_SurfacePerfectlyStirredReactor::GetInstance();
			odeSystemObject->SetSurfacePerfectlyStirredReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_RADAU<ODESystem_RADAU5_SurfacePerfectlyStirredReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters_.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters_.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters_.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters_.SetCPUTime(tEnd-tStart);
			ode_parameters_.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters_.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters_.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters_.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters_.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters_.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters_.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters_.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters_.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif

		#if OPENSMOKE_USE_SUNDIALS == 1
		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_CVODE)
		{
			// ODE system
			ODESystem_CVODE_SurfacePerfectlyStirredReactor *odeSystemObject;
			odeSystemObject = ODESystem_CVODE_SurfacePerfectlyStirredReactor::GetInstance();
			odeSystemObject->SetSurfacePerfectlyStirredReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_CVODE_Sundials<ODESystem_CVODE_SurfacePerfectlyStirredReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters_.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters_.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters_.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters_.SetCPUTime(tEnd-tStart);
			ode_parameters_.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters_.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters_.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters_.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters_.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters_.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters_.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters_.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters_.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif

		#if OPENSMOKE_USE_MEBDF == 1
		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_MEBDF)
		{
			// ODE system
			ODESystem_MEBDF_SurfacePerfectlyStirredReactor *odeSystemObject;
			odeSystemObject = ODESystem_MEBDF_SurfacePerfectlyStirredReactor::GetInstance();
			odeSystemObject->SetSurfacePerfectlyStirredReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_MEBDF<ODESystem_MEBDF_SurfacePerfectlyStirredReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters_.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters_.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters_.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters_.SetCPUTime(tEnd-tStart);
			ode_parameters_.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters_.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters_.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters_.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters_.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters_.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters_.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters_.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters_.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif
	}

	void SurfacePerfectlyStirredReactor::OpenAllFiles()
	{
		if (psr_options_.verbose_output() == true)
		{
			if ( !boost::filesystem::exists(psr_options_.output_path()) )
				OpenSMOKE::CreateDirectory(psr_options_.output_path());

			if (psr_options_.verbose_ascii_file() == true)
			{
				PrepareASCIIFile(fASCIIFinal_, psr_options_.output_path() / "Output.out");
				PrepareASCIIFile(fASCIITime_, psr_options_.output_path() / "Output.history");

				PrepareBulkASCIIFile(fBulkASCIIFinal_, psr_options_.output_path() / "Output.bulk.out");
				PrepareBulkASCIIFile(fBulkASCIITime_, psr_options_.output_path() / "Output.bulk.history");
			}

			if (psr_options_.verbose_xml_file() == true)
				PrepareXMLFile(psr_options_.output_path() / "Output.xml");

			if (on_the_fly_ropa_.is_active() == true)
				PrepareROPAFile(psr_options_.output_path() / "ROPA.out");
		}
	}

	void SurfacePerfectlyStirredReactor::CloseAllFiles()
	{
		if (psr_options_.verbose_output() == true)
			{
				if (psr_options_.verbose_ascii_file() == true)
				{
					fASCIIFinal_.close();
					fASCIITime_.close();
				}

				if (psr_options_.verbose_xml_file() == true)
					CloseXMLFile();

				if (on_the_fly_ropa_.is_active() == true)
					fROPA_.close();
			}

		if (psr_options_.sensitivity_analysis() == true)
			CloseSensitivityXMLFiles();
	}

	void SurfacePerfectlyStirredReactor::FinalStatus(const double tf)
	{
		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_);
		Equations(tf, yf_, dummy);

		if (psr_options_.verbose_video() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << std::setw(30) << std::left << "Short Summary" << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << std::setw(30) << std::left << "Residence time:" << std::setw(12) << std::left << tau_ << "s" << std::endl;
			std::cout << std::setw(30) << std::left << "Mass:" << std::setw(12) << std::left << rho_*V_*1000. << "g" << std::endl;
			std::cout << std::setw(30) << std::left << "Mass flow rate (in):" << std::setw(12) << std::left << mass_flow_rate_in_*1e3 << "g/s" << std::endl;
			std::cout << std::setw(30) << std::left << "Mass flow rate (out):" << std::setw(12) << std::left << mass_flow_rate_out_*1e3 << "g/s" << std::endl;
			std::cout << std::setw(30) << std::left << "Mass flow rate (loss):" << std::setw(12) << std::left << mass_flow_rate_loss_surface_*1e3 << "g/s" << std::endl;
			std::cout << std::setw(30) << std::left << "Volume:" << std::setw(12) << std::left << V_*1.e6 << "cm3" << std::endl;
			std::cout << std::setw(30) << std::left << "Density:" << std::setw(12) << std::left << rho_ / 1000. << "g/cm3" << std::endl;
			std::cout << std::setw(30) << std::left << "Temperature:" << std::setw(12) << std::left << T_ << "K" << std::endl;
			std::cout << std::setw(30) << std::left << "Pressure:" << std::setw(12) << std::left << P_ / 101325. << "atm" << std::endl;
			std::cout << std::setw(30) << std::left << "Integration time:" << std::setw(12) << std::left << tf << "s" << std::endl;
			std::cout << std::setw(30) << std::left << "Max derivative:" << std::setw(12) << std::left << dummy.MaxAbs() << "" << std::endl;


			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << std::setw(30) << std::left << "Status" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << std::setw(30) << std::left << "Gamma[kmol/m2]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Gamma0_.SumElements() << Gamma_.SumElements() << std::endl;
			std::cout << std::setw(30) << std::left << "Bulk[kg]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << massBulk0_.SumElements() << massBulk_.SumElements() << std::endl;

			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << std::setw(30) << std::left << "Status" << std::setw(20) << std::left << "Inlet" << "Outlet" << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << std::setw(30) << std::left << "Time[s]" << std::setw(20) << std::left << 0. << tau_ << std::endl;
			std::cout << std::setw(30) << std::left << "T[K]" << std::setw(20) << std::left << std::fixed << std::setprecision(3) << TInlet_ << T_ << std::endl;
			std::cout << std::setw(30) << std::left << "P[atm]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << PInlet_ / 101325. << P_ / 101325. << std::endl;
			std::cout << std::setw(30) << std::left << "rho[kg/m3]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << rhoInlet_ << rho_ << std::endl;
			std::cout << std::setw(30) << std::left << "MW[kg/kmol]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << MWInlet_ << MW_ << std::endl;

			// Equivalence ratio
			{
				const double phiInlet = thermodynamicsMap_.GetLocalEquivalenceRatioFromMoleFractions(x0_.GetHandle());
				const double phi = thermodynamicsMap_.GetLocalEquivalenceRatioFromMoleFractions(x_.GetHandle());
				std::cout << std::setw(30) << std::left << "Phi[-]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << phiInlet << phi << std::endl;
			}

			// Enthalpy analysis
			{
				const double H_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x_.GetHandle());
				const double U_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(x_.GetHandle());

				std::cout << std::setw(30) << std::left << "H[J/kg]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << HInlet_ / MWInlet_ << H_ / MW_ << std::endl;
				std::cout << std::setw(30) << std::left << "U[J/kg]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << UInlet_ / MWInlet_ << U_ / MW_ << std::endl;
			}

			for (unsigned int i = 1; i <= NC_; i++)
			if (omegaInlet_[i] > 0.)
			{
				std::string label = "Conv.(%) " + thermodynamicsMap_.NamesOfSpecies()[i - 1];
				std::cout << std::setw(30) << std::left << label << std::setw(20) << std::left << std::fixed << std::setprecision(3) << 0. << (omegaInlet_[i] - omega_[i]) / omegaInlet_[i] * 100. << std::endl;
			}

			// Atomic analysis
			{
				std::vector<double> sum_inlet(thermodynamicsMap_.elements().size());
				std::vector<double> sum_final(thermodynamicsMap_.elements().size());
				for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
				{
					sum_inlet[j] = 0.;
					sum_final[j] = 0.;
					for (unsigned int i = 0; i < NC_; i++)
					{
						sum_inlet[j] += thermodynamicsMap_.atomic_composition()(i, j) * xInlet_[i + 1];
						sum_final[j] += thermodynamicsMap_.atomic_composition()(i, j) * x_[i + 1];
					}
				}

				const double moles_inlet = mass_flow_rate_in_ / MWInlet_;
				const double moles_final = mass_flow_rate_out_ / MW_;

				std::cout << std::setw(30) << std::left << "Total species (kmol/s) " << std::setw(20) << std::left << std::scientific << std::setprecision(6) << moles_inlet << moles_final << std::endl;

				for (unsigned int j = 0; j<thermodynamicsMap_.elements().size(); j++)
				if (sum_inlet[j] > 0.)
				{
					std::string label = "Element(kmol/s) " + thermodynamicsMap_.elements()[j];
					std::cout << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(6) << sum_inlet[j] * moles_inlet << sum_final[j] * moles_final << std::endl;
				}
			}

			std::cout << "-----------------------------------------------------------------------------" << std::endl;

			if (psr_options_.sensitivity_analysis() == true)
			{
				const double totalCumulative = sensitivityMap_->cpuTimeFactorization() + sensitivityMap_->cpuTimeAssembling() + sensitivityMap_->cpuTimeSolution();

				std::cout << std::endl;
				std::cout << "-----------------------------------------------------------------------------" << std::endl;
				std::cout << "                              Sensitivity Analysis                           " << std::endl;
				std::cout << "-----------------------------------------------------------------------------" << std::endl;

				std::cout << "Cumulative times" << std::endl;
				std::cout << " * Factorizing:    " << sensitivityMap_->cpuTimeFactorization() << " s"
					<< "(" << sensitivityMap_->cpuTimeFactorization() / totalCumulative*100. << "%)" << std::endl;
				std::cout << " * Assembling rhs: " << sensitivityMap_->cpuTimeAssembling() << " s"
					<< "(" << sensitivityMap_->cpuTimeAssembling() / totalCumulative*100. << "%)" << std::endl;
				std::cout << " * Solving:        " << sensitivityMap_->cpuTimeSolution() << " s"
					<< "(" << sensitivityMap_->cpuTimeSolution() / totalCumulative*100. << "%)" << std::endl;

				std::cout << "Single times" << std::endl;
				std::cout << " * Factorizing:    " << sensitivityMap_->cpuTimeSingleFactorization() << " s" << std::endl;
				std::cout << " * Assembling rhs: " << sensitivityMap_->cpuTimeSingleAssembling() << " s" << std::endl;
				std::cout << " * Solving:        " << sensitivityMap_->cpuTimeSingleSolution() << " s" << std::endl;

				std::cout << "-----------------------------------------------------------------------------" << std::endl;
				std::cout << std::endl;
			}
		}
	}

	void SurfacePerfectlyStirredReactor::FinalSummary(const boost::filesystem::path summary_file, const double tf)
	{
		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_);
		Equations(tf, yf_, dummy);

		if (psr_options_.verbose_output() == true)
		{
			std::ofstream fOut;
			fOut.open(summary_file.c_str(), std::ios::out);

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Status" << std::setw(20) << std::left << "Inlet" << "Outlet" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Time[s]" << std::setw(20) << std::left << 0. << tau_ << std::endl;
			fOut << std::setw(30) << std::left << "T[K]" << std::setw(20) << std::left << std::fixed << std::setprecision(3) << TInlet_ << T_ << std::endl;
			fOut << std::setw(30) << std::left << "P[atm]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << PInlet_ / 101325. << P_ / 101325. << std::endl;
			fOut << std::setw(30) << std::left << "rho[kg/m3]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << rhoInlet_ << rho_ << std::endl;
			fOut << std::setw(30) << std::left << "MW[kg/kmol]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << MWInlet_ << MW_ << std::endl;
			
			// Equivalence ratio
			{
				const double phiInlet = thermodynamicsMap_.GetLocalEquivalenceRatioFromMoleFractions(x0_.GetHandle());
				const double phi = thermodynamicsMap_.GetLocalEquivalenceRatioFromMoleFractions(x_.GetHandle());
				fOut << std::setw(30) << std::left << "Phi[-]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << phiInlet << phi << std::endl;
			}

			// Enthalpy analysis
			{
				const double H_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x_.GetHandle());
				const double U_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(x_.GetHandle());

				fOut << std::setw(30) << std::left << "H[J/kg]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << HInlet_ / MWInlet_ << H_ / MW_ << std::endl;
				fOut << std::setw(30) << std::left << "U[J/kg]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << UInlet_ / MWInlet_ << U_ / MW_ << std::endl;
			}

			for (unsigned int i = 1; i <= NC_; i++)
			if (omega0_[i] > 0.)
			{
				std::string label = "Conv.(%) " + thermodynamicsMap_.NamesOfSpecies()[i - 1];
				fOut << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(6) << 0. << (omegaInlet_[i] - omega_[i]) / omegaInlet_[i] * 100. << std::endl;
			}

			// Atomic analysis
			{
				std::vector<double> sum_inlet(thermodynamicsMap_.elements().size());
				std::vector<double> sum_final(thermodynamicsMap_.elements().size());
				for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
				{
					sum_inlet[j] = 0.;
					sum_final[j] = 0.;
					for (unsigned int i = 0; i < NC_; i++)
					{
						sum_inlet[j] += thermodynamicsMap_.atomic_composition()(i, j) * xInlet_[i + 1];
						sum_final[j] += thermodynamicsMap_.atomic_composition()(i, j) * x_[i + 1];
					}
				}

				const double moles_inlet = mass_flow_rate_in_ / MWInlet_;
				const double moles_final = mass_flow_rate_out_ / MW_;

				fOut << std::setw(30) << std::left << "Total species (kmol/s) " << std::setw(20) << std::left << std::scientific << std::setprecision(6) << moles_inlet << moles_final << std::endl;

				for (unsigned int j = 0; j<thermodynamicsMap_.elements().size(); j++)
				if (sum_inlet[j] > 0.)
				{
					std::string label = "Element(kmol/s) " + thermodynamicsMap_.elements()[j];
					fOut << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(6) << sum_inlet[j] * moles_inlet << sum_final[j] * moles_final << std::endl;
				}
			}

			fOut << "-----------------------------------------------------------------------------" << std::endl;

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Status" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "volume[m3]" << std::setw(20) << std::left << std::left << std::scientific << std::setprecision(6) << V0_ << V_ << std::endl;
			fOut << std::setw(30) << std::left << "mass[kg]" << std::setw(20) << std::left << std::left << std::scientific << std::setprecision(6) << mass0_ << mass_ << std::endl;
			fOut << std::setw(30) << std::left << "mass flow rate (in)[kg/s]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << mass_flow_rate_in_0_ << mass_flow_rate_in_ << std::endl;
			fOut << std::setw(30) << std::left << "mass flow rate (out)[kg/s]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << mass_flow_rate_in_0_ << mass_flow_rate_out_ << std::endl;
			fOut << std::setw(30) << std::left << "mass flow rate (loss)[kg/s]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << 0. << mass_flow_rate_loss_surface_ << std::endl;
			fOut << std::setw(30) << std::left << "Gamma[kmol/m2]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Gamma0_.SumElements() << Gamma0_.SumElements() << std::endl;
			fOut << std::setw(30) << std::left << "Bulk[kg]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << massBulk0_.SumElements() << massBulk_.SumElements() << std::endl;

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Bulk species [kg]" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for (unsigned int j = 0; j<BULK_NC_; j++)
				fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap_.vector_names_bulk_species()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << massBulk0_[j+1] << massBulk_[j + 1] << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Site densities [kmol/m2]" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for (unsigned int j = 0; j<thermodynamicsSurfaceMap_.number_of_site_phases(0); j++)
				fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap_.matrix_names_site_phases()[0][j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Gamma0_[j + 1] << Gamma_[j + 1] << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Surface fractions" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for (unsigned int j = 0; j<thermodynamicsSurfaceMap_.number_of_site_species(); j++)
				fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap_.NamesOfSpecies()[j + NC_] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Z0_[j + 1] << Z_[j + 1] << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;


			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Mass fractions" << std::setw(20) << std::left << "Inlet" << "Outlet" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << omegaInlet_[j + 1] << omega_[j + 1] << std::endl;
			fOut << "--------------------------------------------------------------" << std::endl;

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Mole fractions" << std::setw(20) << std::left << "Inlet" << "Outlet" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << xInlet_[j + 1] << x_[j + 1] << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
		}
	}

	void SurfacePerfectlyStirredReactor::GetFinalStatus(double& T, double& P, OpenSMOKE::OpenSMOKEVectorDouble& omega, OpenSMOKE::OpenSMOKEVectorDouble& Z, OpenSMOKE::OpenSMOKEVectorDouble& Gamma, double& mass, OpenSMOKE::OpenSMOKEVectorDouble& massBulk)
	{
		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_);
		Equations(0, yf_, dummy);

		T = T_;
		P = P_;
		omega = omega_;

		Z = Z_;
		Gamma = Gamma_;
		mass = mass_;
		massBulk = massBulk_;
	}

}

