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
#include "math/OpenSMOKEVector.h"
#include "SurfaceBatchReactor_OdeInterfaces.h"
#include "SurfaceBatchReactor_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"
#include "math/external-ode-solvers/ODE_Parameters.h"

namespace OpenSMOKE
{
	void SurfaceBatchReactor::EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap& sensitivityMap, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfaceBatchReactor_Options& batch_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options)
	{
		FatalErrorMessage("SurfaceBatchReactor::EnableSensitivityAnalysis not yet implemented");

		sensitivityMap_ = &sensitivityMap;

		PrepareSensitivityXMLFiles(thermodynamicsMap, batch_options, sensitivity_options);
		
		ChangeDimensions(NE_, &scaling_Jp, true);
		ChangeDimensions(NE_, NE_, &J, true);
	}

	void SurfaceBatchReactor::MemoryAllocation()
	{
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &RfromGas_, true);
		ChangeDimensions(NC_, &RfromSurface_, true);
		
		ChangeDimensions(NE_, &y0_, true);
		ChangeDimensions(NE_, &yf_, true);

		ChangeDimensions(NC_, &x0_, true);

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

	void SurfaceBatchReactor::NumericalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEMatrixDouble& J)
	{
		FatalErrorMessage("SurfaceBatchReactor::NumericalJacobian not yet implemented");

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

	void SurfaceBatchReactor::PrepareXMLFile(const boost::filesystem::path output_file_xml, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap)
	{
		fXML_.open(output_file_xml.c_str(), std::ios::out);
		OpenSMOKE::SetXMLFile(fXML_);
		fXML_ << "<Type> HomogeneousReactor </Type>" << std::endl;

		unsigned int counter = 2;
		fXML_ << "<additional>" << std::endl;
		fXML_ << 6 << std::endl;
		fXML_ << "time [s] " << counter++ << std::endl;
		fXML_ << "temperature [K] " << counter++ << std::endl;
		fXML_ << "pressure [Pa] " << counter++ << std::endl;
		fXML_ << "mol-weight [kg/kmol] " << counter++ << std::endl;
		fXML_ << "density [kg/m3] " << counter++ << std::endl;
		fXML_ << "heat-release [W/m3] " << counter++ << std::endl;
		fXML_ << "heat-release [W/m2] " << counter++ << std::endl;
		fXML_ << "</additional>" << std::endl;

		fXML_ << "<t-p-mw>" << std::endl;
		fXML_ << 1 << " " << 2 << " " << 3 << std::endl;
		fXML_ << "</t-p-mw>" << std::endl;

		fXML_ << "<mass-fractions>" << std::endl;
		fXML_ << thermodynamicsMap.NumberOfSpecies() << std::endl;
		for (unsigned int j=0;j<NC_;j++)
			fXML_ << thermodynamicsMap.NamesOfSpecies()[j] << " " << thermodynamicsMap.MW(j) << " " << counter++ << std::endl;
		fXML_ << "</mass-fractions>" << std::endl;

		fXML_ << "<surface-fractions>" << std::endl;
		fXML_ << thermodynamicsMap.NumberOfSpecies() << std::endl;
		for (unsigned int j=0;j<SURF_NC_;j++)
			fXML_ << thermodynamicsMap.vector_names_site_species()[j] << " " << thermodynamicsMap.MW(NC_+j) << " " << counter++ << std::endl;
		fXML_ << "</surface-fractions>" << std::endl;

		fXML_ << "<profiles>" << std::endl;
	}

	void SurfaceBatchReactor::CloseXMLFile()
	{
		fXML_ << "</profiles>" << std::endl;
		fXML_ << "<profiles-size> " << std::endl;
		fXML_ << counter_file_XML_ << " " << 1 + (NC_+1) << std::endl;
		fXML_ << "</profiles-size> " << std::endl;
		fXML_ << "</opensmoke>" << std::endl;
	}

	void SurfaceBatchReactor::PrepareSensitivityXMLFiles(OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfaceBatchReactor_Options& batch_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options)
	{
		FatalErrorMessage("SurfaceBatchReactor::PrepareSensitivityXMLFiles not yet implemented");

		indices_of_sensitivity_species_.resize(sensitivity_options.list_of_species().size());
		for(unsigned int i=0;i<indices_of_sensitivity_species_.size();i++)
			indices_of_sensitivity_species_[i] = thermodynamicsMap.IndexOfSpecies(sensitivity_options.list_of_species()[i]);

		const boost::filesystem::path parent_file = batch_options.output_path() / "Sensitivities.xml";
		fSensitivityParentXML_.open(parent_file.c_str(), std::ios::out);
		SetXMLFile(fSensitivityParentXML_);
		fSensitivityParentXML_ << "<variables>" << std::endl;
		
		if (type_ == SURFACEBATCH_REACTOR_NONISOTHERMAL_CONSTANTP || type_ == SURFACEBATCH_REACTOR_NONISOTHERMAL_CONSTANTV)
			fSensitivityParentXML_ << indices_of_sensitivity_species_.size()+1 << std::endl;
		else if (type_ == SURFACEBATCH_REACTOR_ISOTHERMAL_CONSTANTP || type_ == SURFACEBATCH_REACTOR_ISOTHERMAL_CONSTANTV)
			fSensitivityParentXML_ << indices_of_sensitivity_species_.size() << std::endl;

		for (unsigned int j=0;j<indices_of_sensitivity_species_.size();j++)
			fSensitivityParentXML_ << thermodynamicsMap.NamesOfSpecies()[indices_of_sensitivity_species_[j]-1] << " " << j << " " << indices_of_sensitivity_species_[j] << std::endl;
		
		if (type_ == SURFACEBATCH_REACTOR_NONISOTHERMAL_CONSTANTP || type_ == SURFACEBATCH_REACTOR_NONISOTHERMAL_CONSTANTV)
			fSensitivityParentXML_ << "temperature" << " " << indices_of_sensitivity_species_.size() << " " << NC_+1 << std::endl;
		
		fSensitivityParentXML_ << "</variables>" << std::endl;
		fSensitivityParentXML_ << "<n-parameters> " << std::endl;
		fSensitivityParentXML_ << sensitivityMap_->number_of_parameters() << std::endl;
		fSensitivityParentXML_ << "</n-parameters> " << std::endl;

		fSensitivityChildXML_ = new std::ofstream[indices_of_sensitivity_species_.size()+1];
		for (unsigned int j=0;j<indices_of_sensitivity_species_.size();j++)
		{
			const std::string name = "Sensitivities." + thermodynamicsMap.NamesOfSpecies()[indices_of_sensitivity_species_[j]-1] +".xml";
			const boost::filesystem::path child_file = batch_options.output_path() / name;
			fSensitivityChildXML_[j].open(child_file.c_str(), std::ios::out);
		}
		{
			const boost::filesystem::path child_file = batch_options.output_path() / "Sensitivities.temperature.xml";
			fSensitivityChildXML_[indices_of_sensitivity_species_.size()].open(child_file.c_str(), std::ios::out);
		}

		for (unsigned int j=0;j<indices_of_sensitivity_species_.size()+1;j++)
		{
			SetXMLFile(fSensitivityChildXML_[j]);
			fSensitivityChildXML_[j] << std::setprecision(5);
			fSensitivityChildXML_[j] << "<coefficients>" << std::endl;
		}		
	}

	void SurfaceBatchReactor::CloseSensitivityXMLFiles()
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

	void SurfaceBatchReactor::PrepareASCIIFile(const boost::filesystem::path output_file_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, OpenSMOKE::SurfaceBatchReactor_Options& batch_options)
	{
		indices_of_output_species_.resize(batch_options.output_species().size());
		for(unsigned int i=0;i<batch_options.output_species().size();i++)
			indices_of_output_species_[i] = thermodynamicsMap.IndexOfSpecies(batch_options.output_species()[i]);
		
		if (indices_of_output_species_.size() != 0)	
		{		
			widths_of_output_species_.resize(batch_options.output_species().size());
			for(unsigned int i=0;i<batch_options.output_species().size();i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(batch_options.output_species()[i], NC_);
		}
		else
		{
			widths_of_output_species_.resize(NC_);
			for(unsigned int i=0;i<NC_;i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap.NamesOfSpecies()[i], NC_);
		}

		{
			widths_of_output_surface_species_.resize(SURF_NC_);
			for(unsigned int i=0;i<SURF_NC_;i++)
				widths_of_output_surface_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap.NamesOfSpecies()[NC_+i], NC_);
		}
			
		fASCII_.open(output_file_ascii.c_str(), std::ios::out);

		unsigned int counter = 1;
		OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, "t[s]",		counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, "T[K]",		counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, "P[Pa]",		counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, "V[m3]",		counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, "rho[kg/m3]",	counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, "MW[kg/kmol]",	counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, "Qgas[W/m3]",	counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, "Qsurf[W/m2]",	counter);

		if (BULK_NC_ != 0)
		{
			OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, "bulk[kg]", counter);

			for (unsigned int i = 0; i < BULK_NC_; i++)
			{
				const std::string name = thermodynamicsMap.vector_names_bulk_species()[i].substr(0, 11) + "[kg]";
				OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, name, counter);
			}
		}

		for(unsigned int i=0;i<thermodynamicsMap.number_of_site_phases(0);i++)
		{
			const std::string name = thermodynamicsMap.matrix_names_site_phases()[0][i].substr(0,13);
			OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, name,	counter);
		}

		for(unsigned int i=0;i<SURF_NC_;i++)
			OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_surface_species_[i], fASCII_,  thermodynamicsMap.NamesOfSpecies()[NC_+i] + "_z", counter);

		if (indices_of_output_species_.size() != 0)
		{
			for(unsigned int i=0;i<indices_of_output_species_.size();i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fASCII_,  thermodynamicsMap.NamesOfSpecies()[indices_of_output_species_[i]-1] + "_x", counter);
			for(unsigned int i=0;i<indices_of_output_species_.size();i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fASCII_,  thermodynamicsMap.NamesOfSpecies()[indices_of_output_species_[i]-1] + "_w", counter);
		}
		else
		{
			for(unsigned int i=0;i<NC_;i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fASCII_,  thermodynamicsMap.NamesOfSpecies()[i] + "_x", counter);
			for(unsigned int i=0;i<NC_;i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fASCII_,  thermodynamicsMap.NamesOfSpecies()[i] + "_w", counter);
		}

		fASCII_ << std::endl;
	}

	void SurfaceBatchReactor::PrepareBulkASCIIFile(const boost::filesystem::path output_file_bulk_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, OpenSMOKE::SurfaceBatchReactor_Options& batch_options)
	{
		fBulkASCII_.open(output_file_bulk_ascii.c_str(), std::ios::out);

		unsigned int counter = 1;
		OpenSMOKE::PrintTagOnASCIILabel(20, fBulkASCII_, "t[s]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fBulkASCII_, "T[K]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fBulkASCII_, "P[Pa]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fBulkASCII_, "V[m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fBulkASCII_, "rho[kg/m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fBulkASCII_, "MW[kg/kmol]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fBulkASCII_, "Qgas[W/m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fBulkASCII_, "Qsurf[W/m2]", counter);

		OpenSMOKE::PrintTagOnASCIILabel(24, fBulkASCII_, "total[kg]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(24, fBulkASCII_, "total[kg/m2/s]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(24, fBulkASCII_, "total[kg/s]", counter);

		for (unsigned int i = 0; i < BULK_NC_; i++)
		{
			const std::string name_species = thermodynamicsMap.vector_names_bulk_species()[i].substr(0, 9);
			
			OpenSMOKE::PrintTagOnASCIILabel(24, fBulkASCII_, name_species + "[kg]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fBulkASCII_, name_species + "[kg/m2/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fBulkASCII_, name_species + "[kg/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fBulkASCII_, name_species + "[m]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fBulkASCII_, name_species + "[micron]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fBulkASCII_, name_species + "[m/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fBulkASCII_, name_species + "[mm/hr]", counter);
		}

		fBulkASCII_ << std::endl;
	}

	void SurfaceBatchReactor::SolveOpenSourceSolvers(const double tf, OpenSMOKE::ODE_Parameters& ode_parameters)
	{
		#if OPENSMOKE_USE_DVODE == 1
		if (ode_parameters.type() == ODE_Parameters::ODE_INTEGRATOR_DVODE)
		{
			// ODE system
			ODESystem_DVODE_SurfaceBatchReactor *odeSystemObject;
			odeSystemObject = ODESystem_DVODE_SurfaceBatchReactor::GetInstance();
			odeSystemObject->SetBatchReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DVODE<ODESystem_DVODE_SurfaceBatchReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters.SetCPUTime(tEnd-tStart);
			ode_parameters.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif

		#if OPENSMOKE_USE_ODEPACK == 1
		if (ode_parameters.type() == ODE_Parameters::ODE_INTEGRATOR_DLSODE)
		{
			// ODE system
			ODESystem_DLSODE_SurfaceBatchReactor *odeSystemObject;
			odeSystemObject = ODESystem_DLSODE_SurfaceBatchReactor::GetInstance();
			odeSystemObject->SetBatchReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DLSODE<ODESystem_DLSODE_SurfaceBatchReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters.SetCPUTime(tEnd-tStart);
			ode_parameters.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters.SetLastStepUsed(o.GetLastStepUsed());
		}

		if (ode_parameters.type() == ODE_Parameters::ODE_INTEGRATOR_DLSODA)
		{
			// ODE system
			ODESystem_DLSODA_SurfaceBatchReactor *odeSystemObject;
			odeSystemObject = ODESystem_DLSODA_SurfaceBatchReactor::GetInstance();
			odeSystemObject->SetBatchReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DLSODA<ODESystem_DLSODA_SurfaceBatchReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters.SetCPUTime(tEnd-tStart);
			ode_parameters.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif

		#if OPENSMOKE_USE_DASPK == 1
		if (ode_parameters.type() == ODE_Parameters::ODE_INTEGRATOR_DASPK)
		{
			// ODE system
			ODESystem_DASPK_SurfaceBatchReactor *odeSystemObject;
			odeSystemObject = ODESystem_DASPK_SurfaceBatchReactor::GetInstance(NE_);
			odeSystemObject->SetBatchReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DASPK<ODESystem_DASPK_SurfaceBatchReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters.SetCPUTime(tEnd-tStart);
			ode_parameters.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif

		#if OPENSMOKE_USE_RADAU == 1
		if (ode_parameters.type() == ODE_Parameters::ODE_INTEGRATOR_RADAU5)
		{
			// ODE system
			ODESystem_RADAU5_SurfaceBatchReactor *odeSystemObject;
			odeSystemObject = ODESystem_RADAU5_SurfaceBatchReactor::GetInstance();
			odeSystemObject->SetBatchReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_RADAU<ODESystem_RADAU5_SurfaceBatchReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters.SetCPUTime(tEnd-tStart);
			ode_parameters.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif

		#if OPENSMOKE_USE_SUNDIALS == 1
		if (ode_parameters.type() == ODE_Parameters::ODE_INTEGRATOR_CVODE)
		{
			// ODE system
			ODESystem_CVODE_SurfaceBatchReactor *odeSystemObject;
			odeSystemObject = ODESystem_CVODE_SurfaceBatchReactor::GetInstance();
			odeSystemObject->SetBatchReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_CVODE_Sundials<ODESystem_CVODE_SurfaceBatchReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters.SetCPUTime(tEnd-tStart);
			ode_parameters.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif

		#if OPENSMOKE_USE_MEBDF == 1
		if (ode_parameters.type() == ODE_Parameters::ODE_INTEGRATOR_MEBDF)
		{
			// ODE system
			ODESystem_MEBDF_SurfaceBatchReactor *odeSystemObject;
			odeSystemObject = ODESystem_MEBDF_SurfaceBatchReactor::GetInstance();
			odeSystemObject->SetBatchReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_MEBDF<ODESystem_MEBDF_SurfaceBatchReactor>   o(odeSystemObject);
			o.SetDimensions(NE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., y0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yf_.GetHandle());

			ode_parameters.SetCPUTime(tEnd-tStart);
			ode_parameters.SetNumberOfSteps(o.GetNumberOfSteps());
			ode_parameters.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			ode_parameters.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			ode_parameters.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			ode_parameters.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			ode_parameters.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			ode_parameters.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			ode_parameters.SetLastOrderUsed(o.GetLastOrderUsed());
			ode_parameters.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif
	}

	void SurfaceBatchReactor::OpenAllFiles(OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, OpenSMOKE::SurfaceBatchReactor_Options& batch_options)
	{
		if (batch_options.verbose_output() == true)
		{
			if ( !boost::filesystem::exists(batch_options.output_path()) )
				OpenSMOKE::CreateDirectory(batch_options.output_path());

			if (batch_options.verbose_ascii_file() == true)
				PrepareASCIIFile(batch_options.output_path() / "Output.out", thermodynamicsMap, batch_options);
			
			if (batch_options.verbose_xml_file() == true)
				PrepareXMLFile(batch_options.output_path() / "Output.xml", thermodynamicsMap);

			if (batch_options.verbose_ascii_file() == true && BULK_NC_ != 0)
				PrepareBulkASCIIFile(batch_options.output_path() / "Output.bulk.out", thermodynamicsMap, batch_options);
		}
	}

	void SurfaceBatchReactor::CloseAllFiles(OpenSMOKE::SurfaceBatchReactor_Options& batch_options)
	{
		if (batch_options.verbose_output() == true)
			{
				if (batch_options.verbose_ascii_file() == true)
					fASCII_.close();

				if (batch_options.verbose_xml_file() == true)
					CloseXMLFile();
			}

		if (batch_options.sensitivity_analysis() == true)
			CloseSensitivityXMLFiles();
	}

	void SurfaceBatchReactor::FinalStatus(const double tf,	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
															OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap)
	{
		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_);
		Equations(tf, yf_, dummy);

		std::cout << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << std::setw(30) << std::left << "Status" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << std::setw(30) << std::left << "Time[s]"		<< std::setw(20) << std::left << 0.  << tf << std::endl;
		std::cout << std::setw(30) << std::left << "T[K]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(3) << T0_ << T_ << std::endl;
		std::cout << std::setw(30) << std::left << "P[atm]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << P0_/101325. << P_/101325. << std::endl;
		std::cout << std::setw(30) << std::left << "V[m3]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << V0_ << V_ << std::endl;
		std::cout << std::setw(30) << std::left << "A[m2]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << A0_ << A_ << std::endl;
		std::cout << std::setw(30) << std::left << "mass[kg]"		<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << mass0_ << mass_ << std::endl;
		std::cout << std::setw(30) << std::left << "rho[kg/m3]"		<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << mass_/V0_ << mass_/V_ << std::endl;
		std::cout << std::setw(30) << std::left << "MW[kg/kmol]"	<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << MW0_ << MW_ << std::endl;
		std::cout << std::setw(30) << std::left << "Gamma[kmol/m2]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Gamma0_.SumElements() << Gamma0_.SumElements() << std::endl;
		std::cout << std::setw(30) << std::left << "Bulk[kg]"       << std::setw(20) << std::left << std::scientific << std::setprecision(6) << massBulk0_ << massBulk_.SumElements() << std::endl;

		// Equivalence ratio
		{
			const double phi0 = thermodynamicsMap.GetLocalEquivalenceRatioFromMoleFractions(x0_.GetHandle());
			const double phi = thermodynamicsMap.GetLocalEquivalenceRatioFromMoleFractions(x_.GetHandle());
			std::cout << std::setw(30) << std::left << "Phi[-]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << phi0 << phi << std::endl;
		}

		// Enthalpy analysis
		{
			const double H_ = thermodynamicsMap.hMolar_Mixture_From_MoleFractions(x_.GetHandle());
			const double U_ = thermodynamicsMap.uMolar_Mixture_From_MoleFractions(x_.GetHandle());
			
			std::cout << std::setw(30) << std::left << "Gas H [J]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_/MW0_*mass0_ << H_/MW_*mass_ << std::endl;
			std::cout << std::setw(30) << std::left << "Gas U [J]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << U0_/MW0_*mass0_ << U_/MW_*mass_ << std::endl;

			double H_Surface_ = 0;
			double U_Surface_ = 0;
			for(unsigned int i=1;i<=SURF_NC_;i++)
			{
				const unsigned int index_phase = thermodynamicsSurfaceMap.vector_site_phases_belonging()[i-1]+1;
				H_Surface_ += Z_[i] * thermodynamicsSurfaceMap.Species_H_over_RT()[i+NC_-1]*PhysicalConstants::R_J_kmol*T_ *
							  Gamma_[index_phase] * A_;
				U_Surface_ += Z_[i] * (thermodynamicsSurfaceMap.Species_H_over_RT()[i+NC_-1]-1.)*PhysicalConstants::R_J_kmol*T_ *
							  Gamma_[index_phase] * A_;
			}

			std::cout << std::setw(30) << std::left << "Surface H [J]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_Surface_ << H_Surface_ << std::endl;
			std::cout << std::setw(30) << std::left << "Surface U [J]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << U0_Surface_ << U_Surface_ << std::endl;
			
			std::cout << std::setw(30) << std::left << "Total H [J]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_/MW0_*mass0_ + H0_Surface_ << H_/MW_*mass_ + H_Surface_ << std::endl;
			std::cout << std::setw(30) << std::left << "Total U [J]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << U0_/MW0_*mass0_ + U0_Surface_ << U_/MW_*mass_ + U_Surface_ << std::endl;
		}

		// Atomic analysis
		{
			const double gas_moles_initial = P0_*V0_/(PhysicalConstants::R_J_kmol*T0_);
			const double gas_moles_final = P_*V_/(PhysicalConstants::R_J_kmol*T_);

			double surface_moles_initial = 0.;
			double surface_moles_final = 0.;

			std::vector<double> sum_initial(thermodynamicsSurfaceMap.elements().size());
			std::vector<double> sum_final(thermodynamicsSurfaceMap.elements().size());
			for(unsigned int j=0;j<thermodynamicsSurfaceMap.elements().size();j++)
			{
				sum_initial[j]	= 0.;
				sum_final[j]	= 0.;
				for(unsigned int i=0;i<NC_;i++)
				{
					sum_initial[j] += thermodynamicsSurfaceMap.atomic_composition()(i,j) * x0_[i+1] * gas_moles_initial;
					sum_final[j] += thermodynamicsSurfaceMap.atomic_composition()(i,j) * x_[i+1] * gas_moles_final;
				}
				for(unsigned int i=0;i<SURF_NC_;i++)
				{
					const unsigned int index_phase = thermodynamicsSurfaceMap.vector_site_phases_belonging()[i]+1;
					sum_initial[j] += thermodynamicsSurfaceMap.atomic_composition()(i+NC_,j) * Z0_[i+1] * Gamma0_[index_phase] * A0_;
					sum_final[j] += thermodynamicsSurfaceMap.atomic_composition()(i+NC_,j) * Z_[i+1] * Gamma_[index_phase] * A_;
					surface_moles_initial += Gamma0_[index_phase] * A0_;
					surface_moles_final += Gamma_[index_phase] * A_;
				}
			}

			for(unsigned int i=1;i<=NC_;i++)
			if (omega0_[i] > 0.)
			{
				const double m0 = omega0_[i]*mass0_;
				const double mf = omega_[i]*mass_;
				std::string label = "Gas Conversion (%) " + thermodynamicsMap.NamesOfSpecies()[i-1];
				std::cout << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(3) << 0. << (m0-mf)/m0*100. << std::endl;
			}
			

			std::cout << std::setw(30) << std::left << "Total gas species (kmol) "		<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << gas_moles_initial << gas_moles_final << std::endl;
			std::cout << std::setw(30) << std::left << "Total surface species (kmol) "	<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << surface_moles_initial << surface_moles_final << std::endl;

			for(unsigned int j=0;j<thermodynamicsMap.elements().size();j++)
				if (sum_initial[j] > 0.)
				{
					std::string label = "Gas element (kmol) " + thermodynamicsMap.elements()[j];
					std::cout << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(6) << sum_initial[j] << sum_final[j] << std::endl;
				}
		}

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
	}

	void SurfaceBatchReactor::FinalSummary(const boost::filesystem::path summary_file, const double tf, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
																										OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap)
	{
		std::ofstream fOut;
		fOut.open(summary_file.c_str(), std::ios::out);

		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_);
		Equations(tf, yf_, dummy);

		fOut << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		fOut << std::setw(30) << std::left << "Status" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		fOut << std::setw(30) << std::left << "Time[s]"			<< std::setw(20) << std::left << 0.  << tf << std::endl;
		fOut << std::setw(30) << std::left << "T[K]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(3) << T0_ << T_ << std::endl;
		fOut << std::setw(30) << std::left << "P[atm]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << P0_/101325. << P_/101325. << std::endl;
		fOut << std::setw(30) << std::left << "V[m3]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << V0_ << V_ << std::endl;
		fOut << std::setw(30) << std::left << "A[m2]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << A0_ << A_ << std::endl;
		fOut << std::setw(30) << std::left << "mass[kg]"		<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << mass0_ << mass_ << std::endl;
		fOut << std::setw(30) << std::left << "rho[kg/m3]"		<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << rho_ << rho_ << std::endl;
		fOut << std::setw(30) << std::left << "MW[kg/kmol]"		<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << MW0_ << MW_ << std::endl;
		fOut << std::setw(30) << std::left << "Gamma[kmol/m2]"	<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << Gamma0_.SumElements() << Gamma0_.SumElements() << std::endl;
		fOut << std::setw(30) << std::left << "Bulk[kg]"		<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << massBulk0_ << massBulk_.SumElements() << std::endl;

		// Equivalence ratio
		{
			const double phi0 = thermodynamicsMap.GetLocalEquivalenceRatioFromMoleFractions(x0_.GetHandle());
			const double phi = thermodynamicsMap.GetLocalEquivalenceRatioFromMoleFractions(x_.GetHandle());
			fOut << std::setw(30) << std::left << "Phi[-]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << phi0 << phi << std::endl;
		}

		// Enthalpy analysis
		{
			const double H_ = thermodynamicsMap.hMolar_Mixture_From_MoleFractions(x_.GetHandle());
			const double U_ = thermodynamicsMap.uMolar_Mixture_From_MoleFractions(x_.GetHandle());
			
			fOut << std::setw(30) << std::left << "Gas H [J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_/MW0_ << H_/MW_ << std::endl;
			fOut << std::setw(30) << std::left << "Gas U [J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << U0_/MW0_ << U_/MW_ << std::endl;
		}

		// Atomic analysis
		{
			const double gas_moles_initial = P0_*V0_/(PhysicalConstants::R_J_kmol*T0_);
			const double gas_moles_final = P_*V_/(PhysicalConstants::R_J_kmol*T_);

			double surface_moles_initial = 0.;
			double surface_moles_final = 0.;

			std::vector<double> sum_initial(thermodynamicsSurfaceMap.elements().size());
			std::vector<double> sum_final(thermodynamicsSurfaceMap.elements().size());
			for(unsigned int j=0;j<thermodynamicsSurfaceMap.elements().size();j++)
			{
				sum_initial[j]	= 0.;
				sum_final[j]	= 0.;
				for(unsigned int i=0;i<NC_;i++)
				{
					sum_initial[j] += thermodynamicsSurfaceMap.atomic_composition()(i,j) * x0_[i+1] * gas_moles_initial;
					sum_final[j] += thermodynamicsSurfaceMap.atomic_composition()(i,j) * x_[i+1] * gas_moles_final;
				}
				for(unsigned int i=0;i<SURF_NC_;i++)
				{
					const unsigned int index_phase = thermodynamicsSurfaceMap.vector_site_phases_belonging()[i]+1;
					sum_initial[j] += thermodynamicsSurfaceMap.atomic_composition()(i+NC_,j) * Z0_[i+1] * Gamma0_[index_phase] * A0_;
					sum_final[j] += thermodynamicsSurfaceMap.atomic_composition()(i+NC_,j) * Z_[i+1] * Gamma_[index_phase] * A_;
					surface_moles_initial += Gamma0_[index_phase] * A0_;
					surface_moles_final += Gamma_[index_phase] * A_;
				}
			}

			for(unsigned int i=1;i<=NC_;i++)
			if (omega0_[i] > 0.)
			{
				const double m0 = omega0_[i]*mass0_;
				const double mf = omega_[i]*mass_;
				std::string label = "Gas Conversion (%) " + thermodynamicsMap.NamesOfSpecies()[i-1];
				fOut << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(3) << 0. << (m0-mf)/m0*100. << std::endl;
			}
			

			fOut << std::setw(30) << std::left << "Total gas species (kmol) "		<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << gas_moles_initial << gas_moles_final << std::endl;
			fOut << std::setw(30) << std::left << "Total surface species (kmol) "	<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << surface_moles_initial << surface_moles_final << std::endl;

			for(unsigned int j=0;j<thermodynamicsMap.elements().size();j++)
				if (sum_initial[j] > 0.)
				{
					std::string label = "Gas element (kmol) " + thermodynamicsMap.elements()[j];
					fOut << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(6) << sum_initial[j] << sum_final[j] << std::endl;
				}
		}

		fOut << "-----------------------------------------------------------------------------" << std::endl;

		fOut << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		fOut << std::setw(30) << std::left << "Bulk species [kg]" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		for (unsigned int j = 0; j<BULK_NC_; j++)
			fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap.vector_names_bulk_species()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << massBulk0_ << massBulk_[j + 1] << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;


		fOut << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		fOut << std::setw(30) << std::left << "Site densities [kmol/m2]" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		for(unsigned int j=0;j<thermodynamicsSurfaceMap.number_of_site_phases(0);j++)
			fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap.matrix_names_site_phases()[0][j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Gamma0_[j+1] << Gamma_[j+1] << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;

		fOut << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		fOut << std::setw(30) << std::left << "Surface fractions" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		for(unsigned int j=0;j<thermodynamicsSurfaceMap.number_of_site_species();j++)
			fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap.NamesOfSpecies()[j+NC_] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Z0_[j+1] << Z_[j+1] << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;

		fOut << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		fOut << std::setw(30) << std::left << "Mass fractions" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		for(unsigned int j=0;j<thermodynamicsMap.NumberOfSpecies();j++)
			fOut << std::setw(30) << std::left << thermodynamicsMap.NamesOfSpecies()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << omega0_[j+1] << omega_[j+1] << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		
		fOut << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		fOut << std::setw(30) << std::left << "Mole fractions" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		for(unsigned int j=0;j<thermodynamicsMap.NumberOfSpecies();j++)
			fOut << std::setw(30) << std::left << thermodynamicsMap.NamesOfSpecies()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << x0_[j+1] << x_[j+1] << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
	}

	void SurfaceBatchReactor::GetFinalStatus(double& T, double& P, OpenSMOKE::OpenSMOKEVectorDouble& omega, OpenSMOKE::OpenSMOKEVectorDouble& Z, OpenSMOKE::OpenSMOKEVectorDouble& Gamma, double& mass, OpenSMOKE::OpenSMOKEVectorDouble& massBulk)
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

