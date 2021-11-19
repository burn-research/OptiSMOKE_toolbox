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
#include "PlugFlowReactor_OdeInterfaces.h"
#include "PlugFlowReactor_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"
#include "math/external-ode-solvers/ODE_Parameters.h"

namespace OpenSMOKE
{
	PlugFlowReactor::PlugFlowReactor(
		OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
		OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
		OpenSMOKE::ODE_Parameters& ode_parameters,
		OpenSMOKE::PlugFlowReactor_Options& plugflow_options,
		OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa,
		OpenSMOKE::OnTheFlyPostProcessing& on_the_fly_post_processing,
		OpenSMOKE::IgnitionDelayTimes_Analyzer& idts_analyzer,
		OpenSMOKE::PolimiSoot_Analyzer& polimi_soot_analyzer) :
		thermodynamicsMap_(thermodynamicsMap),
		kineticsMap_(kineticsMap),
		ode_parameters_(ode_parameters),
		plugflow_options_(plugflow_options),
		on_the_fly_ropa_(on_the_fly_ropa),
		on_the_fly_post_processing_(on_the_fly_post_processing),
		idts_analyzer_(idts_analyzer),
		polimi_soot_analyzer_(polimi_soot_analyzer)
	{
	}

	void PlugFlowReactor::EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap& sensitivityMap, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options)
	{
		iXmlMaps_ = true;
		
		sensitivityMap_ = &sensitivityMap;

		PrepareSensitivityXMLFiles(sensitivity_options);
		
		ChangeDimensions(NE_, &scaling_Jp, true);
		ChangeDimensions(NE_, NE_, &J, true);
	}

	void PlugFlowReactor::EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap& sensitivityMap, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options, bool iXmlMaps)
	{
		iXmlMaps_ = iXmlMaps;

		sensitivityMap_ = &sensitivityMap;

		if(iXmlMaps == true)
		    PrepareSensitivityXMLFiles(sensitivity_options);
		
		ChangeDimensions(NE_, &scaling_Jp, true);
		ChangeDimensions(NE_, NE_, &J, true);
	}

	void PlugFlowReactor::MemoryAllocation()
	{
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &R_, true);
		
		ChangeDimensions(NE_, &y0_, true);
		ChangeDimensions(NE_, &yf_, true);
		ChangeDimensions(NC_, &x0_, true);
	}

	void PlugFlowReactor::NumericalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEMatrixDouble& J)
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

	void PlugFlowReactor::PrepareROPAFile(const boost::filesystem::path output_file_ropa)
	{
		// Open ROPA file
		fROPA_.open(output_file_ropa.c_str(), std::ios::out);

		// Write head in file
		on_the_fly_ropa_.WriteHead(fROPA_, "Plug Flow Reactor");

		// ROPA at initial conditions
		const double cTot0_ = rho0_ / MW0_;
		OpenSMOKE::OpenSMOKEVectorDouble c0(thermodynamicsMap_.NumberOfSpecies());
		OpenSMOKE::Product(cTot0_, x0_, &c0);
		on_the_fly_ropa_.Analyze(fROPA_, 0, 0., T0_, P0_, c0, omega0_, omega0_);
	}

	void PlugFlowReactor::PreparePolimiSootFiles(const boost::filesystem::path output_file_polimi_soot, const boost::filesystem::path output_file_polimi_soot_distribution)
	{
		// Polimi Soot Analysis
		{
			fPolimiSoot_.open(output_file_polimi_soot.c_str(), std::ios::out);
			fPolimiSoot_.setf(std::ios::scientific);
			fPolimiSoot_.setf(std::ios::left);

			const unsigned int width = 20;
			unsigned int counter = 1;
			OpenSMOKE::PrintTagOnASCIILabel(width, fPolimiSoot_, "t[s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(width, fPolimiSoot_, "csi[m]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(width, fPolimiSoot_, "T0[K]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(width, fPolimiSoot_, "P0[Pa]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(width, fPolimiSoot_, "T[K]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(width, fPolimiSoot_, "P[Pa]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(width, fPolimiSoot_, "v[m/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(width, fPolimiSoot_, "rho[kg/m3]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(width, fPolimiSoot_, "MW[kg/kmol]", counter);

			polimi_soot_analyzer_.WriteLabelIntegralDataFile(fPolimiSoot_, counter, width);
			fPolimiSoot_ << std::endl;
		}

		// Polimi Soot Particle Size Distribution
		if (polimi_soot_analyzer_.write_psdf() == true)
		{
			fPolimiSootDistribution_.open(output_file_polimi_soot_distribution.c_str(), std::ios::out);
			fPolimiSootDistribution_.setf(std::ios::scientific);
			fPolimiSootDistribution_.setf(std::ios::left);

			polimi_soot_analyzer_.WriteDistributionLabel(fPolimiSootDistribution_);
		}
	}

	void PlugFlowReactor::PrepareXMLFile(const boost::filesystem::path output_file_xml)
	{
		fXML_.open(output_file_xml.c_str(), std::ios::out);
		OpenSMOKE::SetXMLFile(fXML_);
		fXML_ << "<Type> PlugFlowReactor </Type>" << std::endl;

		unsigned int counter = 2;
		fXML_ << "<additional>" << std::endl;
		fXML_ << 8 << std::endl;
		if (time_independent_ == true)
		{
			fXML_ << "time [s] " << counter++ << std::endl;
			fXML_ << "axial-coordinate [m] " << counter++ << std::endl;
		}
		else
		{
			fXML_ << "axial-coordinate [m] " << counter++ << std::endl;
			fXML_ << "time [s] " << counter++ << std::endl;
		}
		fXML_ << "temperature [K] " << counter++ << std::endl;
		fXML_ << "pressure [Pa] " << counter++ << std::endl;
		fXML_ << "mol-weight [kg/kmol] " << counter++ << std::endl;
		fXML_ << "density [kg/m3] " << counter++ << std::endl;
		fXML_ << "heat-release [W/m3] " << counter++ << std::endl;
		fXML_ << "velocity [m/s] " << counter++ << std::endl;
		fXML_ << "</additional>" << std::endl;

		fXML_ << "<t-p-mw>" << std::endl;
		fXML_ << 2 << " " << 3 << " " << 4 << std::endl;
		fXML_ << "</t-p-mw>" << std::endl;

		fXML_ << "<mass-fractions>" << std::endl;
		fXML_ << thermodynamicsMap_.NumberOfSpecies() << std::endl;
		for (unsigned int j=0;j<NC_;j++)
			fXML_ << thermodynamicsMap_.NamesOfSpecies()[j] << " " << thermodynamicsMap_.MW(j) << " " << counter++ << std::endl;
		fXML_ << "</mass-fractions>" << std::endl;
		fXML_ << "<profiles>" << std::endl;
	}

	void PlugFlowReactor::CloseXMLFile()
	{
		fXML_ << "</profiles>" << std::endl;
		fXML_ << "<profiles-size> " << std::endl;
		fXML_ << counter_file_XML_ << " " << 1 + (NC_+1) << std::endl;
		fXML_ << "</profiles-size> " << std::endl;
		fXML_ << "</opensmoke>" << std::endl;
	}

	void PlugFlowReactor::PrepareSensitivityXMLFiles(OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options)
	{
		indices_of_sensitivity_species_.resize(sensitivity_options.list_of_species().size());
		for(unsigned int i=0;i<indices_of_sensitivity_species_.size();i++)
			indices_of_sensitivity_species_[i] = thermodynamicsMap_.IndexOfSpecies(sensitivity_options.list_of_species()[i]);

		const boost::filesystem::path parent_file = plugflow_options_.output_path() / "Sensitivities.xml";
		fSensitivityParentXML_.open(parent_file.c_str(), std::ios::out);
		SetXMLFile(fSensitivityParentXML_);
		fSensitivityParentXML_ << "<variables>" << std::endl;
		
		if (type_ == PLUGFLOW_REACTOR_NONISOTHERMAL)
			fSensitivityParentXML_ << indices_of_sensitivity_species_.size()+1 << std::endl;
		else if (type_ == PLUGFLOW_REACTOR_ISOTHERMAL)
			fSensitivityParentXML_ << indices_of_sensitivity_species_.size() << std::endl;

		for (unsigned int j=0;j<indices_of_sensitivity_species_.size();j++)
			fSensitivityParentXML_ << thermodynamicsMap_.NamesOfSpecies()[indices_of_sensitivity_species_[j]-1] << " " << j << " " << indices_of_sensitivity_species_[j] << std::endl;
		
		if (type_ == PLUGFLOW_REACTOR_NONISOTHERMAL)
			fSensitivityParentXML_ << "temperature" << " " << indices_of_sensitivity_species_.size() << " " << NC_+1 << std::endl;
		
		fSensitivityParentXML_ << "</variables>" << std::endl;
		fSensitivityParentXML_ << "<n-parameters> " << std::endl;
		fSensitivityParentXML_ << sensitivityMap_->number_of_parameters() << std::endl;
		fSensitivityParentXML_ << "</n-parameters> " << std::endl;

		fSensitivityChildXML_ = new std::ofstream[indices_of_sensitivity_species_.size()+1];
		for (unsigned int j=0;j<indices_of_sensitivity_species_.size();j++)
		{
			const std::string name = "Sensitivities." + thermodynamicsMap_.NamesOfSpecies()[indices_of_sensitivity_species_[j]-1] +".xml";
			const boost::filesystem::path child_file = plugflow_options_.output_path() / name;
			fSensitivityChildXML_[j].open(child_file.c_str(), std::ios::out);
		}
		{
			const boost::filesystem::path child_file = plugflow_options_.output_path() / "Sensitivities.temperature.xml";
			fSensitivityChildXML_[indices_of_sensitivity_species_.size()].open(child_file.c_str(), std::ios::out);
		}

		for (unsigned int j=0;j<indices_of_sensitivity_species_.size()+1;j++)
		{
			SetXMLFile(fSensitivityChildXML_[j]);
			fSensitivityChildXML_[j] << std::setprecision(5);
			fSensitivityChildXML_[j] << "<coefficients>" << std::endl;
		}		
	}

	void PlugFlowReactor::CloseSensitivityXMLFiles()
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

	void PlugFlowReactor::PrepareASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii)
	{
		indices_of_output_species_.resize(plugflow_options_.output_species().size());
		for (unsigned int i = 0; i<plugflow_options_.output_species().size(); i++)
			indices_of_output_species_[i] = thermodynamicsMap_.IndexOfSpecies(plugflow_options_.output_species()[i]);

		if (indices_of_output_species_.size() != 0)
		{
			widths_of_output_species_.resize(plugflow_options_.output_species().size());
			for (unsigned int i = 0; i<plugflow_options_.output_species().size(); i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(plugflow_options_.output_species()[i], NC_);
		}
		else
		{
			widths_of_output_species_.resize(NC_);
			for (unsigned int i = 0; i<NC_; i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap_.NamesOfSpecies()[i], NC_);
		}

		fOutput.open(output_file_ascii.c_str(), std::ios::out);

		unsigned int counter = 1;
		fOutput.setf(std::ios::scientific);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "t[s]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "csi[m]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T0[K]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P0[Pa]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P[Pa]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "v[m/s]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "rho[kg/m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "MW[kg/kmol]", counter);

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

	void PlugFlowReactor::PrepareParametricASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii)
	{
		indices_of_output_species_.resize(plugflow_options_.output_species().size());
		for (unsigned int i = 0; i<plugflow_options_.output_species().size(); i++)
			indices_of_output_species_[i] = thermodynamicsMap_.IndexOfSpecies(plugflow_options_.output_species()[i]);

		if (indices_of_output_species_.size() != 0)
		{
			widths_of_output_species_.resize(plugflow_options_.output_species().size());
			for (unsigned int i = 0; i<plugflow_options_.output_species().size(); i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(plugflow_options_.output_species()[i], NC_);
		}
		else
		{
			widths_of_output_species_.resize(NC_);
			for (unsigned int i = 0; i<NC_; i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap_.NamesOfSpecies()[i], NC_);
		}

		fOutput.open(output_file_ascii.c_str(), std::ios::out);

		unsigned int counter = 1;
		fOutput.setf(std::ios::scientific);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "t[s]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "csi[m]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T0[K]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P0[Pa]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P[Pa]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "v[m/s]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "rho[kg/m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "MW[kg/kmol]", counter);

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

		if (indices_of_output_species_.size() != 0)
		{
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap_.NamesOfSpecies()[indices_of_output_species_[i] - 1] + "_x0", counter);
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap_.NamesOfSpecies()[indices_of_output_species_[i] - 1] + "_w0", counter);
		}
		else
		{
			for (unsigned int i = 0; i<NC_; i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap_.NamesOfSpecies()[i] + "_x0", counter);
			for (unsigned int i = 0; i<NC_; i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap_.NamesOfSpecies()[i] + "_w0", counter);
		}

		fOutput << std::endl;
	}

	void PlugFlowReactor::PrepareParametricIDTASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii)
	{
		if (idts_analyzer_.is_active() == true)
		{
			fOutput.open(output_file_ascii.c_str(), std::ios::out);

			unsigned int counter = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "t[s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "csi[m]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T0[K]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P0[Pa]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "v0[m/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "P[Pa]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "v[m/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "rho[kg/m3]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "MW[kg/kmol]", counter);

			if (idts_analyzer_.is_active() == true)
				idts_analyzer_.PrintHeaderLine(counter, fOutput);

			fOutput << std::endl;
		}
	}

	void PlugFlowReactor::PrepareASCIIFile(const boost::filesystem::path output_file_ascii)
	{
		PrepareASCIIFile(fASCII_, output_file_ascii);
	}

	void PlugFlowReactor::PrintParametricIDT(std::ostream& fOutput, const double t)
	{
		if (idts_analyzer_.is_active())
		{
			fOutput.setf(std::ios::scientific);
			fOutput << std::setw(20) << std::left << tau_;
			fOutput << std::setw(20) << std::left << csi_;
			fOutput << std::setw(20) << std::left << T0_;
			fOutput << std::setw(20) << std::left << P0_;
			fOutput << std::setw(20) << std::left << v0_;
			fOutput << std::setw(20) << std::left << T_;
			fOutput << std::setw(20) << std::left << P_;
			fOutput << std::setw(20) << std::left << v_;
			fOutput << std::setw(20) << std::left << rho_;
			fOutput << std::setw(20) << std::left << MW_;

			idts_analyzer_.Print(fOutput);

			fOutput << std::endl;
		}
	}


	void PlugFlowReactor::SolveOpenSourceSolvers(const double tf)
	{
		#if OPENSMOKE_USE_DVODE == 1
		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_DVODE)
		{
			// ODE system
			ODESystem_DVODE_PlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_DVODE_PlugFlowReactor::GetInstance();
			odeSystemObject->SetPlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DVODE<ODESystem_DVODE_PlugFlowReactor>   o(odeSystemObject);
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
			ODESystem_DLSODE_PlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_DLSODE_PlugFlowReactor::GetInstance();
			odeSystemObject->SetPlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DLSODE<ODESystem_DLSODE_PlugFlowReactor>   o(odeSystemObject);
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
			ODESystem_DLSODA_PlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_DLSODA_PlugFlowReactor::GetInstance();
			odeSystemObject->SetPlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DLSODA<ODESystem_DLSODA_PlugFlowReactor>   o(odeSystemObject);
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
			ODESystem_DASPK_PlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_DASPK_PlugFlowReactor::GetInstance(NE_);
			odeSystemObject->SetPlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DASPK<ODESystem_DASPK_PlugFlowReactor>   o(odeSystemObject);
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
			ODESystem_RADAU5_PlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_RADAU5_PlugFlowReactor::GetInstance();
			odeSystemObject->SetPlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_RADAU<ODESystem_RADAU5_PlugFlowReactor>   o(odeSystemObject);
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
			ODESystem_CVODE_PlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_CVODE_PlugFlowReactor::GetInstance();
			odeSystemObject->SetPlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_CVODE_Sundials<ODESystem_CVODE_PlugFlowReactor>   o(odeSystemObject);
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
			ODESystem_MEBDF_PlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_MEBDF_PlugFlowReactor::GetInstance();
			odeSystemObject->SetPlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_MEBDF<ODESystem_MEBDF_PlugFlowReactor>   o(odeSystemObject);
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

	void PlugFlowReactor::OpenAllFiles()
	{
		if (plugflow_options_.verbose_output() == true)
		{
			if ( !boost::filesystem::exists(plugflow_options_.output_path()) )
				OpenSMOKE::CreateDirectory(plugflow_options_.output_path());

			if (plugflow_options_.verbose_ascii_file() == true)
				PrepareASCIIFile(plugflow_options_.output_path() / "Output.out");
			
			if (plugflow_options_.verbose_xml_file() == true)
				PrepareXMLFile(plugflow_options_.output_path() / "Output.xml");

			if (on_the_fly_ropa_.is_active() == true)
				PrepareROPAFile(plugflow_options_.output_path() / "ROPA.out");

			if (polimi_soot_analyzer_.is_active() == true)
				PreparePolimiSootFiles(	plugflow_options_.output_path() / "PolimiSoot.out",
										plugflow_options_.output_path() / "PolimiSootDistribution.out");

			// Prepare post-processing output files
			if (on_the_fly_post_processing_.is_active() == true)
				on_the_fly_post_processing_.PrepareOutputFiles();

			// Ignition delay times
			if (idts_analyzer_.is_active())
				idts_analyzer_.Reset();
		}
	}

	void PlugFlowReactor::PrintPolimiSoot(const double t)
	{
		// Integral data
		{
			const unsigned int width = 20;

			fPolimiSoot_ << std::setw(width) << std::left << tau_;
			fPolimiSoot_ << std::setw(width) << std::left << csi_;
			fPolimiSoot_ << std::setw(width) << std::left << T0_;
			fPolimiSoot_ << std::setw(width) << std::left << P0_;
			fPolimiSoot_ << std::setw(width) << std::left << T_;
			fPolimiSoot_ << std::setw(width) << std::left << P_;
			fPolimiSoot_ << std::setw(width) << std::left << v_;
			fPolimiSoot_ << std::setw(width) << std::left << rho_;
			fPolimiSoot_ << std::setw(width) << std::left << MW_;

			polimi_soot_analyzer_.WriteIntegralDataFile(fPolimiSoot_, width);

			fPolimiSoot_ << std::endl;
		}

		// Soot Particle Size Distribution
		if (polimi_soot_analyzer_.write_psdf() == true)
		{
			if (polimi_soot_analyzer_.fv_large() > polimi_soot_analyzer_.threshold_for_psdf())
			{
				polimi_soot_analyzer_.WriteDistribution(fPolimiSootDistribution_, t, csi_, 0., 0., T_);
				fPolimiSootDistribution_ << std::endl;
			}
		}
	}

	void PlugFlowReactor::CloseAllFiles()
	{
		if (plugflow_options_.verbose_output() == true)
		{
			if (plugflow_options_.verbose_ascii_file() == true)
				fASCII_.close();

			if (plugflow_options_.verbose_xml_file() == true)
				CloseXMLFile();

			if (on_the_fly_ropa_.is_active() == true)
				fROPA_.close();

			if (on_the_fly_post_processing_.is_active() == true)
				on_the_fly_post_processing_.CloseOutputFiles();

			if (polimi_soot_analyzer_.is_active() == true)
			{
				fPolimiSoot_.close();

				if (polimi_soot_analyzer_.write_psdf() == true)
					fPolimiSootDistribution_.close();
			}
		}

		if (plugflow_options_.sensitivity_analysis() == true && iXmlMaps_ == true)
			CloseSensitivityXMLFiles();
	}

	void PlugFlowReactor::FinalStatus(const double tf)
	{
		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_);
		Equations(tf, yf_, dummy);

		if (plugflow_options_.verbose_video() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << std::setw(30) << std::left << "Status" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;	
			std::cout << std::setw(30) << std::left << "Time[s]"		<< std::setw(20) << std::left << 0.  << tau_ << std::endl;
			std::cout << std::setw(30) << std::left << "Length[m]"		<< std::setw(20) << std::left << 0.  << csi_ << std::endl;
			std::cout << std::setw(30) << std::left << "T[K]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(3) << T0_ << T_ << std::endl;
			std::cout << std::setw(30) << std::left << "P[atm]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << P0_/101325. << P_/101325. << std::endl;
			std::cout << std::setw(30) << std::left << "v[m/s]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << v0_ << v_ << std::endl;
			std::cout << std::setw(30) << std::left << "rho[kg/m3]"		<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << rho0_ << rho_ << std::endl;
			std::cout << std::setw(30) << std::left << "MW[kg/kmol]"	<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << MW0_ << MW_ << std::endl;
			
			// Equivalence ratio
			{
				const double phi0 = thermodynamicsMap_.GetLocalEquivalenceRatioFromMoleFractions(x0_.GetHandle());
				const double phi = thermodynamicsMap_.GetLocalEquivalenceRatioFromMoleFractions(x_.GetHandle());
				std::cout << std::setw(30) << std::left << "Phi[-]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << phi0 << phi << std::endl;
			}

			// Enthalpy analysis
			{
				const double H_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x_.GetHandle());
				const double U_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(x_.GetHandle());
			
				std::cout << std::setw(30) << std::left << "H[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_/MW0_ << H_/MW_ << std::endl;
				std::cout << std::setw(30) << std::left << "U[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << U0_/MW0_ << U_/MW_  << std::endl;
				std::cout << std::setw(30) << std::left << "E[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_/MW0_+v0_*v0_/2. << H_/MW_+v_*v_/2. << std::endl;
			}

			for(unsigned int i=1;i<=NC_;i++)
				if (omega0_[i] > 0.)
				{
					std::string label = "Conv.(%) " + thermodynamicsMap_.NamesOfSpecies()[i-1];
					std::cout << std::setw(30) << std::left << label << std::setw(20) << std::left << std::fixed << std::setprecision(3) << 0. << (omega0_[i]-omega_[i])/omega0_[i]*100. << std::endl;
				}

			// Atomic analysis
			{
				std::vector<double> sum_initial(thermodynamicsMap_.elements().size());
				std::vector<double> sum_final(thermodynamicsMap_.elements().size());
				for(unsigned int j=0;j<thermodynamicsMap_.elements().size();j++)
				{
					sum_initial[j]	= 0.;
					sum_final[j]	= 0.;
					for(unsigned int i=0;i<NC_;i++)
					{
						sum_initial[j] += thermodynamicsMap_.atomic_composition()(i,j) * x0_[i+1];
						sum_final[j] += thermodynamicsMap_.atomic_composition()(i,j) * x_[i+1];
					}
				}

				const double moles_initial = P0_*v0_/(PhysicalConstants::R_J_kmol*T0_);
				const double moles_final = P_*v_/(PhysicalConstants::R_J_kmol*T_);

				std::cout << std::setw(30) << std::left << "Total species (kmol/m2/s) " << std::setw(20) << std::left << std::scientific << std::setprecision(6) << moles_initial << moles_final << std::endl;

				for(unsigned int j=0;j<thermodynamicsMap_.elements().size();j++)
					if (sum_initial[j] > 0.)
					{
						std::string label = "Element(kmol/m2/s) " + thermodynamicsMap_.elements()[j];
						std::cout << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(6) << sum_initial[j]*moles_initial << sum_final[j]*moles_final << std::endl;
					}
			}

			std::cout << "-----------------------------------------------------------------------------" << std::endl;

			if (plugflow_options_.sensitivity_analysis() == true)
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

	void PlugFlowReactor::FinalSummary(const boost::filesystem::path summary_file, const double tf)
	{
		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_);
		Equations(tf, yf_, dummy);

		if (plugflow_options_.verbose_output() == true)
		{
			std::ofstream fOut;
			fOut.open(summary_file.c_str(), std::ios::out);

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Status" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Time[s]"			<< std::setw(20) << std::left << 0.  << tau_ << std::endl;
			fOut << std::setw(30) << std::left << "Length[m]"		<< std::setw(20) << std::left << 0.  << csi_ << std::endl;
			fOut << std::setw(30) << std::left << "T[K]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(3) << T0_ << T_ << std::endl;
			fOut << std::setw(30) << std::left << "P[atm]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << P0_/101325. << P_/101325. << std::endl;
			fOut << std::setw(30) << std::left << "V[m3]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << v0_ << v_ << std::endl;
			fOut << std::setw(30) << std::left << "rho[kg/m3]"		<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << rho_ << rho_ << std::endl;
			fOut << std::setw(30) << std::left << "MW[kg/kmol]"		<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << MW0_ << MW_ << std::endl;
		
			// Equivalence ratio
			{
				const double phi0 = thermodynamicsMap_.GetLocalEquivalenceRatioFromMoleFractions(x0_.GetHandle());
				const double phi = thermodynamicsMap_.GetLocalEquivalenceRatioFromMoleFractions(x_.GetHandle());
				fOut << std::setw(30) << std::left << "Phi[-]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << phi0 << phi << std::endl;
			}

			// Enthalpy analysis
			{
				const double H_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x_.GetHandle());
				const double U_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(x_.GetHandle());
			
				fOut << std::setw(30) << std::left << "H[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_/MW0_ << H_/MW_ << std::endl;
				fOut << std::setw(30) << std::left << "U[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << U0_/MW0_ << U_/MW_ << std::endl;
				fOut << std::setw(30) << std::left << "E[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_/MW0_+v0_*v0_/2. << H_/MW_+v_*v_/2. << std::endl;
			}

			for(unsigned int i=1;i<=NC_;i++)
				if (omega0_[i] > 0.)
				{
					std::string label = "Conv.(%) " + thermodynamicsMap_.NamesOfSpecies()[i-1];
					fOut << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(3) << 0. << (omega0_[i]-omega_[i])/omega0_[i]*100. << std::endl;
				}

			// Atomic analysis
			{
				std::vector<double> sum_initial(thermodynamicsMap_.elements().size());
				std::vector<double> sum_final(thermodynamicsMap_.elements().size());
				for(unsigned int j=0;j<thermodynamicsMap_.elements().size();j++)
				{
					sum_initial[j]	= 0.;
					sum_final[j]	= 0.;
					for(unsigned int i=0;i<NC_;i++)
					{
						sum_initial[j] += thermodynamicsMap_.atomic_composition()(i,j) * x0_[i+1];
						sum_final[j] += thermodynamicsMap_.atomic_composition()(i,j) * x_[i+1];
					}
				}

				const double moles_initial = P0_*v0_/(PhysicalConstants::R_J_kmol*T0_);
				const double moles_final = P_*v_/(PhysicalConstants::R_J_kmol*T_);

				fOut << std::setw(30) << std::left << "Total species (kmol/m2/s) " << std::setw(20) << std::left << std::scientific << std::setprecision(6) << moles_initial << moles_final << std::endl;

				for(unsigned int j=0;j<thermodynamicsMap_.elements().size();j++)
					if (sum_initial[j] > 0.)
					{
						std::string label = "Element(kmol/m2/s) " + thermodynamicsMap_.elements()[j];
						fOut << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(6) << sum_initial[j]*moles_initial << sum_final[j]*moles_final << std::endl;
					}
			}

			fOut << "-----------------------------------------------------------------------------" << std::endl;

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Mass fractions" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for(unsigned int j=0;j<thermodynamicsMap_.NumberOfSpecies();j++)
				fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << omega0_[j+1] << omega_[j+1] << std::endl;
			fOut << "--------------------------------------------------------------" << std::endl;
		
			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Mole fractions" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for(unsigned int j=0;j<thermodynamicsMap_.NumberOfSpecies();j++)
				fOut << std::setw(30) << std::left << thermodynamicsMap_.NamesOfSpecies()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << x0_[j+1] << x_[j+1] << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;

		}
	}

	void PlugFlowReactor::GetFinalStatus(double& T, double& P, OpenSMOKE::OpenSMOKEVectorDouble& omega)
	{
		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_);
		Equations(0, yf_, dummy);

		T = T_;
		P = P_;
		omega = omega_;
	}

	void PlugFlowReactor::PolimiSootAnalysis(const double t)
	{
		Eigen::VectorXd Yeigen(thermodynamicsMap_.NumberOfSpecies());
		Eigen::VectorXd Xeigen(thermodynamicsMap_.NumberOfSpecies());

		for (unsigned int j = 1; j <= thermodynamicsMap_.NumberOfSpecies(); j++)
		{
			Yeigen(j - 1) = omega_[j];
			Xeigen(j - 1) = x_[j];
		}

		// Analysis of soot
		polimi_soot_analyzer_.Analysis(T_, P_, rho_, Yeigen, Xeigen);

		// Particle size distribution function
		if (polimi_soot_analyzer_.write_psdf() == true)
			polimi_soot_analyzer_.Distribution();

		// Print on file
		PrintPolimiSoot(t);
	}

}

