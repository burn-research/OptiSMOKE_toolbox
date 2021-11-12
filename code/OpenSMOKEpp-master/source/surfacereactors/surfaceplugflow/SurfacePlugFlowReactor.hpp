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
#include "SurfacePlugFlowReactor_OdeInterfaces.h"
#include "SurfacePlugFlowReactor_DaeInterfaces.h"
#include "SurfacePlugFlowReactor_Options.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap.h"
#include "math/external-ode-solvers/ODE_Parameters.h"
#include "math/external-dae-solvers/DAE_Parameters.h"

namespace OpenSMOKE
{
	#if OPENSMOKE_USE_BZZMATH == 1
	
		ODESystem_BzzOde_SurfacePlugFlowReactor* ptOde_SurfacePlugFlowReactor;
		void ODE_Print_SurfacePlugFlowReactor(BzzVector &Y, double t)
		{
			ptOde_SurfacePlugFlowReactor->MyPrint(Y,t);
		}
	
	#endif

	void SurfacePlugFlowReactor::EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap& sensitivityMap, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& plugflow_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options)
	{
		iXmlMaps_ = true;
		
		sensitivityMap_ = &sensitivityMap;

		PrepareSensitivityXMLFiles(thermodynamicsMap, plugflow_options, sensitivity_options);
		
		ChangeDimensions(NE_DAE_, &scaling_Jp, true);
		ChangeDimensions(NE_DAE_, NE_DAE_, &J, true);
	}

	void SurfacePlugFlowReactor::EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap& sensitivityMap, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& plugflow_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options, bool iXmlMaps)
	{
		iXmlMaps_ = iXmlMaps;

		sensitivityMap_ = &sensitivityMap;

		if(iXmlMaps == true)
		    PrepareSensitivityXMLFiles(thermodynamicsMap, plugflow_options, sensitivity_options);
		
		ChangeDimensions(NE_DAE_, &scaling_Jp, true);
		ChangeDimensions(NE_DAE_, NE_DAE_, &J, true);
	}

	void SurfacePlugFlowReactor::MemoryAllocation()
	{
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &RfromGas_, true);
		ChangeDimensions(NC_, &RfromSurface_, true);
		
		ChangeDimensions(NE_DAE_, &yDae0_, true);
		ChangeDimensions(NE_DAE_, &yDaef_, true);
		ChangeDimensions(NE_DAE_, &dyDae0_, true);
		ChangeDimensions(NE_DAE_, &dyDaef_, true);

		ChangeDimensions(NE_ODE_, &yOde0_, true);
		ChangeDimensions(NE_ODE_, &yOdef_, true);

		ChangeDimensions(NC_, &x0_, true);

		ChangeDimensions(SURF_NC_, &Z_, true);
		ChangeDimensions(SURF_NP_, &Gamma_, true);
		ChangeDimensions(SURF_NC_, &Rsurface_, true);
		ChangeDimensions(SURF_NP_, &RsurfacePhases_, true);

		ChangeDimensions(BULK_NC_, &thicknessBulk_, true);
		ChangeDimensions(BULK_NC_, &a_, true);
		ChangeDimensions(BULK_NC_, &Rbulk_, true);

		// TODO: in the current implementation activities are always equal to 1
		a_ = 1.;
	}

	void SurfacePlugFlowReactor::NumericalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEMatrixDouble& J)
	{
		FatalErrorMessage("SurfacePlugFlowReactor::NumericalJacobian not yet implemented");

		// Calculated as suggested by Buzzi (private communication)

		const double ZERO_DER = sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);			
		const double TOLR = 100. * OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT;
		const double TOLA = 1.e-10;
	
		OpenSMOKE::OpenSMOKEVectorDouble y_plus = y;
		OpenSMOKE::OpenSMOKEVectorDouble dy_original(y.Size());
		OpenSMOKE::OpenSMOKEVectorDouble dy_plus(y.Size());

		DaeEquations(t, y, dy_original);

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
			DaeEquations(t, y_plus, dy_plus);

			for(int j=1;j<=y.Size();j++)
				J[j][kd] = (dy_plus[j]-dy_original[j]) * udy;

			y_plus[kd] = y[kd];
		}
	}

	void SurfacePlugFlowReactor::PrepareXMLFile(const boost::filesystem::path output_file_xml, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap)
	{
		fXML_.open(output_file_xml.c_str(), std::ios::out);
		OpenSMOKE::SetXMLFile(fXML_);
		fXML_ << "<Type> SurfacePlugFlowReactor </Type>" << std::endl;

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
		fXML_ << thermodynamicsMap.NumberOfSpecies() << std::endl;
		for (unsigned int j=0;j<NC_;j++)
			fXML_ << thermodynamicsMap.NamesOfSpecies()[j] << " " << thermodynamicsMap.MW(j) << " " << counter++ << std::endl;
		fXML_ << "</mass-fractions>" << std::endl;

		fXML_ << "<surface-fractions>" << std::endl;
		fXML_ << thermodynamicsMap.NumberOfSpecies() << std::endl;
		for (unsigned int j = 0; j<SURF_NC_; j++)
			fXML_ << thermodynamicsMap.vector_names_site_species()[j] << " " << thermodynamicsMap.MW(NC_ + j) << " " << counter++ << std::endl;
		fXML_ << "</surface-fractions>" << std::endl;

		fXML_ << "<profiles>" << std::endl;
	}

	void SurfacePlugFlowReactor::CloseXMLFile()
	{
		fXML_ << "</profiles>" << std::endl;
		fXML_ << "<profiles-size> " << std::endl;
		fXML_ << counter_file_XML_ << " " << 1 + (NC_+1) << std::endl;
		fXML_ << "</profiles-size> " << std::endl;
		fXML_ << "</opensmoke>" << std::endl;
	}

	void SurfacePlugFlowReactor::PrepareSensitivityXMLFiles(OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, SurfacePlugFlowReactor_Options& plugflow_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options)
	{
		FatalErrorMessage("SurfaceBatchReactor::PrepareSensitivityXMLFiles not yet implemented");

		indices_of_sensitivity_species_.resize(sensitivity_options.list_of_species().size());
		for(unsigned int i=0;i<indices_of_sensitivity_species_.size();i++)
			indices_of_sensitivity_species_[i] = thermodynamicsMap.IndexOfSpecies(sensitivity_options.list_of_species()[i]);

		const boost::filesystem::path parent_file = plugflow_options.output_path() / "Sensitivities.xml";
		fSensitivityParentXML_.open(parent_file.c_str(), std::ios::out);
		SetXMLFile(fSensitivityParentXML_);
		fSensitivityParentXML_ << "<variables>" << std::endl;
		
		if (type_ == SURFACEPLUGFLOW_REACTOR_NONISOTHERMAL)
			fSensitivityParentXML_ << indices_of_sensitivity_species_.size()+1 << std::endl;
		else if (type_ == SURFACEPLUGFLOW_REACTOR_ISOTHERMAL)
			fSensitivityParentXML_ << indices_of_sensitivity_species_.size() << std::endl;

		for (unsigned int j=0;j<indices_of_sensitivity_species_.size();j++)
			fSensitivityParentXML_ << thermodynamicsMap.NamesOfSpecies()[indices_of_sensitivity_species_[j]-1] << " " << j << " " << indices_of_sensitivity_species_[j] << std::endl;
		
		if (type_ == SURFACEPLUGFLOW_REACTOR_NONISOTHERMAL)
			fSensitivityParentXML_ << "temperature" << " " << indices_of_sensitivity_species_.size() << " " << NC_+1 << std::endl;
		
		fSensitivityParentXML_ << "</variables>" << std::endl;
		fSensitivityParentXML_ << "<n-parameters> " << std::endl;
		fSensitivityParentXML_ << sensitivityMap_->number_of_parameters() << std::endl;
		fSensitivityParentXML_ << "</n-parameters> " << std::endl;

		fSensitivityChildXML_ = new std::ofstream[indices_of_sensitivity_species_.size()+1];
		for (unsigned int j=0;j<indices_of_sensitivity_species_.size();j++)
		{
			const std::string name = "Sensitivities." + thermodynamicsMap.NamesOfSpecies()[indices_of_sensitivity_species_[j]-1] +".xml";
			const boost::filesystem::path child_file = plugflow_options.output_path() / name;
			fSensitivityChildXML_[j].open(child_file.c_str(), std::ios::out);
		}
		{
			const boost::filesystem::path child_file = plugflow_options.output_path() / "Sensitivities.temperature.xml";
			fSensitivityChildXML_[indices_of_sensitivity_species_.size()].open(child_file.c_str(), std::ios::out);
		}

		for (unsigned int j=0;j<indices_of_sensitivity_species_.size()+1;j++)
		{
			SetXMLFile(fSensitivityChildXML_[j]);
			fSensitivityChildXML_[j] << std::setprecision(5);
			fSensitivityChildXML_[j] << "<coefficients>" << std::endl;
		}		
	}

	void SurfacePlugFlowReactor::CloseSensitivityXMLFiles()
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

	void SurfacePlugFlowReactor::PrepareASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, OpenSMOKE::SurfacePlugFlowReactor_Options& plugflow_options)
	{
		indices_of_output_species_.resize(plugflow_options.output_species().size());
		for (unsigned int i = 0; i<plugflow_options.output_species().size(); i++)
			indices_of_output_species_[i] = thermodynamicsMap.IndexOfSpecies(plugflow_options.output_species()[i]);

		if (indices_of_output_species_.size() != 0)
		{
			widths_of_output_species_.resize(plugflow_options.output_species().size());
			for (unsigned int i = 0; i<plugflow_options.output_species().size(); i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(plugflow_options.output_species()[i], NC_);
		}
		else
		{
			widths_of_output_species_.resize(NC_);
			for (unsigned int i = 0; i<NC_; i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap.NamesOfSpecies()[i], NC_);
		}

		{
			widths_of_output_surface_species_.resize(SURF_NC_);
			for (unsigned int i = 0; i<SURF_NC_; i++)
				widths_of_output_surface_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap.NamesOfSpecies()[NC_ + i], NC_);
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
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "QRGas[W/m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "QRSurf[W/m2]", counter);

		if (BULK_NC_ != 0)
		{
			OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, "thickness[m]", counter);

			for (unsigned int i = 0; i < BULK_NC_; i++)
			{
				const std::string name = thermodynamicsMap.vector_names_bulk_species()[i].substr(0, 11) + "[m]";
				OpenSMOKE::PrintTagOnASCIILabel(20, fASCII_, name, counter);
			}
		}

		for (unsigned int i = 0; i<thermodynamicsMap.number_of_site_phases(0); i++)
		{
			const std::string name = thermodynamicsMap.matrix_names_site_phases()[0][i].substr(0, 13) + "[kmol/m2]";
			OpenSMOKE::PrintTagOnASCIILabel(30, fASCII_, name, counter);
		}

		for (unsigned int i = 0; i<SURF_NC_; i++)
			OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_surface_species_[i], fASCII_, thermodynamicsMap.NamesOfSpecies()[NC_ + i] + "_z", counter);

		if (indices_of_output_species_.size() != 0)
		{
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[indices_of_output_species_[i] - 1] + "_x", counter);
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[indices_of_output_species_[i] - 1] + "_w", counter);
		}
		else
		{
			for (unsigned int i = 0; i<NC_; i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[i] + "_x", counter);
			for (unsigned int i = 0; i<NC_; i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[i] + "_w", counter);
		}

		fOutput << std::endl;
	}

	void SurfacePlugFlowReactor::PrepareBulkASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_bulk_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, OpenSMOKE::SurfacePlugFlowReactor_Options& plugflow_options)
	{
		fBulkASCII_.open(output_file_bulk_ascii.c_str(), std::ios::out);

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
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "QRGas[W/m3]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "QRSurf[W/m2]", counter);

		OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, "total[m]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, "total[kg/m2/s]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, "total[kg/s]", counter);

		for (unsigned int i = 0; i < BULK_NC_; i++)
		{
			const std::string name_species = thermodynamicsMap.vector_names_bulk_species()[i].substr(0, 9);

			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[m]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[kg/m2/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[kg/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[m/s]", counter);
			OpenSMOKE::PrintTagOnASCIILabel(24, fOutput, name_species + "[mm/hr]", counter);
		}

		fOutput << std::endl;
	}

	void SurfacePlugFlowReactor::PrepareParametricASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::SurfacePlugFlowReactor_Options& plugflow_options)
	{
		indices_of_output_species_.resize(plugflow_options.output_species().size());
		for (unsigned int i = 0; i<plugflow_options.output_species().size(); i++)
			indices_of_output_species_[i] = thermodynamicsMap.IndexOfSpecies(plugflow_options.output_species()[i]);

		if (indices_of_output_species_.size() != 0)
		{
			widths_of_output_species_.resize(plugflow_options.output_species().size());
			for (unsigned int i = 0; i<plugflow_options.output_species().size(); i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(plugflow_options.output_species()[i], NC_);
		}
		else
		{
			widths_of_output_species_.resize(NC_);
			for (unsigned int i = 0; i<NC_; i++)
				widths_of_output_species_[i] = CalculateSpeciesFieldWidth(thermodynamicsMap.NamesOfSpecies()[i], NC_);
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
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[indices_of_output_species_[i] - 1] + "_x", counter);
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[indices_of_output_species_[i] - 1] + "_w", counter);
		}
		else
		{
			for (unsigned int i = 0; i<NC_; i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[i] + "_x", counter);
			for (unsigned int i = 0; i<NC_; i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[i] + "_w", counter);
		}

		if (indices_of_output_species_.size() != 0)
		{
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[indices_of_output_species_[i] - 1] + "_x0", counter);
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[indices_of_output_species_[i] - 1] + "_w0", counter);
		}
		else
		{
			for (unsigned int i = 0; i<NC_; i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[i] + "_x0", counter);
			for (unsigned int i = 0; i<NC_; i++)
				OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMap.NamesOfSpecies()[i] + "_w0", counter);
		}

		fOutput << std::endl;
	}

	void SurfacePlugFlowReactor::PrepareASCIIFile(const boost::filesystem::path output_file_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap, OpenSMOKE::SurfacePlugFlowReactor_Options& plugflow_options)
	{
		PrepareASCIIFile(fASCII_, output_file_ascii, thermodynamicsSurfaceMap, plugflow_options);
	}

	void SurfacePlugFlowReactor::PrepareBulkASCIIFile(const boost::filesystem::path output_file_bulk_ascii, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap, OpenSMOKE::SurfacePlugFlowReactor_Options& plugflow_options)
	{
		PrepareBulkASCIIFile(fBulkASCII_, output_file_bulk_ascii, thermodynamicsSurfaceMap, plugflow_options);
	}

	void SurfacePlugFlowReactor::SolveInletConditions(const double tInf, OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap, SurfacePlugFlowReactor_Options& plugflow_options, OpenSMOKE::ODE_Parameters& ode_parameters)
	{
		if (plugflow_options.verbose_video() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << " Solving for inlet conditions to the plug flow reactor...                    " << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
		}

		// Initial conditions
		unsigned int k = 1;
		for (unsigned int i = 1; i <= SURF_NP_; ++i)
			yOde0_[k++] = Gamma0_[i];
		for (unsigned int i = 1; i <= SURF_NC_; ++i)
			yOde0_[k++] = Z0_[i];

		// Print intial conditions
		{
			OpenSMOKE::OpenSMOKEVectorDouble dy0(yOde0_.Size());
			OdeEquations(0., yOde0_, dy0);
			OdePrint(0., yOde0_);
		}

		if (ode_parameters.type() == ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			// Min and max values
			Eigen::VectorXd yMin(NE_ODE_); for (unsigned int i = 0; i<NE_ODE_; i++) yMin(i) = 0.;
			Eigen::VectorXd yMax(NE_ODE_); for (unsigned int i = 0; i<NE_ODE_; i++) yMax(i) = 1.;

			// Initial conditions
			Eigen::VectorXd y0_eigen(yOde0_.Size());
			yOde0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Create the solver
			typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE_SurfacePlugFlowReactor> denseOde;
			typedef OdeSMOKE::MethodGear<denseOde> methodGear;
			OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
			ode_solver.SetReactor(this);

			// Set initial conditions
			ode_solver.SetInitialConditions(0., y0_eigen);

			// Set linear algebra options
			ode_solver.SetLinearAlgebraSolver(ode_parameters.linear_algebra());
			ode_solver.SetFullPivoting(ode_parameters.full_pivoting());

			// Set relative and absolute tolerances
			ode_solver.SetAbsoluteTolerances(ode_parameters.absolute_tolerance());
			ode_solver.SetRelativeTolerances(ode_parameters.relative_tolerance());

			// Set minimum and maximum values
			ode_solver.SetMinimumValues(yMin);
			ode_solver.SetMaximumValues(yMax);

			// Set maximum number of steps
			if (ode_parameters.maximum_number_of_steps() > 0)
				ode_solver.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());

			// Set maximum integration order
			if (ode_parameters.maximum_order() > 0)
				ode_solver.SetMaximumOrder(ode_parameters.maximum_order());

			// Set maximum step size allowed
			if (ode_parameters.maximum_step() > 0)
				ode_solver.SetMaximumStepSize(ode_parameters.maximum_step());

			// Set minimum step size allowed
			if (ode_parameters.minimum_step() > 0)
				ode_solver.SetMinimumStepSize(ode_parameters.minimum_step());

			// Set initial step size
			if (ode_parameters.initial_step() > 0)
				ode_solver.SetFirstStepSize(ode_parameters.initial_step());

			// Solve the system
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			OdeSMOKE::OdeStatus status = ode_solver.Solve(tInf);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			// Check the solution
			if (status > 0)
			{
				ode_solver.Solution(yf_eigen);
				yOdef_.CopyFrom(yf_eigen.data());
				ode_parameters.TransferDataFromOdeSolver(ode_solver, tEnd - tStart);
			}
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else if (ode_parameters.type() == ODE_Parameters::ODE_INTEGRATOR_BZZODE)
		{
			// Min and max values
			BzzVector yMin(NE_ODE_); yMin=0.;
			BzzVector yMax(NE_ODE_); yMax=1.;

			// Initial conditions
			BzzVector y0_bzz(yOde0_.Size());
			yOde0_.CopyTo(y0_bzz.GetHandle());
			
			// Final solution
			BzzVector yf_bzz(y0_bzz.Size());

			ODESystem_BzzOde_SurfacePlugFlowReactor odeplugflow(*this);
			BzzOdeStiffObject o(y0_bzz, 0., &odeplugflow);
			
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
            o.SetTolAbs(ode_parameters.absolute_tolerance());
            o.SetTolRel(ode_parameters.relative_tolerance());

			ptOde_SurfacePlugFlowReactor = &odeplugflow;
			o.StepPrint(ODE_Print_SurfacePlugFlowReactor);

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			yf_bzz = o(tInf,tInf);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			yOdef_.CopyFrom(yf_bzz.GetHandle());

			ode_parameters.SetCPUTime(tEnd-tStart);
			ode_parameters.SetNumberOfFunctionCalls(o.GetNumFunction());
			ode_parameters.SetNumberOfJacobians(o.GetNumNumericalJacobian());
			ode_parameters.SetNumberOfFactorizations(o.GetNumFactorization());
			ode_parameters.SetNumberOfSteps(o.GetNumStep());
			ode_parameters.SetLastOrderUsed(o.GetOrderUsed());
			ode_parameters.SetLastStepUsed(o.GetHUsed());
		}
		#endif
		else 
		{
			OdeSolveOpenSourceSolvers(tInf, ode_parameters);
		}

		if (plugflow_options.verbose_video() == true)
			ode_parameters.Status(std::cout);

		OdeFinalSummary(plugflow_options.output_path() / "InletSummary.out", tInf, thermodynamicsSurfaceMap, plugflow_options);

		SetInletConditions();

		// Reset counters
		counter_file_video_ = 0;
		counter_file_ASCII_ = 0;
		counter_file_XML_ = 0;
		counter_sensitivity_XML_ = 0; 
	}

	void SurfacePlugFlowReactor::OdeSolveOpenSourceSolvers(const double tf, OpenSMOKE::ODE_Parameters& ode_parameters)
	{
		#if OPENSMOKE_USE_DVODE == 1
		if (ode_parameters.type() == ODE_Parameters::ODE_INTEGRATOR_DVODE)
		{
			// ODE system
			ODESystem_DVODE_SurfacePlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_DVODE_SurfacePlugFlowReactor::GetInstance();
			odeSystemObject->SetSurfacePlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DVODE<ODESystem_DVODE_SurfacePlugFlowReactor>   o(odeSystemObject);
			o.SetDimensions(NE_ODE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., yOde0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yOdef_.GetHandle());

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
			ODESystem_DLSODE_SurfacePlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_DLSODE_SurfacePlugFlowReactor::GetInstance();
			odeSystemObject->SetSurfacePlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DLSODE<ODESystem_DLSODE_SurfacePlugFlowReactor>   o(odeSystemObject);
			o.SetDimensions(NE_ODE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., yOde0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yOdef_.GetHandle());

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
			ODESystem_DLSODA_SurfacePlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_DLSODA_SurfacePlugFlowReactor::GetInstance();
			odeSystemObject->SetSurfacePlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DLSODA<ODESystem_DLSODA_SurfacePlugFlowReactor>   o(odeSystemObject);
			o.SetDimensions(NE_ODE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., yOde0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yOdef_.GetHandle());

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
			ODESystem_DASPK_SurfacePlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_DASPK_SurfacePlugFlowReactor::GetInstance(NE_ODE_);
			odeSystemObject->SetSurfacePlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_DASPK<ODESystem_DASPK_SurfacePlugFlowReactor>   o(odeSystemObject);
			o.SetDimensions(NE_ODE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., yOde0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yOdef_.GetHandle());

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
			ODESystem_RADAU5_SurfacePlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_RADAU5_SurfacePlugFlowReactor::GetInstance();
			odeSystemObject->SetSurfacePlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_RADAU<ODESystem_RADAU5_SurfacePlugFlowReactor>   o(odeSystemObject);
			o.SetDimensions(NE_ODE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., yOde0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yOdef_.GetHandle());

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
			ODESystem_CVODE_SurfacePlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_CVODE_SurfacePlugFlowReactor::GetInstance();
			odeSystemObject->SetSurfacePlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_CVODE_Sundials<ODESystem_CVODE_SurfacePlugFlowReactor>   o(odeSystemObject);
			o.SetDimensions(NE_ODE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., yOde0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yOdef_.GetHandle());

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
			ODESystem_MEBDF_SurfacePlugFlowReactor *odeSystemObject;
			odeSystemObject = ODESystem_MEBDF_SurfacePlugFlowReactor::GetInstance();
			odeSystemObject->SetSurfacePlugFlowReactor(this);

			// Ode solver
			OpenSMOKE::OpenSMOKE_MEBDF<ODESystem_MEBDF_SurfacePlugFlowReactor>   o(odeSystemObject);
			o.SetDimensions(NE_ODE_);
			o.SetAbsoluteTolerance(ode_parameters.absolute_tolerance());
			o.SetRelativeTolerance(ode_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(ode_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., yOde0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yOdef_.GetHandle());

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
	
	void SurfacePlugFlowReactor::DaeSolveOpenSourceSolvers(const double tf, const bool* algebraic_equations, OpenSMOKE::DAE_Parameters& dae_parameters)
	{
		#if OPENSMOKE_USE_SUNDIALS == 1
		if (dae_parameters.type() == DAE_Parameters::DAE_INTEGRATOR_IDA)
		{
			// DAE system
			DAESystem_IDA_SurfacePlugFlowReactor *daeSystemObject;
			daeSystemObject = DAESystem_IDA_SurfacePlugFlowReactor::GetInstance(NE_DAE_, algebraic_equations);
			daeSystemObject->SetSurfacePlugFlowReactor(this);

			// Daesolver
			OpenSMOKE::OpenSMOKE_IDA_Sundials<DAESystem_IDA_SurfacePlugFlowReactor>   o(daeSystemObject);
			o.SetDimensions(NE_DAE_);
			o.SetAbsoluteTolerance(dae_parameters.absolute_tolerance());
			o.SetRelativeTolerance(dae_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(dae_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., yDae0_.GetHandle(), dyDae0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf, algebraic_equations);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yDaef_.GetHandle(), dyDaef_.GetHandle());

			dae_parameters.SetCPUTime(tEnd-tStart);
			dae_parameters.SetNumberOfSteps(o.GetNumberOfSteps());
			dae_parameters.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			dae_parameters.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			dae_parameters.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			dae_parameters.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			dae_parameters.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			dae_parameters.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			dae_parameters.SetLastOrderUsed(o.GetLastOrderUsed());
			dae_parameters.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif

		#if OPENSMOKE_USE_DASPK == 1
		if (dae_parameters.type() == DAE_Parameters::DAE_INTEGRATOR_DASPK)
		{
			// DAE system
			DAESystem_DASPK_SurfacePlugFlowReactor *daeSystemObject;
			daeSystemObject = DAESystem_DASPK_SurfacePlugFlowReactor::GetInstance(NE_DAE_, algebraic_equations);
			daeSystemObject->SetSurfacePlugFlowReactor(this);

			// Dae solver
			OpenSMOKE::OpenSMOKE_DASPK_DAE<DAESystem_DASPK_SurfacePlugFlowReactor>   o(daeSystemObject);
			o.SetDimensions(NE_DAE_);
			o.SetAbsoluteTolerance(dae_parameters.absolute_tolerance());
			o.SetRelativeTolerance(dae_parameters.relative_tolerance());
			o.SetMaximumNumberOfSteps(dae_parameters.maximum_number_of_steps());
			o.SetAnalyticalJacobian(false);
			o.SetInitialValues(0., yDae0_.GetHandle(), dyDae0_.GetHandle());

			const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solve(tf, algebraic_equations);
			const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			o.Solution(yDaef_.GetHandle(), dyDaef_.GetHandle());

			dae_parameters.SetCPUTime(tEnd - tStart);
			dae_parameters.SetNumberOfSteps(o.GetNumberOfSteps());
			dae_parameters.SetNumberOfFunctionCalls(o.GetNumberOfFunctionEvaluations());
			dae_parameters.SetNumberOfJacobians(o.GetNumberOfJacobianEvaluations());
			dae_parameters.SetNumberOfFactorizations(o.GetNumberOfLUFactorizations());
			dae_parameters.SetNumberOfNonLinearIterations(o.GetNumberOfNonLinearIterations());
			dae_parameters.SetNumberOfConvergenceFailures(o.GetNumberOfConvergenceFailures());
			dae_parameters.SetNumberOfErrorTestFailures(o.GetNumberOfErrorTestFailures());
			dae_parameters.SetLastOrderUsed(o.GetLastOrderUsed());
			dae_parameters.SetLastStepUsed(o.GetLastStepUsed());
		}
		#endif
	}

	void SurfacePlugFlowReactor::OpenAllFiles(OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsMap, OpenSMOKE::SurfacePlugFlowReactor_Options& plugflow_options)
	{
		if (plugflow_options.verbose_output() == true)
		{
			if ( !boost::filesystem::exists(plugflow_options.output_path()) )
				OpenSMOKE::CreateDirectory(plugflow_options.output_path());

			if (plugflow_options.verbose_ascii_file() == true)
				PrepareASCIIFile(plugflow_options.output_path() / "Output.out", thermodynamicsMap, plugflow_options);
			
			if (plugflow_options.verbose_xml_file() == true)
				PrepareXMLFile(plugflow_options.output_path() / "Output.xml", thermodynamicsMap);

			if (plugflow_options.verbose_ascii_file() == true && BULK_NC_ != 0)
				PrepareBulkASCIIFile(plugflow_options.output_path() / "Output.bulk.out", thermodynamicsMap, plugflow_options);
		}
	}

	void SurfacePlugFlowReactor::CloseAllFiles(OpenSMOKE::SurfacePlugFlowReactor_Options& plugflow_options)
	{
		if (plugflow_options.verbose_output() == true)
			{
				if (plugflow_options.verbose_ascii_file() == true)
					fASCII_.close();

				if (plugflow_options.verbose_xml_file() == true)
					CloseXMLFile();
			}

		if (plugflow_options.sensitivity_analysis() == true && iXmlMaps_ == true)
			CloseSensitivityXMLFiles();
	}

	void SurfacePlugFlowReactor::DaeFinalStatus(const double tf,	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
																OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap, 
																SurfacePlugFlowReactor_Options& plugflow_options)
	{
		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_DAE_);
		DaeEquations(tf, yDaef_, dummy);

		if (plugflow_options.verbose_video() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << std::setw(30) << std::left << "Status" << std::setw(20) << std::left << "Inlet" << "Outlet" << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;	
			std::cout << std::setw(30) << std::left << "Time[s]"		<< std::setw(20) << std::left << 0.  << tau_ << std::endl;
			std::cout << std::setw(30) << std::left << "Length[m]"		<< std::setw(20) << std::left << 0.  << csi_ << std::endl;
			std::cout << std::setw(30) << std::left << "T[K]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(3) << T0_ << T_ << std::endl;
			std::cout << std::setw(30) << std::left << "P[atm]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << P0_/101325. << P_/101325. << std::endl;
			std::cout << std::setw(30) << std::left << "v[m/s]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << v0_ << v_ << std::endl;
			std::cout << std::setw(30) << std::left << "rho[kg/m3]"		<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << rho0_ << rho_ << std::endl;
			std::cout << std::setw(30) << std::left << "MW[kg/kmol]"	<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << MW0_ << MW_ << std::endl;
			std::cout << std::setw(30) << std::left << "Gamma[kmol/m2]" << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Gamma0_.SumElements() << Gamma0_.SumElements() << std::endl;
			std::cout << std::setw(30) << std::left << "Bulk[m]"		<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << thicknessBulk0_ << thicknessBulk_.SumElements() << std::endl;

			// Equivalence ratio
			{
				const double phi0 = thermodynamicsMap.GetLocalEquivalenceRatioFromMoleFractions(x0_.GetHandle());
				const double phi = thermodynamicsMap.GetLocalEquivalenceRatioFromMoleFractions(x_.GetHandle());
				std::cout << std::setw(30) << std::left << "Phi[-]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << phi0 << phi << std::endl;
			}

			// Enthalpy analysis
			{
				const double H_ = thermodynamicsMap.hMolar_Mixture_From_MoleFractions(x_.GetHandle());
				const double U_ = thermodynamicsMap.uMolar_Mixture_From_MoleFractions(x_.GetHandle());
			
				std::cout << std::setw(30) << std::left << "H[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_/MW0_ << H_/MW_ << std::endl;
				std::cout << std::setw(30) << std::left << "U[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << U0_/MW0_ << U_/MW_  << std::endl;
				std::cout << std::setw(30) << std::left << "E[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_/MW0_+v0_*v0_/2. << H_/MW_+v_*v_/2. << std::endl;
			}

			for(unsigned int i=1;i<=NC_;i++)
				if (omega0_[i] > 0.)
				{
					std::string label = "Conv.(%) " + thermodynamicsMap.NamesOfSpecies()[i-1];
					std::cout << std::setw(30) << std::left << label << std::setw(20) << std::left << std::fixed << std::setprecision(3) << 0. << (omega0_[i]-omega_[i])/omega0_[i]*100. << std::endl;
				}

			// Atomic analysis
			{
				std::vector<double> sum_initial(thermodynamicsMap.elements().size());
				std::vector<double> sum_final(thermodynamicsMap.elements().size());
				for(unsigned int j=0;j<thermodynamicsMap.elements().size();j++)
				{
					sum_initial[j]	= 0.;
					sum_final[j]	= 0.;
					for(unsigned int i=0;i<NC_;i++)
					{
						sum_initial[j] += thermodynamicsMap.atomic_composition()(i,j) * x0_[i+1];
						sum_final[j] += thermodynamicsMap.atomic_composition()(i,j) * x_[i+1];
					}
				}

				const double moles_initial = P0_*v0_/(PhysicalConstants::R_J_kmol*T0_);
				const double moles_final = P_*v_/(PhysicalConstants::R_J_kmol*T_);

				std::cout << std::setw(30) << std::left << "Total species (kmol/m2/s) " << std::setw(20) << std::left << std::scientific << std::setprecision(6) << moles_initial << moles_final << std::endl;

				for(unsigned int j=0;j<thermodynamicsMap.elements().size();j++)
					if (sum_initial[j] > 0.)
					{
						std::string label = "Element(kmol/m2/s) " + thermodynamicsMap.elements()[j];
						std::cout << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(6) << sum_initial[j]*moles_initial << sum_final[j]*moles_final << std::endl;
					}
			}

			std::cout << "-----------------------------------------------------------------------------" << std::endl;

			if (plugflow_options.sensitivity_analysis() == true)
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

	void SurfacePlugFlowReactor::OdeFinalSummary(	const boost::filesystem::path summary_file, const double tf,
													OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap,
													SurfacePlugFlowReactor_Options& plugflow_options)
	{
		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_ODE_);
		OdeEquations(tf, yOdef_, dummy);

		if (plugflow_options.verbose_output() == true)
		{
			std::ofstream fOut;
			fOut.open(summary_file.c_str(), std::ios::out);

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Bulk species [m]" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for (unsigned int j = 0; j<BULK_NC_; j++)
				fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap.vector_names_bulk_species()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << thicknessBulk0_ << thicknessBulk_[j + 1] << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;


			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Site densities [kmol/m2]" << std::setw(20) << std::left << "First guess" << "Calculated" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for (unsigned int j = 0; j<thermodynamicsSurfaceMap.number_of_site_phases(0); j++)
				fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap.matrix_names_site_phases()[0][j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Gamma0_[j+1] << Gamma_[j + 1] << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;

			for (unsigned int i = 0; i<thermodynamicsSurfaceMap.number_of_site_phases(0); i++)
			{
				fOut << std::endl;
				fOut << "-----------------------------------------------------------------------------" << std::endl;
				fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap.matrix_names_site_phases()[0][i] << std::endl;
				fOut << std::setw(30) << std::left << "Surface fractions" << std::setw(20) << std::left << "First guess" << "Calculated" << std::endl;
				fOut << "-----------------------------------------------------------------------------" << std::endl;
				for (unsigned int j = 0; j<thermodynamicsSurfaceMap.number_of_site_species(); j++)
					if (thermodynamicsSurfaceMap.vector_site_phases_belonging()[j] == i) 
					fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap.NamesOfSpecies()[j + NC_] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Z0_[j+1] << Z_[j + 1] << std::endl;
				fOut << "-----------------------------------------------------------------------------" << std::endl;
			}
		}
	}

	void SurfacePlugFlowReactor::DaeFinalSummary(	const boost::filesystem::path summary_file, const double tf, 
													OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
													OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap, 
													SurfacePlugFlowReactor_Options& plugflow_options)
	{
		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_DAE_);
		DaeEquations(tf, yDaef_, dummy);

		if (plugflow_options.verbose_output() == true)
		{
			std::ofstream fOut;
			fOut.open(summary_file.c_str(), std::ios::out);

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Status" << std::setw(20) << std::left << "Inlet" << "Outlet" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Time[s]"			<< std::setw(20) << std::left << 0.  << tau_ << std::endl;
			fOut << std::setw(30) << std::left << "Length[m]"		<< std::setw(20) << std::left << 0.  << csi_ << std::endl;
			fOut << std::setw(30) << std::left << "T[K]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(3) << T0_ << T_ << std::endl;
			fOut << std::setw(30) << std::left << "P[atm]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << P0_/101325. << P_/101325. << std::endl;
			fOut << std::setw(30) << std::left << "v[m/s]"			<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << v0_ << v_ << std::endl;
			fOut << std::setw(30) << std::left << "G[kg/m2/s]"		<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << G0_ << G_ << std::endl;
			fOut << std::setw(30) << std::left << "rho[kg/m3]"		<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << rho0_ << rho_ << std::endl;
			fOut << std::setw(30) << std::left << "MW[kg/kmol]"		<< std::setw(20) << std::left << std::fixed << std::setprecision(5) << MW0_ << MW_ << std::endl;
			fOut << std::setw(30) << std::left << "Gamma[kmol/m2]"  << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Gamma0_.SumElements() << Gamma0_.SumElements() << std::endl;
			fOut << std::setw(30) << std::left << "Bulk[m]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << thicknessBulk0_ << thicknessBulk_.SumElements() << std::endl;

			// Equivalence ratio
			{
				const double phi0 = thermodynamicsMap.GetLocalEquivalenceRatioFromMoleFractions(x0_.GetHandle());
				const double phi = thermodynamicsMap.GetLocalEquivalenceRatioFromMoleFractions(x_.GetHandle());
				fOut << std::setw(30) << std::left << "Phi[-]" << std::setw(20) << std::left << std::fixed << std::setprecision(5) << phi0 << phi << std::endl;
			}

			// Enthalpy analysis
			{
				const double H_ = thermodynamicsMap.hMolar_Mixture_From_MoleFractions(x_.GetHandle());
				const double U_ = thermodynamicsMap.uMolar_Mixture_From_MoleFractions(x_.GetHandle());
			
				fOut << std::setw(30) << std::left << "H[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_/MW0_ << H_/MW_ << std::endl;
				fOut << std::setw(30) << std::left << "U[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << U0_/MW0_ << U_/MW_ << std::endl;
				fOut << std::setw(30) << std::left << "E[J/kg]"			<< std::setw(20) << std::left << std::scientific << std::setprecision(6) << H0_/MW0_+v0_*v0_/2. << H_/MW_+v_*v_/2. << std::endl;
			}

			for(unsigned int i=1;i<=NC_;i++)
				if (omega0_[i] > 0.)
				{
					std::string label = "Conv.(%) " + thermodynamicsMap.NamesOfSpecies()[i-1];
					fOut << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(3) << 0. << (omega0_[i]-omega_[i])/omega0_[i]*100. << std::endl;
				}

			// Atomic analysis
			{
				std::vector<double> sum_initial(thermodynamicsMap.elements().size());
				std::vector<double> sum_final(thermodynamicsMap.elements().size());
				for(unsigned int j=0;j<thermodynamicsMap.elements().size();j++)
				{
					sum_initial[j]	= 0.;
					sum_final[j]	= 0.;
					for(unsigned int i=0;i<NC_;i++)
					{
						sum_initial[j] += thermodynamicsMap.atomic_composition()(i,j) * x0_[i+1];
						sum_final[j] += thermodynamicsMap.atomic_composition()(i,j) * x_[i+1];
					}
				}

				const double moles_initial = P0_*v0_/(PhysicalConstants::R_J_kmol*T0_);
				const double moles_final = P_*v_/(PhysicalConstants::R_J_kmol*T_);

				fOut << std::setw(30) << std::left << "Total species (kmol/m2/s) " << std::setw(20) << std::left << std::scientific << std::setprecision(6) << moles_initial << moles_final << std::endl;

				for(unsigned int j=0;j<thermodynamicsMap.elements().size();j++)
					if (sum_initial[j] > 0.)
					{
						std::string label = "Element(kmol/m2/s) " + thermodynamicsMap.elements()[j];
						fOut << std::setw(30) << std::left << label << std::setw(20) << std::left << std::scientific << std::setprecision(6) << sum_initial[j]*moles_initial << sum_final[j]*moles_final << std::endl;
					}
			}

			fOut << "-----------------------------------------------------------------------------" << std::endl;

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Bulk species [m]" << std::setw(20) << std::left << "Initial" << "Final" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for (unsigned int j = 0; j<BULK_NC_; j++)
				fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap.vector_names_bulk_species()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << thicknessBulk0_ << thicknessBulk_[j + 1] << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Site densities [kmol/m2]" << std::setw(20) << std::left << "Inlet" << "Outlet" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for (unsigned int j = 0; j<thermodynamicsSurfaceMap.number_of_site_phases(0); j++)
				fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap.matrix_names_site_phases()[0][j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Gamma0_[j + 1] << Gamma_[j + 1] << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;

			for (unsigned int i = 0; i<thermodynamicsSurfaceMap.number_of_site_phases(0); i++)
			{
				fOut << std::endl;
				fOut << "-----------------------------------------------------------------------------" << std::endl;
				fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap.matrix_names_site_phases()[0][i] << std::endl;
				fOut << std::setw(30) << std::left << "Surface fractions" << std::setw(20) << std::left << "Inlet" << "Outlet" << std::endl;
				fOut << "-----------------------------------------------------------------------------" << std::endl;
				for (unsigned int j = 0; j<thermodynamicsSurfaceMap.number_of_site_species(); j++)
					if (thermodynamicsSurfaceMap.vector_site_phases_belonging()[j] == i) 
						fOut << std::setw(30) << std::left << thermodynamicsSurfaceMap.NamesOfSpecies()[j + NC_] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << Z0_[j + 1] << Z_[j + 1] << std::endl;
				fOut << "-----------------------------------------------------------------------------" << std::endl;
			}

			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Mass fractions" << std::setw(20) << std::left << "Inlet" << "Outlet" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for(unsigned int j=0;j<thermodynamicsMap.NumberOfSpecies();j++)
				fOut << std::setw(30) << std::left << thermodynamicsMap.NamesOfSpecies()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << omega0_[j+1] << omega_[j+1] << std::endl;
			fOut << "--------------------------------------------------------------" << std::endl;
		
			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(30) << std::left << "Mole fractions" << std::setw(20) << std::left << "Inlet" << "Outlet" << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;
			for(unsigned int j=0;j<thermodynamicsMap.NumberOfSpecies();j++)
				fOut << std::setw(30) << std::left << thermodynamicsMap.NamesOfSpecies()[j] << std::setw(20) << std::left << std::scientific << std::setprecision(6) << x0_[j+1] << x_[j+1] << std::endl;
			fOut << "-----------------------------------------------------------------------------" << std::endl;

		}
	}

	void SurfacePlugFlowReactor::DaeGetFinalStatus(double& T, double& P, double& G, OpenSMOKE::OpenSMOKEVectorDouble& omega, OpenSMOKE::OpenSMOKEVectorDouble& Z, OpenSMOKE::OpenSMOKEVectorDouble& Gamma, OpenSMOKE::OpenSMOKEVectorDouble& thicknessBulk)
	{
		OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_DAE_);
		DaeEquations(0, yDaef_, dummy);

		T = T_;
		P = P_;
		omega = omega_;
		G = G_;
		Z = Z_;
		Gamma = Gamma_;
		thicknessBulk = thicknessBulk_;
	}

}

