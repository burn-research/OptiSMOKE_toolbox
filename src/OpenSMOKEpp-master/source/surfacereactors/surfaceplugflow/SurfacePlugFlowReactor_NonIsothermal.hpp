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

#include "SurfacePlugFlowReactor_OdeInterfaces.h"
#include "math/PhysicalConstants.h"
#include "math/native-ode-solvers/MultiValueSolver"

namespace OpenSMOKE
{
	#if OPENSMOKE_USE_BZZMATH == 1

		DAESystem_BzzDae_SurfacePlugFlowReactor* ptDae_SurfacePlugFlowReactor_NonIsothermal;
		void DAE_Print_SurfacePlugFlowReactor_NonIsothermal(BzzVector &Y, double t)
		{
			ptDae_SurfacePlugFlowReactor_NonIsothermal->MyPrint(Y, t);
		}
	
	#endif


	SurfacePlugFlowReactor_NonIsothermal::SurfacePlugFlowReactor_NonIsothermal(	
								OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
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
                                const double global_thermal_exchange_coefficient,
                                const double cross_section_over_perimeter,
                                const double T_environment) :
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap), 
	thermodynamicsSurfaceMap_(thermodynamicsSurfaceMap), 
	kineticsSurfaceMap_(kineticsSurfaceMap), 
	ode_parameters_(ode_parameters),
	dae_parameters_(dae_parameters),	
	plugflow_options_(plugflow_options)
	{
		type_ = SURFACEPLUGFLOW_REACTOR_NONISOTHERMAL;

		artificial_speedup_coefficient_ = 1.e4;
		time_independent_ = time_independent;
		constant_pressure_ = constant_pressure;
		ode_iteration_ = 0;
		dae_iteration_ = 0;
		counter_file_video_ = 0;
		counter_file_ASCII_ = 0;
		counter_file_XML_ = 0;
		counter_sensitivity_XML_ = 0;

		v0_ = v0;
		T0_ = T0 ;
		P0_ = P0;
		omega0_ = omega0;
		Z0_ = Z0;
		Gamma0_ = Gamma0;
		thicknessBulk0_ = 0.;

		cross_section_over_perimeter_ = cross_section_over_perimeter;
        global_thermal_exchange_coefficient_ = global_thermal_exchange_coefficient;
        cross_section_over_perimeter_ = cross_section_over_perimeter;
        T_environment_ = T_environment;

		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();
		SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
		SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
		SURF_NR_ = kineticsSurfaceMap_.NumberOfReactions();

		// Bulk phases
		BULK_NP_ = thermodynamicsSurfaceMap_.number_of_bulk_phases(0);
		BULK_NC_ = thermodynamicsSurfaceMap_.number_of_bulk_species();

		NE_ODE_ = SURF_NC_ + SURF_NP_ + BULK_NC_;
		NE_DAE_ = NC_ + SURF_NC_ + SURF_NP_ + 3 + BULK_NC_;

		MemoryAllocation();

		MW0_ = thermodynamicsMap_.MolecularWeight_From_MassFractions(omega0_.GetHandle());
		rho0_ = P0_ * MW0_ / (PhysicalConstants::R_J_kmol*T0_);
		G0_ = rho0_ * v0_;

		site_non_conservation_ = site_non_conservation;

		T_ = T0_;
		P_ = P0_;
		v_ = v0_;
		omega_ = omega0_;
		Z_ = Z0_;
		Gamma_ = Gamma0_;
		GammaFromEqn_ = Gamma0_;
		MW_ = MW0_;
		QRGas_ = 0.;
		QRSurface_ = 0.;
		rho_ = rho0_;
		G_ = G0_;
		thicknessBulk_ = thicknessBulk0_;

		// Inlet enthalpy and internal energy
		thermodynamicsMap_.MoleFractions_From_MassFractions(x0_.GetHandle(), MW0_, omega0_.GetHandle());
		thermodynamicsMap_.SetPressure(P0_);
		thermodynamicsMap_.SetTemperature(T0_);
		H0_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x0_.GetHandle());
		U0_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(x0_.GetHandle());

		// Inlet enthalpy and internal energy for surfaces
		thermodynamicsSurfaceMap_.SetPressure(P0_);
		thermodynamicsSurfaceMap_.SetTemperature(T0_);
		H0_Surface_over_Volume_ = 0.;
		U0_Surface_over_Volume_ = 0.;
		for (unsigned int i = 1; i <= SURF_NC_; i++)
		{
			const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i - 1] + 1;
			H0_Surface_over_Volume_ += Z0_[i] * thermodynamicsSurfaceMap_.Species_H_over_RT()[i + NC_-1] * PhysicalConstants::R_J_kmol*T0_ *
				Gamma0_[index_phase] / cross_section_over_perimeter_;
			U0_Surface_over_Volume_ += Z0_[i] * (thermodynamicsSurfaceMap_.Species_H_over_RT()[i + NC_-1] - 1.)*PhysicalConstants::R_J_kmol*T0_ *
				Gamma0_[index_phase] / cross_section_over_perimeter_;
		}

		OpenAllFiles(thermodynamicsSurfaceMap_, plugflow_options_);
	}

	int SurfacePlugFlowReactor_NonIsothermal::OdeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
	{
		// Recover unknowns
		unsigned int k=1;
		#ifdef CHECK_MASSFRACTIONS
		for(unsigned int i=1;i<=SURF_NP_;++i)
			GammaFromEqn_[i] = std::max(y[k++], 0.);
		for(unsigned int i=1;i<=SURF_NC_;++i)
			Z_[i] = std::max(y[k++], 0.);
		for (unsigned int i = 1; i <= BULK_NC_; ++i)
			thicknessBulk_[i] = y[k++];
		#else
		for(unsigned int i=1;i<=SURF_NP_;++i)
			GammaFromEqn_[i] = y[k++];
		for(unsigned int i=1;i<=SURF_NC_;++i)
			Z_[i] = y[k++];
		for (unsigned int i = 1; i <= BULK_NC_; ++i)
			thicknessBulk_[i] = y[k++];
		#endif

		for (unsigned int i = 1; i <= SURF_NP_; ++i)
			if (site_non_conservation_[i - 1] == true)
				Gamma_[i] = GammaFromEqn_[i];

		// Calculates the volume and the concentrations of species
		cTot_ = P0_ / (PhysicalConstants::R_J_kmol * T0_);
		Product(cTot_, x0_, &c_);

		// Calculates heterogeneous kinetics
		{
			thermodynamicsSurfaceMap_.SetPressure(P0_);
			thermodynamicsSurfaceMap_.SetTemperature(T0_);
			kineticsSurfaceMap_.SetPressure(P0_);
			kineticsSurfaceMap_.SetTemperature(T0_);
			kineticsSurfaceMap_.ReactionEnthalpiesAndEntropies();
			kineticsSurfaceMap_.KineticConstants();
			kineticsSurfaceMap_.ReactionRates(c_.GetHandle(), Z_.GetHandle(), a_.GetHandle(), Gamma_.GetHandle());
			kineticsSurfaceMap_.FormationRates(RfromSurface_.GetHandle(), Rsurface_.GetHandle(), Rbulk_.GetHandle(), RsurfacePhases_.GetHandle());
		}

		// Recovering residuals
		{
			unsigned int k = 1;

			// Phases
			for (unsigned int i = 1; i <= SURF_NP_; ++i)
				dy[k++] = RsurfacePhases_[i];

			// Heterogeneous species
			for (unsigned int i = 1; i <= SURF_NC_; ++i)
			{
				const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i - 1] + 1;
				dy[k++] = (thermodynamicsSurfaceMap_.vector_occupancies_site_species()[i - 1] * Rsurface_[i] -
								Z_[i] * dy[index_phase]) / Gamma_[index_phase];
			}

			// Bulk species
			for (unsigned int i = 1; i <= BULK_NC_; ++i)
				dy[k++] = 0.;
		}

		return 0;
	}
	
	int SurfacePlugFlowReactor_NonIsothermal::DaeEquations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
	{
		unsigned int k=1;
		#ifdef CHECK_MASSFRACTIONS
		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = std::max(y[k++], 0.);
		G_ = std::max(y[k++], 0.);
		for(unsigned int i=1;i<=SURF_NP_;++i)
			Gamma_[i] = std::max(y[k++], 0.);
		for(unsigned int i=1;i<=SURF_NC_;++i)
			Z_[i] = std::max(y[k++], 0.);
		for (unsigned int i = 1; i <= BULK_NC_; ++i)
			thicknessBulk_[i] = y[k++];
		T_ = std::max(y[k++], 0.);
		#else
		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = y[k++];
		G_ = y[k++];
		for(unsigned int i=1;i<=SURF_NP_;++i)
			Gamma_[i] = y[k++];
		for(unsigned int i=1;i<=SURF_NC_;++i)
			Z_[i] = y[k++];
		for (unsigned int i = 1; i <= BULK_NC_; ++i)
			thicknessBulk_[i] = y[k++];
		T_ = y[k++];
		#endif

		if (time_independent_ == false)
		{
			csi_ = t;
			tau_ = y[k++];
		}
		else if (time_independent_ == true)
		{
			tau_ = t;
			csi_ = y[k++];
		}

		// Calculates the volume and the concentrations of species
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
		if (constant_pressure_ == true)
		{
			cTot_ = P_ / (PhysicalConstants::R_J_kmol * T_);
			Product(cTot_, x_, &c_);
			rho_ = cTot_*MW_;
			v_ = G_ / rho_;
		}
		else
		{
			P_ = G_ / v_ * PhysicalConstants::R_J_kmol * T_ / MW_;
			cTot_ = P_ / (PhysicalConstants::R_J_kmol * T_);
			Product(cTot_, x_, &c_);
			rho_ = cTot_*MW_;
		}

		// Calculates thermodynamic properties
		thermodynamicsMap_.SetTemperature(T_);
		thermodynamicsMap_.SetPressure(P_);
		
		const double CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
		CvMixMass_ = (CpMixMolar - PhysicalConstants::R_J_kmol) / MW_;
		CpMixMass_ = CpMixMolar / MW_;
		
		// Calculates kinetics
		{
			kineticsMap_.SetTemperature(T_);
			kineticsMap_.SetPressure(P_);

			kineticsMap_.KineticConstants();
			kineticsMap_.ReactionRates(c_.GetHandle());
			kineticsMap_.FormationRates(RfromGas_.GetHandle());
			QRGas_ = kineticsMap_.HeatRelease(RfromGas_.GetHandle());
		}

		// Calculates heterogeneous kinetics
		{
			thermodynamicsSurfaceMap_.SetPressure(P_);
			thermodynamicsSurfaceMap_.SetTemperature(T_);
			kineticsSurfaceMap_.SetPressure(P_);
			kineticsSurfaceMap_.SetTemperature(T_);
			kineticsSurfaceMap_.ReactionEnthalpiesAndEntropies();
			kineticsSurfaceMap_.KineticConstants();
			kineticsSurfaceMap_.ReactionRates(c_.GetHandle(), Z_.GetHandle(), a_.GetHandle(), Gamma_.GetHandle());
			kineticsSurfaceMap_.FormationRates(RfromSurface_.GetHandle(), Rsurface_.GetHandle(), Rbulk_.GetHandle(), RsurfacePhases_.GetHandle());
			QRSurface_ = kineticsSurfaceMap_.HeatRelease(RfromSurface_.GetHandle(), Rsurface_.GetHandle(), Rbulk_.GetHandle());
		}

		double dG_over_dz = Dot(RfromSurface_.Size(), RfromSurface_.GetHandle(), thermodynamicsMap_.MWs().data()) / cross_section_over_perimeter_;

		// Recovering residuals
		{
			unsigned int k = 1;

			// Homogeneous species
			for (unsigned int i = 1; i <= NC_; ++i)
				dy[k++] = (thermodynamicsMap_.MW(i-1) * RfromGas_[i] +
					-omega_[i] * dG_over_dz +
					RfromSurface_[i] / cross_section_over_perimeter_ * thermodynamicsMap_.MW(i-1)) / rho_/v_;

			// Mass flow rate
			dy[k++] = dG_over_dz;

			// Phases
			for (unsigned int i = 1; i <= SURF_NP_; ++i)
				dy[k++] = RsurfacePhases_[i]/v_ * artificial_speedup_coefficient_;

			// Heterogeneous species
			for (unsigned int i = 1; i <= SURF_NC_; ++i)
			{
				const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i - 1] + 1;
				dy[k++] = (thermodynamicsSurfaceMap_.vector_occupancies_site_species()[i - 1] * Rsurface_[i] -
					Z_[i] * dy[NC_ + 1 + index_phase]) / Gamma_[index_phase];
			}

			// Looking for the most abundant surface species to impose
			// the algebraic constaint
			if (SURF_NP_ == 1)
			{
				int index_local = 0;
				Z_.Max(&index_local);
				const int index_global = (NC_ + 1 + SURF_NP_) + index_local;
				dy[index_global] = 1. - Z_.SumElements();
			}
			else
			{
				Eigen::VectorXi index_local(SURF_NP_);
				Eigen::VectorXd maxZ(SURF_NP_);
				Eigen::VectorXd sumZ(SURF_NP_);
				sumZ.setZero();
				maxZ.setZero();

				for (unsigned int i = 1; i <= SURF_NC_; ++i)
				{
					const unsigned int iphase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i - 1];
					sumZ(iphase) += Z_[i];
					if (Z_[i] > maxZ(iphase))
					{
						maxZ(iphase) = Z_[i];
						index_local(iphase) = i;
					}
				}

				for (unsigned int i = 1; i <= SURF_NP_; ++i)
				{
					const unsigned int index_global = (NC_ + 1 + SURF_NP_) + index_local(i - 1);
					dy[index_global] = sumZ(i - 1) - 1.;
				}
			}

			// Bulk species
			for (unsigned int i = 1; i <= BULK_NC_; ++i)
			{
				const unsigned int index_bulk = NC_ + SURF_NC_ + i;
				dy[k++] = Rbulk_[i] * thermodynamicsSurfaceMap_.MW(index_bulk-1) / thermodynamicsSurfaceMap_.vector_densities_bulk_species()[i - 1];
			}

			// Energy equation
			{
				// Exchange of heat with external environment
				double Qexchange = 0.;
				if (global_thermal_exchange_coefficient_ != 0.)
					Qexchange = global_thermal_exchange_coefficient_ / cross_section_over_perimeter_ * (T_environment_ - T_);		// [W/m3]

				const double v2 = v_*v_;
				const double Beta = 1.;
				const double delta = RfromGas_.SumElements()*MW_ / rho_ / v_;	// [1/m]
				const double T_coefficient = CpMixMass_ + v2 / Beta / T_;       // [m2/s2/K]
				const double T_extra_terms = -v2 / Beta*delta;					// [m/s2]
				
				const double Qsurface = QRSurface_ / cross_section_over_perimeter_;		// [W/m3]

				dy[k++] = ((QRGas_ + Qexchange + Qsurface) / rho_ / v_ + T_extra_terms) / T_coefficient;
			}

			// Time-axial coordinate
			if (time_independent_ == true)	dy[k++] = 1.;
			else                            dy[k++] = 1./v_;

			if (time_independent_ == true)
				dy *= v_;
		}

		return 0;
	}

	int SurfacePlugFlowReactor_NonIsothermal::OdePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		ode_iteration_++;

		// Video output
		if (plugflow_options_.verbose_video() == true)
		{
			if (ode_iteration_%plugflow_options_.n_step_video() == 1 || plugflow_options_.n_step_video() == 1)
			{
				counter_file_video_++;
				if (counter_file_video_ % 100 == 1)
				{
					std::cout << std::endl;
					std::cout << std::setw(9) << std::left << "#Step";
					std::cout << std::setw(18) << std::left << "Time[s]";
					std::cout << std::endl;
				}
				std::cout << std::fixed << std::setw(9) << std::left << ode_iteration_;
				std::cout << std::scientific << std::setw(18) << std::setprecision(6) << std::left << t;
				std::cout << std::endl;
			}
		}

		return 0;
	}

	int SurfacePlugFlowReactor_NonIsothermal::DaePrint(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		dae_iteration_++;

		if (plugflow_options_.verbose_video() == true)
		{
			// Video output
			if (dae_iteration_%plugflow_options_.n_step_video() == 1 || plugflow_options_.n_step_video() == 1)
			{
				counter_file_video_++;
				if (counter_file_video_ % 100 == 1)
				{
					std::cout << std::endl;
					std::cout << std::setw(8) << std::left << "#Step";
					if (time_independent_ == true)
						std::cout << std::setw(16) << std::left << "Time[s]";
					else
						std::cout << std::setw(16) << std::left << "Space[m]";
					std::cout << std::setw(10) << std::left << "T[K]";
					std::cout << std::setw(10) << std::left << "P[atm]";
					if (BULK_NC_ != 0)
						std::cout << std::setw(16) << std::left << "bulk[micron]";
					std::cout << std::endl;
				}
				std::cout << std::fixed << std::setw(8) << std::left << dae_iteration_;
				std::cout << std::scientific << std::setw(16) << std::setprecision(6) << std::left << t;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << T_;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << P_ / 101325.;
				if (BULK_NC_ != 0)
					std::cout << std::scientific << std::setw(16) << std::setprecision(6) << std::left << thicknessBulk_.SumElements()*1e6;
				std::cout << std::endl;
			}
		}

		if (plugflow_options_.verbose_output() == true)
		{
			// ASCII file output
			if (plugflow_options_.verbose_ascii_file() == true)
			{
				if (dae_iteration_%plugflow_options_.n_step_file() == 1 || plugflow_options_.n_step_file() == 1)
				{
					counter_file_ASCII_++;
					PrintCurrentStatus(fASCII_, fBulkASCII_);
				}
			}

			// XML file output
			if (plugflow_options_.verbose_xml_file() == true)
			{
				if (dae_iteration_%plugflow_options_.n_step_file() == 1 || plugflow_options_.n_step_file() == 1)
				{
					counter_file_XML_++;
					if (time_independent_ == true)
					{
						fXML_ << tau_ << " ";
						fXML_ << csi_ << " ";
					}
					else
					{
						fXML_ << csi_ << " ";
						fXML_ << tau_ << " ";
					}
					fXML_ << T_ << " ";
					fXML_ << P_ << " ";
					fXML_ << MW_ << " ";
					fXML_ << rho_ << " ";
					fXML_ << QRGas_ << " ";
					fXML_ << v_ << " ";

					for (unsigned int i=1;i<=NC_;i++)
						fXML_ << std::setprecision(12) << omega_[i] << " ";
					fXML_ << std::endl;
				}
			}
		}

		if (plugflow_options_.sensitivity_analysis() == true)
			SensitivityAnalysis(t, y);
	
		return 0;
	}

	void SurfacePlugFlowReactor_NonIsothermal::PrintCurrentStatus(std::ostream& fOutput, std::ostream& fBulkOutput)
	{
		fOutput.setf(std::ios::scientific);
		fOutput << std::setw(20) << std::left << tau_;
		fOutput << std::setw(20) << std::left << csi_;
		fOutput << std::setw(20) << std::left << T0_;
		fOutput << std::setw(20) << std::left << P0_;
		fOutput << std::setw(20) << std::left << T_;
		fOutput << std::setw(20) << std::left << P_;
		fOutput << std::setw(20) << std::left << v_;
		fOutput << std::setw(20) << std::left << rho_;
		fOutput << std::setw(20) << std::left << MW_;
		fOutput << std::setw(20) << std::left << QRGas_;
		fOutput << std::setw(20) << std::left << QRSurface_;

		if (BULK_NC_ != 0)
		{
			fOutput << std::setw(20) << std::left << thicknessBulk_.SumElements();
			for (unsigned int i = 1; i <= BULK_NC_; i++)
				fOutput << std::setw(20) << std::left << thicknessBulk_[i];
		}

		for (unsigned int i = 1; i <= SURF_NP_; i++)
			fOutput << std::setw(30) << std::left << Gamma_[i];

		for (unsigned int i = 1; i <= SURF_NC_; i++)
			fOutput << std::setw(widths_of_output_surface_species_[i - 1]) << std::left << Z_[i];

		if (indices_of_output_species_.size() != 0)
		{
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				fOutput << std::setw(widths_of_output_species_[i]) << std::left << x_[indices_of_output_species_[i]];
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				fOutput << std::setw(widths_of_output_species_[i]) << std::left << omega_[indices_of_output_species_[i]];
		}
		else
		{
			for (unsigned int i = 1; i <= NC_; i++)
				fOutput << std::setw(widths_of_output_species_[i - 1]) << std::left << x_[i];
			for (unsigned int i = 1; i <= NC_; i++)
				fOutput << std::setw(widths_of_output_species_[i - 1]) << std::left << omega_[i];
		}
		fOutput << std::endl;

		if (BULK_NC_ != 0)
		{
			fBulkOutput.setf(std::ios::scientific);
			fBulkOutput << std::setw(20) << std::left << tau_;
			fBulkOutput << std::setw(20) << std::left << csi_;
			fBulkOutput << std::setw(20) << std::left << T0_;
			fBulkOutput << std::setw(20) << std::left << P0_;
			fBulkOutput << std::setw(20) << std::left << T_;
			fBulkOutput << std::setw(20) << std::left << P_;
			fBulkOutput << std::setw(20) << std::left << v_;
			fBulkOutput << std::setw(20) << std::left << rho_;
			fBulkOutput << std::setw(20) << std::left << MW_;
			fBulkOutput << std::setw(20) << std::left << QRGas_;
			fBulkOutput << std::setw(20) << std::left << QRSurface_;

			OpenSMOKE::OpenSMOKEVectorDouble omegaBulk(BULK_NC_);	// [kg/m2/s]
			for (unsigned int i = 1; i <= BULK_NC_; i++)
			{
				unsigned int index = NC_ + SURF_NC_ + i;
				omegaBulk[i] = Rbulk_[i] * thermodynamicsSurfaceMap_.MW(index-1);
			}

			fBulkOutput << std::setw(24) << std::left << thicknessBulk_.SumElements();	// [m]
			fBulkOutput << std::setw(24) << std::left << omegaBulk.SumElements();		// [kg/m2/s]
			fBulkOutput << std::setw(24) << std::left << omegaBulk.SumElements();		// [kg/s]

			for (unsigned int i = 1; i <= BULK_NC_; i++)
			{
				const double rho = thermodynamicsSurfaceMap_.vector_densities_bulk_species()[i - 1];	// [kg/m3]
				fBulkOutput << std::setw(24) << std::left << thicknessBulk_[i];							// [m]
				fBulkOutput << std::setw(24) << std::left << omegaBulk[i];								// [kg/m2/s]
				fBulkOutput << std::setw(24) << std::left << omegaBulk[i];								// [kg/s]
				fBulkOutput << std::setw(24) << std::left << omegaBulk[i] / rho;						// [m/s]
				fBulkOutput << std::setw(24) << std::left << omegaBulk[i] / rho * 1e3 * 3600;		// [mm/hr]
			}

			fBulkOutput << std::endl;
		}
	}


	void SurfacePlugFlowReactor_NonIsothermal::PrintParametricFinalStatus(std::ostream& fOutput)
	{
		fOutput.setf(std::ios::scientific);
		fOutput << std::setw(20) << std::left << tau_;
		fOutput << std::setw(20) << std::left << csi_;
		fOutput << std::setw(20) << std::left << T0_;
		fOutput << std::setw(20) << std::left << P0_;
		fOutput << std::setw(20) << std::left << T_;
		fOutput << std::setw(20) << std::left << P_;
		fOutput << std::setw(20) << std::left << v_;
		fOutput << std::setw(20) << std::left << rho_;
		fOutput << std::setw(20) << std::left << MW_;

		// Final values
		if (indices_of_output_species_.size() != 0)
		{
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				fOutput << std::setw(widths_of_output_species_[i]) << std::left << x_[indices_of_output_species_[i]];
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				fOutput << std::setw(widths_of_output_species_[i]) << std::left << omega_[indices_of_output_species_[i]];
		}
		else
		{
			for (unsigned int i = 1; i <= NC_; i++)
				fOutput << std::setw(widths_of_output_species_[i - 1]) << std::left << x_[i];
			for (unsigned int i = 1; i <= NC_; i++)
				fOutput << std::setw(widths_of_output_species_[i - 1]) << std::left << omega_[i];
		}

		// Initial values
		if (indices_of_output_species_.size() != 0)
		{
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				fOutput << std::setw(widths_of_output_species_[i]) << std::left << x0_[indices_of_output_species_[i]];
			for (unsigned int i = 0; i<indices_of_output_species_.size(); i++)
				fOutput << std::setw(widths_of_output_species_[i]) << std::left << omega0_[indices_of_output_species_[i]];
		}
		else
		{
			for (unsigned int i = 1; i <= NC_; i++)
				fOutput << std::setw(widths_of_output_species_[i - 1]) << std::left << x0_[i];
			for (unsigned int i = 1; i <= NC_; i++)
				fOutput << std::setw(widths_of_output_species_[i - 1]) << std::left << omega0_[i];
		}
		fOutput << std::endl;
	}

	void SurfacePlugFlowReactor_NonIsothermal::SetInletConditions()
	{
		Gamma0_ = Gamma_;
		Z0_ = Z_;
		thicknessBulk0_ = 0.;

		// Inlet enthalpy and internal energy for surfaces
		thermodynamicsSurfaceMap_.SetPressure(P0_);
		thermodynamicsSurfaceMap_.SetTemperature(T0_);
		H0_Surface_over_Volume_ = 0.;
		U0_Surface_over_Volume_ = 0.;
		for (unsigned int i = 1; i <= SURF_NC_; i++)
		{
			const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i - 1] + 1;
			H0_Surface_over_Volume_ += Z0_[i] * thermodynamicsSurfaceMap_.Species_H_over_RT()[i + NC_-1] * PhysicalConstants::R_J_kmol*T0_ *
				Gamma0_[index_phase] / cross_section_over_perimeter_;
			U0_Surface_over_Volume_ += Z0_[i] * (thermodynamicsSurfaceMap_.Species_H_over_RT()[i + NC_-1] - 1.)*PhysicalConstants::R_J_kmol*T0_ *
				Gamma0_[index_phase] / cross_section_over_perimeter_;
		}
	}

	void SurfacePlugFlowReactor_NonIsothermal::Solve(const double tf)
	{
		if (plugflow_options_.verbose_video() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << " Solving the non-isothermal plug flow reactor...                             " << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
		}

		// Initial conditions
		{
			unsigned int k = 1;
			for (unsigned int i = 1; i <= NC_; ++i)
				yDae0_[k++] = omega0_[i];
			yDae0_[k++] = G0_;
			for (unsigned int i = 1; i <= SURF_NP_; ++i)
				yDae0_[k++] = Gamma0_[i];
			for (unsigned int i = 1; i <= SURF_NC_; ++i)
				yDae0_[k++] = Z0_[i];
			for (unsigned int i = 1; i <= BULK_NC_; ++i)
				yDae0_[k++] = thicknessBulk0_;
			yDae0_[k++] = T0_;
			yDae0_[k++] = 0.;
		}

		// Equation types
		std::vector<OpenSMOKE::EquationType> equation_type(NE_DAE_);
		{
			unsigned int k = 0;
			for (unsigned int i = 1; i <= NC_; ++i)
				equation_type[k++] = OpenSMOKE::EquationType::EQUATION_TYPE_DIFFERENTIAL;
			equation_type[k++] = OpenSMOKE::EquationType::EQUATION_TYPE_DIFFERENTIAL;
			for (unsigned int i = 1; i <= SURF_NP_; ++i)
				equation_type[k++] = OpenSMOKE::EquationType::EQUATION_TYPE_DIFFERENTIAL;
			for (unsigned int i = 1; i <= SURF_NC_; ++i)
				equation_type[k++] = OpenSMOKE::EquationType::EQUATION_TYPE_ALGEBRAIC;
			for (unsigned int i = 1; i <= BULK_NC_; ++i)
				equation_type[k++] = OpenSMOKE::EquationType::EQUATION_TYPE_DIFFERENTIAL;
			equation_type[k++] = OpenSMOKE::EquationType::EQUATION_TYPE_DIFFERENTIAL;
			equation_type[k++] = OpenSMOKE::EquationType::EQUATION_TYPE_DIFFERENTIAL;
		}

		// Print initial conditions
		{
			DaeEquations(0., yDae0_, dyDae0_);
			DaePrint(0., yDae0_);
		}

		if (dae_parameters_.type() == DAE_Parameters::DAE_INTEGRATOR_OPENSMOKE)
		{
			// Minimum values
			Eigen::VectorXd yMin(NE_DAE_);
			{
				unsigned int k = 0;
				for (unsigned int i = 1; i <= NC_; ++i)
					yMin(k++) = 0.;
				yMin(k++) = 0.;
				for (unsigned int i = 1; i <= SURF_NP_; ++i)
					yMin(k++) = 0.;
				for (unsigned int i = 1; i <= SURF_NC_; ++i)
					yMin(k++) = 0.;
				for (unsigned int i = 1; i <= BULK_NC_; ++i)
					yMin(k++) = -1.e16;
				yMin(k++) = 0.;
				yMin(k++) = 0.;
			}

			// Maximum values
			Eigen::VectorXd yMax(NE_DAE_);
			{
				unsigned int k = 0;
				for (unsigned int i = 1; i <= NC_; ++i)
					yMax(k++) = 1.;
				yMax(k++) = 1.e16;
				for (unsigned int i = 1; i <= SURF_NP_; ++i)
					yMax(k++) = 1.e16;
				for (unsigned int i = 1; i <= SURF_NC_; ++i)
					yMax(k++) = 1.;
				for (unsigned int i = 1; i <= BULK_NC_; ++i)
					yMax(k++) = 1.e16;
				yMax(k++) = 1.e4;
				yMax(k++) = 1.e16;
			}

			// Initial conditions
			Eigen::VectorXd y0_eigen(yDae0_.Size());
			yDae0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Create the solver 
			typedef DaeSMOKE::KernelDense<OpenSMOKE::DAESystem_OpenSMOKE_SurfacePlugFlowReactor> denseDae;
			typedef DaeSMOKE::MethodGear<denseDae> methodGear;
			DaeSMOKE::MultiValueSolver<methodGear> dae_solver;
			dae_solver.SetReactor(this);

			// Set initial conditions
			dae_solver.SetInitialConditions(0., y0_eigen, equation_type);

			// Set linear algebra options
			dae_solver.SetLinearAlgebraSolver(dae_parameters_.linear_algebra());
			dae_solver.SetFullPivoting(dae_parameters_.full_pivoting());

			// Set relative and absolute tolerances
			dae_solver.SetAbsoluteTolerances(dae_parameters_.absolute_tolerance());
			dae_solver.SetRelativeTolerances(dae_parameters_.relative_tolerance());

			// Set minimum and maximum values
			dae_solver.SetMinimumValues(yMin);
			dae_solver.SetMaximumValues(yMax);

			// Set maximum number of steps
			if (dae_parameters_.maximum_number_of_steps() > 0)
				dae_solver.SetMaximumNumberOfSteps(dae_parameters_.maximum_number_of_steps());

			// Set maximum integration order
			if (dae_parameters_.maximum_order() > 0)
				dae_solver.SetMaximumOrder(dae_parameters_.maximum_order());

			// Set maximum step size allowed
			if (dae_parameters_.maximum_step() > 0)
				dae_solver.SetMaximumStepSize(dae_parameters_.maximum_step());

			// Set minimum step size allowed
			if (dae_parameters_.minimum_step() > 0)
				dae_solver.SetMinimumStepSize(dae_parameters_.minimum_step());

			// Set initial step size
			if (dae_parameters_.initial_step() > 0)
				dae_solver.SetFirstStepSize(dae_parameters_.initial_step());

			// Solve the system
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			DaeSMOKE::DaeStatus status = dae_solver.Solve(tf);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			// Check the solution
			if (status > 0)
			{
				dae_solver.Solution(yf_eigen);
				yDaef_.CopyFrom(yf_eigen.data());
				dae_parameters_.TransferDataFromDaeSolver(dae_solver, tEnd - tStart);
			}
		}
		#if OPENSMOKE_USE_BZZMATH == 1 
		else if (dae_parameters_.type() == DAE_Parameters::DAE_INTEGRATOR_BZZDAE)
		{
			// Minimum values
			BzzVector yMin(NE_DAE_);
			{
				unsigned int k = 1;
				for (unsigned int i = 1; i <= NC_; ++i)
					yMin[k++] = 0.;
				yMin[k++] = 0.;
				for (unsigned int i = 1; i <= SURF_NP_; ++i)
					yMin[k++] = 0.;
				for (unsigned int i = 1; i <= SURF_NC_; ++i)
					yMin[k++] = 0.;
				for (unsigned int i = 1; i <= BULK_NC_; ++i)
					yMin[k++] = -1.e16;
				yMin[k++] = 0.;
				yMin[k++] = 0.;
			}

			// Maximum values
			BzzVector yMax(NE_DAE_);
			{
				unsigned int k = 1;
				for (unsigned int i = 1; i <= NC_; ++i)
					yMax[k++] = 1.;
				yMax[k++] = 1.e16;
				for (unsigned int i = 1; i <= SURF_NP_; ++i)
					yMax[k++] = 1.e16;
				for (unsigned int i = 1; i <= SURF_NC_; ++i)
					yMax[k++] = 1.;
				for (unsigned int i = 1; i <= BULK_NC_; ++i)
					yMax[k++] = 1.e16;
				yMax[k++] = 1.e4;
				yMax[k++] = 1.e16;
			}
			
			// Initial conditions
			BzzVector yDae0_bzz(yDae0_.Size());
			yDae0_.CopyTo(yDae0_bzz.GetHandle());
			
			// Final solution
			BzzVector yDaef_bzz(yDae0_bzz.Size());

			// Differential vs Algebraic
			BzzVectorInt iDerAlg(equation_type.size());
			for (unsigned int i = 1; i <= equation_type.size(); i++)
				iDerAlg[i] = (equation_type[i - 1] == OpenSMOKE::EquationType::EQUATION_TYPE_DIFFERENTIAL) ? 1 : 0;

			DAESystem_BzzDae_SurfacePlugFlowReactor daeplugflow(*this);
			BzzDaeObject o(yDae0_bzz, 0., iDerAlg, &daeplugflow);
			
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			o.SetTolAbs(dae_parameters_.absolute_tolerance());
			o.SetTolRel(dae_parameters_.relative_tolerance());

			ptDae_SurfacePlugFlowReactor_NonIsothermal = &daeplugflow;
			o.StepPrint(DAE_Print_SurfacePlugFlowReactor_NonIsothermal);

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			yDaef_bzz = o(tf,tf);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			yDaef_.CopyFrom(yDaef_bzz.GetHandle());

			dae_parameters_.SetCPUTime(tEnd-tStart);
			dae_parameters_.SetNumberOfFunctionCalls(o.GetNumFunction());
			dae_parameters_.SetNumberOfJacobians(o.GetNumNumericalJacobian());
			dae_parameters_.SetNumberOfFactorizations(o.GetNumFactorization());
			dae_parameters_.SetNumberOfSteps(o.GetNumStep());
			dae_parameters_.SetLastOrderUsed(o.GetOrderUsed());
			dae_parameters_.SetLastStepUsed(o.GetHUsed());
		}
		#endif
		else 
		{
			// Differential vs Algebraic
			bool* algebraic_equations = new bool[equation_type.size()];
			for (unsigned int i = 0; i < equation_type.size(); i++)
				algebraic_equations[i] = (equation_type[i] == OpenSMOKE::EquationType::EQUATION_TYPE_DIFFERENTIAL) ? false : true;

			// Correct the initial derivatives for algebraic equations
			for (unsigned int i = 0; i < equation_type.size(); i++)
				if (algebraic_equations[i] == true) dyDae0_[i + 1] = 0.;

			DaeSolveOpenSourceSolvers(tf, algebraic_equations, dae_parameters_);
		}

		if (plugflow_options_.verbose_video() == true)
			dae_parameters_.Status(std::cout);

		DaeFinalStatus(tf, thermodynamicsMap_, thermodynamicsSurfaceMap_, plugflow_options_);
		DaeFinalSummary(plugflow_options_.output_path() / "FinalSummary.out", tf, thermodynamicsMap_, thermodynamicsSurfaceMap_, plugflow_options_);

		CloseAllFiles(plugflow_options_);
	}

	void SurfacePlugFlowReactor_NonIsothermal::SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		if (dae_iteration_ == 1)
		{
			// Writes the coefficients on file (only on request)
			if (dae_iteration_%plugflow_options_.n_step_file() == 1 || plugflow_options_.n_step_file() == 1)
			{
				counter_sensitivity_XML_++;
				for (unsigned int k=0;k<indices_of_sensitivity_species_.size();k++)
				{
					for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
						fSensitivityChildXML_[k] << 0. << " ";		
					fSensitivityChildXML_[k] << std::endl;	
				}
			}

			for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
				fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << 0. << " ";
			fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << std::endl;
		}
		else
		{
			// Scaling factors
			for(unsigned int j=1;j<=NC_;j++)
				scaling_Jp[j] = thermodynamicsMap_.MW(j-1)/rho_/v_;
			
			const double Beta = 1.;
			scaling_Jp[NC_+1] = 1./( rho_*v_* (CpMixMass_+v_*v_/Beta/T_));

			scaling_Jp *= v_;

			// Calculates the current Jacobian
			NumericalJacobian(t, y, J);

			// Recover variables
			{
				// Recover mass fractions
				#ifdef CHECK_MASSFRACTIONS
				for(unsigned int i=1;i<=NC_;++i)
					omega_[i] = std::max(y[i], 0.);
				#else
				for(unsigned int i=1;i<=NC_;++i)
					omega_[i] = y[i];
				#endif

				T_ = y[NC_+1];

				if (time_independent_ == false)
				{
					csi_ = t;
					tau_ = y[NC_+2];
				}
				else if (time_independent_ == true)
				{
					tau_ = t;
					csi_ = y[NC_+2];
				}

				// Calculates the pressure and the concentrations of species
				thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
				if (constant_pressure_ == true)
				{
					cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
					Product(cTot_, x_, &c_);
					rho_ = cTot_*MW_;
					v_ = G_/rho_;
				}
				else
				{
					P_ = G_/v_ * PhysicalConstants::R_J_kmol * T_ / MW_;
					cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
					Product(cTot_, x_, &c_);
					rho_ = cTot_*MW_;
				}
			}

			// Calculates the current sensitivity coefficients
			sensitivityMap_->Calculate(t, T_, P_, c_, J, scaling_Jp);

			// Writes the coefficients on file (only on request)
			if (dae_iteration_%plugflow_options_.n_step_file() == 1 || plugflow_options_.n_step_file() == 1)
			{
				counter_sensitivity_XML_++;
				for (unsigned int k=0;k<indices_of_sensitivity_species_.size();k++)
				{
					const unsigned int i = indices_of_sensitivity_species_[k];
					for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
					{
						double sum = 0.;
						for (unsigned int kk=1;kk<=NC_;kk++)
							sum += sensitivityMap_->sensitivity_coefficients()(kk-1,j-1)/thermodynamicsMap_.MW(kk-1);
						sum *= x_[i]*MW_;

						double coefficient = sensitivityMap_->sensitivity_coefficients()(i-1,j-1)*MW_/thermodynamicsMap_.MW(i-1) - sum;
						fSensitivityChildXML_[k] << coefficient << " ";
					}		
					fSensitivityChildXML_[k] << std::endl;	
				}

				for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
					fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << sensitivityMap_->sensitivity_coefficients()(NC_,j-1) << " ";
				fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << std::endl;
			}
		}
	}
}
