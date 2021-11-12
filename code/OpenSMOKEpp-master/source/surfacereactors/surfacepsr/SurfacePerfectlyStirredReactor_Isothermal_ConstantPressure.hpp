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

#include "SurfacePerfectlyStirredReactor_OdeInterfaces.h"
#include "math/PhysicalConstants.h"
#include "math/native-ode-solvers/MultiValueSolver"

namespace OpenSMOKE
{
	#if OPENSMOKE_USE_BZZMATH == 1
	ODESystem_BzzOde_SurfacePerfectlyStirredReactor* ptOde_SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure;
	void ODE_Print_SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure(BzzVector &Y, double t)
	{
		ptOde_SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure->MyPrint(Y,t);
	}
	#endif
	
	SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure::SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure(	
								OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap,
								OpenSMOKE::KineticsMap_Surface_CHEMKIN& kineticsSurfaceMap,
								OpenSMOKE::ODE_Parameters& ode_parameters,
								OpenSMOKE::SurfacePerfectlyStirredReactor_Options& psr_options,
								OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa, 
								const double T0, const double P0, const OpenSMOKE::OpenSMOKEVectorDouble& omega0,
								const double TInlet, const double PInlet, const OpenSMOKE::OpenSMOKEVectorDouble& omegaI,
								const double residence_time, const double volume, const double mass_flow_rate,
								const OpenSMOKE::OpenSMOKEVectorDouble& Z0, const OpenSMOKE::OpenSMOKEVectorDouble& Gamma0, 
								const std::vector<bool>& site_non_conservation, const double internal_area) :

		SurfacePerfectlyStirredReactor(thermodynamicsMap, kineticsMap, thermodynamicsSurfaceMap, kineticsSurfaceMap, ode_parameters, psr_options, on_the_fly_ropa)
	{
		type_ = SURFACEPERFECTLYSTIRRED_REACTOR_ISOTHERMAL_CONSTANTP;

		T0_ = T0;
		P0_ = P0;
		omega0_ = omega0;
		Z0_ = Z0;
		Gamma0_ = Gamma0;
		ChangeDimensions(thermodynamicsSurfaceMap_.number_of_bulk_species(), &massBulk0_, true);
		massBulk0_ = 0.;

		TInlet_ = TInlet;
		PInlet_ = PInlet;
		omegaInlet_ = omegaI;

		tau0_ = residence_time;
		V0_ = volume;
		mass_flow_rate_in_0_ = mass_flow_rate;
		mass_flow_rate_in_ = mass_flow_rate;
		mass_flow_rate_out_ = mass_flow_rate;
		mass_flow_rate_loss_surface_ = 0.;

		global_thermal_exchange_coefficient_ = 0.;
		exchange_area_ = 0.;
		T_environment_ = 0.;

		internal_area_ = internal_area;

		site_non_conservation_ = site_non_conservation;

		iteration_ = 0;
		counter_file_video_ = 0;
		counter_file_ASCII_ = 0;
		counter_file_XML_ = 0;
		counter_sensitivity_XML_ = 0;

		// Gas phase
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();

		// Surface phases
		SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
		SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
		SURF_NR_ = kineticsSurfaceMap_.NumberOfReactions();

		// Bulk phases
		BULK_NP_ = thermodynamicsSurfaceMap_.number_of_bulk_phases(0);
		BULK_NC_ = thermodynamicsSurfaceMap_.number_of_bulk_species();

		// Total number of equations
		NE_ = NC_ + 1 + SURF_NP_ + SURF_NC_ + BULK_NC_;

		MemoryAllocation();

		// Constraints Analysis
		if (V0_<0.)
			constraint_type_ = MASSFLOWRATE_RESIDENCETIME;
		else if (tau0_ < 0.)
			constraint_type_ = MASSFLOWRATE_VOLUME;
		else if (mass_flow_rate_in_0_ < 0.)
			constraint_type_ = VOLUME_RESIDENCETIME;

		// Inlet stream: Mole fractions and molecular weight
		MWInlet_ = thermodynamicsMap_.MolecularWeight_From_MassFractions(omegaInlet_.GetHandle());
		thermodynamicsMap_.MoleFractions_From_MassFractions(xInlet_.GetHandle(), MWInlet_, omegaInlet_.GetHandle());
		rhoInlet_ = PInlet_ * MWInlet_ / (PhysicalConstants::R_J_kmol*TInlet_);

		// Initial composition: Mole fractions and molecular weight
		MW0_ = thermodynamicsMap_.MolecularWeight_From_MassFractions(omega0_.GetHandle());
		thermodynamicsMap_.MoleFractions_From_MassFractions(x0_.GetHandle(), MW0_, omega0_.GetHandle());
		rho0_ = P0_ * MW0_ / (PhysicalConstants::R_J_kmol*T0_);

		// Constraints analysis
		if (constraint_type_ == MASSFLOWRATE_RESIDENCETIME)
			V0_ = mass_flow_rate_in_0_*tau0_/rho0_;
		else if (constraint_type_ == MASSFLOWRATE_VOLUME)
			tau0_ = rho0_*V0_/mass_flow_rate_in_0_;
		else if (constraint_type_ == VOLUME_RESIDENCETIME)
			mass_flow_rate_in_0_ = rho0_*V0_/tau0_;

		// Current values
		mass0_				= rho0_ * V0_;
		mass_flow_rate_in_	= mass_flow_rate_in_0_;
		tau_				= tau0_;
		rho_				= rho0_;
		T_					= T0_;
		P_					= P0_;
		V_					= V0_;
		mass_				= mass0_;
		omega_				= omega0_;
		MW_					= MW0_;
		QRGas_				= 0.;
		QRSurface_			= 0.;
		Z_					= Z0_;
		Gamma_				= Gamma0_;
		GammaFromEqn_		= Gamma0_;

		// Inlet properties
		thermodynamicsMap_.SetPressure(PInlet_);
		thermodynamicsMap_.SetTemperature(TInlet_);
		HInlet_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(xInlet_.GetHandle());
		UInlet_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(xInlet_.GetHandle());

		// Initial properties
		thermodynamicsMap_.SetPressure(P0_);
		thermodynamicsMap_.SetTemperature(T0_);
		H0_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x0_.GetHandle());
		U0_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(x0_.GetHandle());

		// Open all files
		OpenAllFiles();
	}

	int SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
	{
		// Recover unknowns
		{
			unsigned int k = 1;

			#ifdef CHECK_MASSFRACTIONS
			for (unsigned int i = 1; i <= NC_; ++i)
				omega_[i] = std::max(y[k++], 0.);
			mass_ = std::max(y[k++], 0.);
			for (unsigned int i = 1; i <= SURF_NP_; ++i)
				GammaFromEqn_[i] = std::max(y[k++], 0.);
			for (unsigned int i = 1; i <= SURF_NC_; ++i)
				Z_[i] = std::max(y[k++], 0.);
			for (unsigned int i = 1; i <= BULK_NC_; ++i)
				massBulk_[i] = y[k++];
			#else
			for (unsigned int i = 1; i <= NC_; ++i)
				omega_[i] = y[k++];
			mass_ = y[k++];
			for (unsigned int i = 1; i <= SURF_NP_; ++i)
				GammaFromEqn_[i] = y[k++];
			for (unsigned int i = 1; i <= SURF_NC_; ++i)
				Z_[i] = y[k++];
			for (unsigned int i = 1; i <= BULK_NC_; ++i)
				massBulk_[i] = y[k++];
			#endif

			for (unsigned int i = 1; i <= SURF_NP_; ++i)
				if (site_non_conservation_[i-1] == true)
					Gamma_[i] = GammaFromEqn_[i];
		}

		// Calculates the volume and the concentrations of species
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
		cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
		Product(cTot_, x_, &c_);
		rho_ = cTot_*MW_;
		
		// Constraints analysis
		if (constraint_type_ == MASSFLOWRATE_RESIDENCETIME)
			V_ = mass_flow_rate_in_*tau_/rho_;
		else if (constraint_type_ == MASSFLOWRATE_VOLUME)
			tau_ = rho_*V_/mass_flow_rate_in_;
		else if (constraint_type_ == VOLUME_RESIDENCETIME)
			mass_flow_rate_in_ = rho_*V_/tau_;

		// Update the mass
		mass_ = rho_*V_;

		// Calculates thermodynamic properties
		thermodynamicsMap_.SetTemperature(T_);
		thermodynamicsMap_.SetPressure(P_);
		
		const double CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
		CvMixMass_ = (CpMixMolar - PhysicalConstants::R_J_kmol) / MW_;
		CpMixMass_ = CpMixMolar / MW_;

		// Calculates homogeneous kinetics
		{
			kineticsMap_.SetTemperature(T_);
			kineticsMap_.SetPressure(P_);

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

		// Mass flow rates
		mass_flow_rate_loss_surface_ = -internal_area_*Dot(RfromSurface_.Size(), RfromSurface_.GetHandle(), thermodynamicsMap_.MWs().data());
		mass_flow_rate_out_ = mass_flow_rate_in_ - mass_flow_rate_loss_surface_;
		const double dm_over_dt = 0.;

		// Recovering residuals
		{
			unsigned int k = 1;

			// Gas species
			for (unsigned int i = 1; i <= NC_; ++i)
			dy[k++] = (omegaInlet_[i]*mass_flow_rate_in_ - omega_[i]*mass_flow_rate_out_) / mass_ 
							-omega_[i]/mass_*dm_over_dt
							+ thermodynamicsMap_.MW(i-1) * RfromGas_[i] / rho_
							+ thermodynamicsMap_.MW(i-1) * RfromSurface_[i] * internal_area_ / mass_;

			// Gas total mass
			dy[k++] = dm_over_dt;

			// Surface phases
			for (unsigned int i = 1; i <= SURF_NP_; ++i)
				dy[k++] = RsurfacePhases_[i];

			// Surface species
			for (unsigned int i = 1; i <= SURF_NC_; ++i)
			{
				const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i - 1] + 1;
				dy[k++] = ( thermodynamicsSurfaceMap_.vector_occupancies_site_species()[i - 1] * Rsurface_[i] -
								Z_[i] * dy[NC_ + 1 + index_phase] ) / Gamma_[index_phase];
			}

			// Bulk species
			for (unsigned int i = 1; i <= BULK_NC_; ++i)
			{
				const unsigned int index_bulk = NC_ + SURF_NC_ + i;
				dy[k++] = Rbulk_[i] * thermodynamicsSurfaceMap_.MW(index_bulk-1) * internal_area_;
			}
		}

		return 0;
	}

	int SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		iteration_++;

		if (psr_options_.verbose_video() == true)
		{
			// Video output
			if (iteration_%psr_options_.n_step_video() == 1 || psr_options_.n_step_video() == 1)
			{
				counter_file_video_++;
				if (counter_file_video_%100 == 1)
				{
					std::cout << std::endl;
					std::cout << std::setw(7)  << std::left << "#Step";
					std::cout << std::setw(15) << std::left << "Time[s]";
					std::cout << std::setw(10) << std::left << "T[K]";
					std::cout << std::setw(10) << std::left << "P[atm]";
					std::cout << std::setw(14) << std::left << "V[cm3]";
					std::cout << std::setw(14) << std::left << "min[g/s]";
					std::cout << std::setw(14) << std::left << "mout[g/s]";
					std::cout << std::setw(14) << std::left << "mloss[g/s]";
					std::cout << std::endl;
				}
				std::cout << std::fixed << std::setw(7) << std::left << iteration_;
				std::cout << std::scientific << std::setw(15) << std::setprecision(6) << std::left << t;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(2) << T_;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(4) << P_/101325.;
				std::cout << std::setw(14) << std::left << std::scientific << std::setprecision(3) << V_*1.e6;
				std::cout << std::setw(14) << std::left << std::scientific << std::setprecision(3) << mass_flow_rate_in_*1000.;
				std::cout << std::setw(14) << std::left << std::scientific << std::setprecision(3) << mass_flow_rate_out_*1000.;
				std::cout << std::setw(14) << std::left << std::scientific << std::setprecision(3) << mass_flow_rate_loss_surface_*1000.;
				std::cout << std::endl;
			}
		}

		// ASCII file output
		if (psr_options_.verbose_output() == true)
		{
			if (psr_options_.verbose_ascii_file() == true)
			{
				if (iteration_%psr_options_.n_step_file() == 1 || psr_options_.n_step_file() == 1)
				{
					counter_file_ASCII_++;
					fASCIITime_.setf(std::ios::scientific);
					fASCIITime_ << std::setw(20) << std::left << t;
					fASCIITime_ << std::setw(20) << std::left << T0_;
					fASCIITime_ << std::setw(20) << std::left << P0_;
					fASCIITime_ << std::setw(20) << std::left << V0_;
					fASCIITime_ << std::setw(20) << std::left << T_;
					fASCIITime_ << std::setw(20) << std::left << P_;
					fASCIITime_ << std::setw(20) << std::left << V_;
					fASCIITime_ << std::setw(20) << std::left << rho_;
					fASCIITime_ << std::setw(20) << std::left << MW_;
					fASCIITime_ << std::setw(20) << std::left << QRGas_;
					fASCIITime_ << std::setw(20) << std::left << QRSurface_;

					if (BULK_NC_ != 0)
					{
						fASCIITime_ << std::setw(20) << std::left << massBulk_.SumElements();
						for (unsigned int i = 1; i <= BULK_NC_; i++)
							fASCIITime_ << std::setw(20) << std::left << massBulk_[i];
					}

					for (unsigned int i = 1; i <= SURF_NP_; i++)
						fASCIITime_ << std::setw(20) << std::left << Gamma_[i];

					for (unsigned int i = 1; i <= SURF_NC_; i++)
						fASCIITime_ << std::setw(widths_of_output_surface_species_[i - 1]) << std::left << Z_[i];

					if (indices_of_output_species_.size() != 0)
					{
						for (unsigned int i = 0; i < indices_of_output_species_.size(); i++)
							fASCIITime_ << std::setw(widths_of_output_species_[i]) << std::left << x_[indices_of_output_species_[i]];
						for (unsigned int i = 0; i < indices_of_output_species_.size(); i++)
							fASCIITime_ << std::setw(widths_of_output_species_[i]) << std::left << omega_[indices_of_output_species_[i]];
					}
					else
					{
						for (unsigned int i = 1; i <= NC_; i++)
							fASCIITime_ << std::setw(widths_of_output_species_[i - 1]) << std::left << x_[i];
						for (unsigned int i = 1; i <= NC_; i++)
							fASCIITime_ << std::setw(widths_of_output_species_[i - 1]) << std::left << omega_[i];
					}
					fASCIITime_ << std::endl;

					if (BULK_NC_ != 0)
					{
						fBulkASCIITime_.setf(std::ios::scientific);
						fBulkASCIITime_ << std::setw(20) << std::left << t;
						fBulkASCIITime_ << std::setw(20) << std::left << T_;
						fBulkASCIITime_ << std::setw(20) << std::left << P_;
						fBulkASCIITime_ << std::setw(20) << std::left << V_;
						fBulkASCIITime_ << std::setw(20) << std::left << rho_;
						fBulkASCIITime_ << std::setw(20) << std::left << MW_;
						fBulkASCIITime_ << std::setw(20) << std::left << QRGas_;
						fBulkASCIITime_ << std::setw(20) << std::left << QRSurface_;

						OpenSMOKE::OpenSMOKEVectorDouble omegaBulk(BULK_NC_);	// [kg/m2/s]
						for (unsigned int i = 1; i <= BULK_NC_; i++)
						{
							unsigned int index = NC_ + SURF_NC_ + i;
							omegaBulk[i] = Rbulk_[i] * thermodynamicsSurfaceMap_.MW(index-1);
						}

						fBulkASCIITime_ << std::setw(24) << std::left << massBulk_.SumElements();					// [kg]
						fBulkASCIITime_ << std::setw(24) << std::left << omegaBulk.SumElements();					// [kg/m2/s]
						fBulkASCIITime_ << std::setw(24) << std::left << omegaBulk.SumElements()*internal_area_;	// [kg/s]

						for (unsigned int i = 1; i <= BULK_NC_; i++)
						{
							const double rho = thermodynamicsSurfaceMap_.vector_densities_bulk_species()[i - 1];		// [kg/m3]
							fBulkASCIITime_ << std::setw(24) << std::left << massBulk_[i];									// [kg]
							fBulkASCIITime_ << std::setw(24) << std::left << omegaBulk[i];									// [kg/m2/s]
							fBulkASCIITime_ << std::setw(24) << std::left << omegaBulk[i] * internal_area_;					// [kg/s]
							fBulkASCIITime_ << std::setw(24) << std::left << massBulk_[i] / rho / internal_area_;			// [m]
							fBulkASCIITime_ << std::setw(24) << std::left << massBulk_[i] / rho / internal_area_ * 1e6;		// [micron]
							fBulkASCIITime_ << std::setw(24) << std::left << omegaBulk[i] / rho;							// [m/s]
							fBulkASCIITime_ << std::setw(24) << std::left << omegaBulk[i] / rho * 1e3 * 3600;				// [mm/hr]
						}

						fBulkASCIITime_ << std::endl;
					}
				}
			}
		}

		return 0;
	}

	void SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure::PrintFinalStatus(std::ostream& fOutput, std::ostream& fOutputBulk)
	{
		fOutput << std::setw(20) << std::left << tau_;
		fOutput << std::setw(20) << std::left << T0_;
		fOutput << std::setw(20) << std::left << P0_;
		fOutput << std::setw(20) << std::left << V0_;
		fOutput << std::setw(20) << std::left << T_;
		fOutput << std::setw(20) << std::left << P_;
		fOutput << std::setw(20) << std::left << V_;
		fOutput << std::setw(20) << std::left << rho_;
		fOutput << std::setw(20) << std::left << MW_;
		fOutput << std::setw(20) << std::left << QRGas_;
		fOutput << std::setw(20) << std::left << QRSurface_;

		if (BULK_NC_ != 0)
		{
			fOutput << std::setw(20) << std::left << massBulk_.SumElements();
			for (unsigned int i = 1; i <= BULK_NC_; i++)
				fOutput << std::setw(20) << std::left << massBulk_[i];
		}

		for (unsigned int i = 1; i <= SURF_NP_; i++)
			fOutput << std::setw(20) << std::left << Gamma_[i];

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
			fOutputBulk.setf(std::ios::scientific);
			fOutputBulk << std::setw(20) << std::left << tau_;
			fOutputBulk << std::setw(20) << std::left << T_;
			fOutputBulk << std::setw(20) << std::left << P_;
			fOutputBulk << std::setw(20) << std::left << V_;
			fOutputBulk << std::setw(20) << std::left << rho_;
			fOutputBulk << std::setw(20) << std::left << MW_;
			fOutputBulk << std::setw(20) << std::left << QRGas_;
			fOutputBulk << std::setw(20) << std::left << QRSurface_;

			OpenSMOKE::OpenSMOKEVectorDouble omegaBulk(BULK_NC_);	// [kg/m2/s]
			for (unsigned int i = 1; i <= BULK_NC_; i++)
			{
				unsigned int index = NC_ + SURF_NC_ + i;
				omegaBulk[i] = Rbulk_[i] * thermodynamicsSurfaceMap_.MW(index-1);
			}

			fOutputBulk << std::setw(24) << std::left << massBulk_.SumElements();					// [kg]
			fOutputBulk << std::setw(24) << std::left << omegaBulk.SumElements();					// [kg/m2/s]
			fOutputBulk << std::setw(24) << std::left << omegaBulk.SumElements()*internal_area_;	// [kg/s]

			for (unsigned int i = 1; i <= BULK_NC_; i++)
			{
				const double rho = thermodynamicsSurfaceMap_.vector_densities_bulk_species()[i - 1];		// [kg/m3]
				fOutputBulk << std::setw(24) << std::left << massBulk_[i];									// [kg]
				fOutputBulk << std::setw(24) << std::left << omegaBulk[i];									// [kg/m2/s]
				fOutputBulk << std::setw(24) << std::left << omegaBulk[i] * internal_area_;					// [kg/s]
				fOutputBulk << std::setw(24) << std::left << massBulk_[i] / rho / internal_area_;			// [m]
				fOutputBulk << std::setw(24) << std::left << massBulk_[i] / rho / internal_area_ * 1e6;		// [micron]
				fOutputBulk << std::setw(24) << std::left << omegaBulk[i] / rho;							// [m/s]
				fOutputBulk << std::setw(24) << std::left << omegaBulk[i] / rho * 1e3 * 3600;				// [mm/hr]
			}

			fOutputBulk << std::endl;
		}
	}

	void SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure::PrintParametricFinalStatus(std::ostream& fOutput)
	{
		fOutput << std::setw(20) << std::left << tau_;
		fOutput << std::setw(20) << std::left << T0_;
		fOutput << std::setw(20) << std::left << P0_;
		fOutput << std::setw(20) << std::left << V0_;
		fOutput << std::setw(20) << std::left << T_;
		fOutput << std::setw(20) << std::left << P_;
		fOutput << std::setw(20) << std::left << V_;
		fOutput << std::setw(20) << std::left << rho_;
		fOutput << std::setw(20) << std::left << MW_;

		// Outlet values
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

		// Inlet values
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

	void SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure::PrintSolution()
	{
		if (psr_options_.verbose_output() == true)
		{
			// ASCII file output
			if (psr_options_.verbose_ascii_file() == true)
			{
				counter_file_ASCII_++;
				fASCIIFinal_.setf(std::ios::scientific);
				PrintFinalStatus(fASCIIFinal_, fBulkASCIIFinal_);
			}

			// XML file output
			for (unsigned int k = 0; k < 3; k++)
			{
				if (psr_options_.verbose_xml_file() == true)
				{
					counter_file_XML_++;
					fXML_ << double(k) / 2. << " ";
					fXML_ << T_ << " ";
					fXML_ << P_ << " ";
					fXML_ << MW_ << " ";
					fXML_ << rho_ << " ";
					fXML_ << QRGas_ << " ";
					for (unsigned int i = 1; i <= NC_; i++)
						fXML_ << std::setprecision(12) << omega_[i] << " ";
					fXML_ << std::endl;
				}
			}

			// Rate of Production Analysis (on the fly)
			if (on_the_fly_ropa_.is_active() == true)
				on_the_fly_ropa_.Analyze(fROPA_, 0, tau_, T_, P_, c_, omega_, omegaInlet_);
		}
	}

	void SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure::Solve(const double tf)
	{
		if (psr_options_.verbose_video() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << " Solving the steady-state, isothermal perfectly stirred reactor...           " << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
		}

		// Initial conditions
		{
			unsigned int k = 1;
			for (unsigned int i = 1; i <= NC_; ++i)
				y0_[k++] = omega0_[i];
			y0_[k++] = mass0_;
			for (unsigned int i = 1; i <= SURF_NP_; ++i)
				y0_[k++] = Gamma0_[i];
			for (unsigned int i = 1; i <= SURF_NC_; ++i)
				y0_[k++] = Z0_[i];
			for (unsigned int i = 1; i <= BULK_NC_; ++i)
				y0_[k++] = massBulk0_[i];
		}

		// Print intial conditions
		{
			OpenSMOKE::OpenSMOKEVectorDouble dy0(y0_.Size());
			Equations(0., y0_, dy0);
			Print(0., y0_);
		}

		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			// Minimum values
			Eigen::VectorXd yMin(NE_);
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
			}

			// Maximum values
			Eigen::VectorXd yMax(NE_);
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
			}

			// Initial conditions
			Eigen::VectorXd y0_eigen(y0_.Size());
			y0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Create the solver
			typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE_SurfacePerfectlyStirredReactor> denseOde;
			typedef OdeSMOKE::MethodGear<denseOde> methodGear;
			OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
			ode_solver.SetReactor(this);

			// Set initial conditions
			ode_solver.SetInitialConditions(0., y0_eigen);

			// Set linear algebra options
			ode_solver.SetLinearAlgebraSolver(ode_parameters_.linear_algebra());
			ode_solver.SetFullPivoting(ode_parameters_.full_pivoting());

			// Set relative and absolute tolerances
			ode_solver.SetAbsoluteTolerances(ode_parameters_.absolute_tolerance());
			ode_solver.SetRelativeTolerances(ode_parameters_.relative_tolerance());

			// Set minimum and maximum values
			ode_solver.SetMinimumValues(yMin);
			ode_solver.SetMaximumValues(yMax);

			// Set maximum number of steps
			if (ode_parameters_.maximum_number_of_steps() > 0)
				ode_solver.SetMaximumNumberOfSteps(ode_parameters_.maximum_number_of_steps());

			// Set maximum integration order
			if (ode_parameters_.maximum_order() > 0)
				ode_solver.SetMaximumOrder(ode_parameters_.maximum_order());

			// Set maximum step size allowed
			if (ode_parameters_.maximum_step() > 0)
				ode_solver.SetMaximumStepSize(ode_parameters_.maximum_step());

			// Set minimum step size allowed
			if (ode_parameters_.minimum_step() > 0)
				ode_solver.SetMinimumStepSize(ode_parameters_.minimum_step());

			// Set initial step size
			if (ode_parameters_.initial_step() > 0)
				ode_solver.SetFirstStepSize(ode_parameters_.initial_step());

			// Solve the system
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			OdeSMOKE::OdeStatus status = ode_solver.Solve(tf);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			// Check the solution
			if (status > 0)
			{
				ode_solver.Solution(yf_eigen);
				yf_.CopyFrom(yf_eigen.data());
				ode_parameters_.TransferDataFromOdeSolver(ode_solver, tEnd - tStart);
			}
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_BZZODE)
		{
			// Mininum values
			BzzVector yMin(NE_);
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
			}

			// Maximum values
			BzzVector yMax(NE_);
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
			}
			
			// Initial conditions
			BzzVector y0_bzz(y0_.Size());
			y0_.CopyTo(y0_bzz.GetHandle());
			
			// Final solution
			BzzVector yf_bzz(y0_bzz.Size());

			ODESystem_BzzOde_SurfacePerfectlyStirredReactor odepsr(*this);
			BzzOdeStiffObject o(y0_bzz, 0., &odepsr);
			
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			o.SetTolAbs(ode_parameters_.absolute_tolerance());
			o.SetTolRel(ode_parameters_.relative_tolerance());

			ptOde_SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure = &odepsr;
			o.StepPrint(ODE_Print_SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure);

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			yf_bzz = o(tf,tf);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			yf_.CopyFrom(yf_bzz.GetHandle());

			ode_parameters_.SetCPUTime(tEnd-tStart);
			ode_parameters_.SetNumberOfFunctionCalls(o.GetNumFunction());
			ode_parameters_.SetNumberOfJacobians(o.GetNumNumericalJacobian());
			ode_parameters_.SetNumberOfFactorizations(o.GetNumFactorization());
			ode_parameters_.SetNumberOfSteps(o.GetNumStep());
			ode_parameters_.SetLastOrderUsed(o.GetOrderUsed());
			ode_parameters_.SetLastStepUsed(o.GetHUsed());
		}
		#endif
		else 
		{
			SolveOpenSourceSolvers(tf);
		}

		if (psr_options_.sensitivity_analysis() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << " Performing the sensitivity analysis...                                      " << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			OpenSMOKE::OpenSMOKEVectorDouble dummy(NE_);
			Equations(tf, yf_, dummy);
			SensitivityAnalysis(tf,yf_);
		}

		if (psr_options_.verbose_video() == true)
			ode_parameters_.Status(std::cout);
	
		FinalStatus(tf);
		FinalSummary(psr_options_.output_path() / "FinalSummary.out", tf);
		PrintSolution();
		
		CloseAllFiles();
	}

	void SurfacePerfectlyStirredReactor_Isothermal_ConstantPressure::SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		// Scaling factors
		for(unsigned int j=1;j<=NC_;j++)
			scaling_Jp[j] = thermodynamicsMap_.MW(j-1)/rho_;

		// Calculates the current Jacobian
		NumericalJacobian(t, y, J);

		// Recover variables
		{
			#ifdef CHECK_MASSFRACTIONS
			for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = std::max(y[i], 0.);
			#else
				omega_ = y;
			#endif

			// Calculates the pressure and the concentrations of species
			thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
			cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
			Product(cTot_, x_, &c_);
			rho_ = cTot_*MW_;

			// Constraints analysis
			if (constraint_type_ == MASSFLOWRATE_RESIDENCETIME)
				V_ = mass_flow_rate_in_*tau_/rho_;
			else if (constraint_type_ == MASSFLOWRATE_VOLUME)
				tau_ = rho_*V_/mass_flow_rate_in_;
			else if (constraint_type_ == VOLUME_RESIDENCETIME)
				mass_flow_rate_in_ = rho_*V_/tau_;
		} 

		// Calculates the current sensitivity coefficients
		sensitivityMap_->Calculate(T_, P_, c_, J, scaling_Jp);

		// Writes the coefficients on file (only on request)
		for(unsigned int it=0;it<3;it++)	
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
		}
	}
}
