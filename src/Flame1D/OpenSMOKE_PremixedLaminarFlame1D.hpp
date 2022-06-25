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

#include "math/OpenSMOKEUtilities.h"
#include "utilities/sensitivityanalysis/SensitivityAnalysisMap_BlockTridiagonal_SteadyState.h"

// OpenSMOKE++ Solvers
#include "interfaces/Interface_OpenSMOKEppDae.h"
#include "interfaces/Interface_OpenSMOKEppNls.h"
#include "interfaces/Interface_FalseTransient_OpenSMOKEpp.h"

// Interfaces (Sundials)
#if OPENSMOKE_USE_SUNDIALS == 1
#include "interfaces/sundials_header.h"
#include "interfaces/Interface_Ida.h"
#include "interfaces/Interface_KinSol.h"
#include "interfaces/Interface_FalseTransient_KinSol.h"
#endif

// Interfaces (BzzMath)
#if OPENSMOKE_USE_BZZMATH == 1
#include "BzzMath.hpp"
#include "interfaces/Interface_BzzDae.h"
#include "interfaces/Interface_BzzNls.h"
#include "interfaces/Interface_FalseTransient_BzzNls.h"
#endif

// Interfaces (DASPK)
#if OPENSMOKE_USE_DASPK == 1
#include "interfaces/Interface_Daspk.h"
#endif

namespace OpenSMOKE
{
	double SutherlandViscosity(const double T);

	OpenSMOKE_PremixedLaminarFlame1D::OpenSMOKE_PremixedLaminarFlame1D(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
																		OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
																		OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap,
																		OpenSMOKE::Grid1D& grid) :

		thermodynamicsMap_(thermodynamicsMap),
		kineticsMap_(kineticsMap),
		transportMap_(transportMap),
		grid_(grid)
	{
		sensitivity_analysis_ = false;
		type_ = SIMULATION_TYPE_Y;
		solver_type_ = SOLVER_TYPE_FLAMESPEED;
		soret_effect_ = false;
		n_steps_video_ = 10;

		count_video_ = n_steps_video_;

		output_folder_ = "Output";
		use_dae_solver_ = true;
		use_nls_solver_ = true;

		timeFixedTemperature_ = 1.;
		timeFixedComposition_ = 1.;
		timeComplete_ = 1.;

		is_hmom_soot_ = false;
		is_polimi_soot_ = false;
		is_on_the_fly_post_processing_ = false;
		is_fixed_temperature_profile_ = false;
		is_fixed_specific_mass_flow_rate_profile_ = false;

		is_v_inlet_ = false;
		v_inlet_ = 0.;
		m_inlet_ = 0.;
		
		is_fixed_outlet_temperature_ = false;
		fixed_outlet_temperature_ = 298.15;

		is_wall_heat_exchange_ = false;
		wall_heat_exchange_coefficient_ = 0.;
		wall_heat_nusselt_number_ = 0.;
		wall_heat_internal_diameter_ = 0.01;

		gas_temperature_1st_derivative_type_ = OpenSMOKE::DERIVATIVE_1ST_BACKWARD;
		gas_mass_fractions_1st_derivative_type_ = OpenSMOKE::DERIVATIVE_1ST_BACKWARD;

		// Diffusion coefficients (by default they are calculated using the molecular theory of gases)
		mass_diffusion_coefficients_type_ = MASS_DIFFUSION_COEFFICIENTS_TYPE_MOLECULAR_THEORY_GASES;

		radiative_heat_transfer_ = false;
		environment_temperature_ = 298.15;

		MemoryAllocation();
		SetAlgebraicDifferentialEquations();
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetSolverType(const std::string solver_type)
	{
		if (solver_type == "FlameSpeed")			solver_type_ = SOLVER_TYPE_FLAMESPEED;
		else if (solver_type == "BurnerStabilized")	solver_type_ = SOLVER_TYPE_BURNERSTABILIZED;
		else FatalErrorMessage("Unknown solver type. The allowed solver types are: FlameSpeed | BurnerStabilized");
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetSoret(const bool soret)
	{
		soret_effect_ = soret;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetRadiativeHeatTransfer(const bool radiative_heat_transfer)
	{
		radiative_heat_transfer_ = radiative_heat_transfer;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetEnvironmentTemperature(const double environment_temperature)
	{
		environment_temperature_ = environment_temperature;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetPolimiSoot(OpenSMOKE::PolimiSoot_Analyzer* polimi_soot_analyzer)
	{
		polimi_soot_analyzer_ = polimi_soot_analyzer;
		is_polimi_soot_ = true;

		if (polimi_soot_analyzer_->thermophoretic_effect() == true)
		{
			j_thermophoretic_star_.resize(grid_.Np() - 1);
			for (int i = 0; i < grid_.Np() - 1; i++)
				j_thermophoretic_star_[i].resize(polimi_soot_analyzer_->bin_indices().size());
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetOnTheFlyPostProcessing(OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing)
	{
		on_the_fly_post_processing_ = on_the_fly_post_processing;
		is_on_the_fly_post_processing_ = true;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetFixedTemperatureProfile(const OpenSMOKE::OpenSMOKEVectorDouble& x, const OpenSMOKE::OpenSMOKEVectorDouble& T)
	{
		if (solver_type_ == SOLVER_TYPE_FLAMESPEED)
			FatalErrorMessage("The @FixedTemperatureProfile option cannot be used for calculating laminar flames speed");

		is_fixed_temperature_profile_ = true;
		fixed_temperature_profile_ = new FixedProfile(x.Size(), x.GetHandle(), T.GetHandle());
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetFixedSpecificMassFlowRateProfile(const OpenSMOKE::OpenSMOKEVectorDouble& x, const OpenSMOKE::OpenSMOKEVectorDouble& m)
	{
		if (solver_type_ == SOLVER_TYPE_FLAMESPEED)
			FatalErrorMessage("The @FixedSpecificMassFlowRateProfile option cannot be used for calculating laminar flames speed");

		is_fixed_specific_mass_flow_rate_profile_ = true;
		fixed_specific_mass_flow_rate_profile_ = new FixedProfile(x.Size(), x.GetHandle(), m.GetHandle());
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetFixedOutletTemperature(const double fixed_outlet_temperature)
	{
		is_fixed_outlet_temperature_ = true;
		fixed_outlet_temperature_ = fixed_outlet_temperature;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetWallHeatExchange(	const double wall_heat_exchange_coefficient, const double wall_heat_nusselt_number, const double wall_heat_internal_diameter,
																const OpenSMOKE::OpenSMOKEVectorDouble& x, const OpenSMOKE::OpenSMOKEVectorDouble& T)
	{
		is_wall_heat_exchange_ = true;

		wall_heat_exchange_coefficient_ = wall_heat_exchange_coefficient;
		wall_heat_nusselt_number_ = wall_heat_nusselt_number;
		wall_heat_internal_diameter_ = wall_heat_internal_diameter;
		wall_heat_temperature_profile_ = new FixedProfile(x.Size(), x.GetHandle(), T.GetHandle());
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetLewisNumbers(const std::vector<double> lewis_numbers)
	{
		if (lewis_numbers.size() != thermodynamicsMap_.NumberOfSpecies())
			FatalErrorMessage("The provided Lewis numbers are not correct");

		mass_diffusion_coefficients_type_ = MASS_DIFFUSION_COEFFICIENTS_TYPE_LEWIS_NUMBERS;
		lewis_numbers_ = lewis_numbers;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetupForFlameSpeed(const Eigen::VectorXd& w)
	{
		// First guess composition
		for (int i = 0; i < grid_.Np(); i++)
			Y_[i] = Y_inlet_ + (Y_[grid_.Np() - 1] - Y_inlet_)*w(i);

		// First guess temperature
		for (int i = 0; i < grid_.Np(); i++)
			T_(i) = T_inlet_ + (T_outlet_ - T_inlet_)*w(i);

		// First guess mass flow rate (i.e. velocity)
		{
			// Molecular weight
			aux_Y.CopyFrom(Y_[0].data());
			thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(0), aux_Y.GetHandle());

			// Density
			const double rho = P_*mw_(0) / PhysicalConstants::R_J_kmol / T_inlet_;

			// Specific mass flow rate [kg/m2/s]
			if (is_v_inlet_ == true)	m_inlet_ = rho * v_inlet_;
			if (is_v_inlet_ == false)	v_inlet_ = m_inlet_ / rho;
			m_.setConstant(m_inlet_);

			// Fixed specific mass flow rate [kg/m2/s]
			if (is_fixed_specific_mass_flow_rate_profile_ == true)
			{
				fixed_specific_mass_flow_rate_profile_->Interpolate(grid_.x(), m_);
				m_inlet_ = m_(0);
				v_inlet_ = m_inlet_ / rho;
			}
		}

		Properties();
		DiffusionFluxes();

		fixed_T_ = T_(grid_.fixed_point());
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetupForBurnerStabilized(const Eigen::VectorXd& w)
	{
		if (is_fixed_temperature_profile_ == false)
		{
			SetupForFlameSpeed(w);
		}
		else
		{
			fixed_temperature_profile_->Interpolate(grid_.x(), T_);

			// First guess composition
			for (int i = 0; i < grid_.Np(); i++)
				Y_[i] = Y_inlet_ + (Y_[grid_.Np() - 1] - Y_inlet_)*w(i);

			// Mass flow rate (i.e. velocity)
			{
				// Molecular weight
				aux_Y.CopyFrom(Y_[0].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(0), aux_Y.GetHandle());

				// Density
				const double rho = P_*mw_(0) / PhysicalConstants::R_J_kmol / T_inlet_;

				// Specific mass flow rate [kg/m2/s]
				if (is_v_inlet_ == true)	m_inlet_ = rho * v_inlet_;
				if (is_v_inlet_ == false)	v_inlet_ = m_inlet_ / rho;
				m_.setConstant(m_inlet_);

				// Fixed specific mass flow rate [kg/m2/s]
				if (is_fixed_specific_mass_flow_rate_profile_ == true)
				{
					fixed_specific_mass_flow_rate_profile_->Interpolate(grid_.x(), m_);
					m_inlet_ = m_(0);
					v_inlet_ = m_inlet_ / rho;
				}
			}

			Properties();
			DiffusionFluxes();

			fixed_T_ = T_(grid_.fixed_point());
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::ChangeInletConditions(const double Tinlet, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omegaInlet)
	{
		SetInlet(Tinlet, P_Pa, omegaInlet);
		Y_[0] = Y_inlet_;
		T_(0) = T_inlet_;

		// Adjust mass flow rate (i.e. velocity)
		{
			// Molecular weight
			aux_Y.CopyFrom(Y_[0].data());
			thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(0), aux_Y.GetHandle());

			// Density
			const double rho = P_*mw_(0) / PhysicalConstants::R_J_kmol / T_inlet_;

			// Specific mass flow rate [kg/m2/s]
			m_inlet_ = rho/rho_(0)*m_(0);
			m_.setConstant(m_inlet_);

			// Fixed specific mass flow rate [kg/m2/s]
			if (is_fixed_specific_mass_flow_rate_profile_ == true)
			{
				fixed_specific_mass_flow_rate_profile_->Interpolate(grid_.x(), m_);
				m_inlet_ = m_(0);
				v_inlet_ = m_inlet_ / rho;
			}
		}

		Properties();
		DiffusionFluxes();

		fixed_T_ = T_(grid_.fixed_point());
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetInlet(const double Tinlet, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& omegaInlet)
	{
		Y_inlet_.resize(thermodynamicsMap_.NumberOfSpecies());
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			Y_inlet_(j) = omegaInlet[j+1];

		P_ = P_Pa;

		T_inlet_ = Tinlet;
	}
	
	void OpenSMOKE_PremixedLaminarFlame1D::SetOutlet(const double Toutlet, const OpenSMOKE::OpenSMOKEVectorDouble& omegaOutlet)
	{
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			Y_[grid_.Np()-1](j) = omegaOutlet[j + 1];

		T_outlet_ = Toutlet;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetInletVelocity(const double inlet_velocity)
	{
		is_v_inlet_ = true;
		v_inlet_ = inlet_velocity;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetInletMassFlux(const double inlet_mass_flux)
	{
		is_v_inlet_ = false;
		m_inlet_ = inlet_mass_flux;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetUseNlsSolver(const bool use_nls_solver)
	{
		use_nls_solver_ = use_nls_solver;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetUseDaeSolver(const bool use_dae_solver)
	{
		use_dae_solver_ = use_dae_solver;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetOutputFolder(const boost::filesystem::path output_folder)
	{
		output_folder_ = output_folder;
		if (!boost::filesystem::exists(output_folder_))
			OpenSMOKE::CreateDirectory(output_folder_);
	}

	int OpenSMOKE_PremixedLaminarFlame1D::BlockDimensions() const
	{
		if (type_ == SIMULATION_TYPE_YTM)
			return thermodynamicsMap_.NumberOfSpecies() + 2;
		else if (type_ == SIMULATION_TYPE_YT)
			return thermodynamicsMap_.NumberOfSpecies() + 1;
		else if (type_ == SIMULATION_TYPE_Y)
			return thermodynamicsMap_.NumberOfSpecies();
		else if (type_ == SIMULATION_TYPE_TM)
			return 2;
		else if (type_ == SIMULATION_TYPE_HMOM)
			return hmom_->n_moments();
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
			return (thermodynamicsMap_.NumberOfSpecies()+hmom_->n_moments());
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
			return (thermodynamicsMap_.NumberOfSpecies() + 1 + hmom_->n_moments());
		else
		{
			ErrorMessage("int OpenSMOKE_PremixedLaminarFlame1D::BlockDimensions() const", "Type not provided");
			return 0;
		}
	}

	int OpenSMOKE_PremixedLaminarFlame1D::NumberOfEquations() const
	{
		return BlockDimensions()*grid_.Np();
	}

	void OpenSMOKE_PremixedLaminarFlame1D::MemoryAllocation()
	{
		// Must be allocated only the first time
		if (T_.size() == 0)
		{
			OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &aux_X, true);
			OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &aux_Y, true);
			OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &aux_C, true);
			OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &aux_R, true);
			OpenSMOKE::ChangeDimensions(thermodynamicsMap_.NumberOfSpecies(), &aux_prov, true);
		}

		T_.resize(grid_.Np());
		U_.resize(grid_.Np());
		m_.resize(grid_.Np());

		rho_.resize(grid_.Np());
		rho_.resize(grid_.Np());
		mw_.resize(grid_.Np());
		cp_.resize(grid_.Np());
		lambda_.resize(grid_.Np());
		Q_.resize(grid_.Np());

		cp_species_.resize(grid_.Np());
		for (int i = 0; i <grid_.Np(); i++)
			cp_species_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		gamma_fick_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			gamma_fick_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		gamma_fick_star_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			gamma_fick_star_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		gamma_soret_star_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			gamma_soret_star_[i].resize(transportMap_.iThermalDiffusionRatios().size());

		X_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			X_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		Y_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			Y_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		omega_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			omega_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		j_star_.resize(grid_.Np() - 1);
		for (int i = 0; i < grid_.Np() - 1; i++)
			j_star_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		j_fick_star_.resize(grid_.Np() - 1);
		for (int i = 0; i < grid_.Np() - 1; i++)
			j_fick_star_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		j_soret_star_.resize(grid_.Np() - 1);
		for (int i = 0; i < grid_.Np() - 1; i++)
			j_soret_star_[i].resize(transportMap_.iThermalDiffusionRatios().size());

		if (is_polimi_soot_ == true)
		{
			if (polimi_soot_analyzer_->thermophoretic_effect() == true)
			{
				j_thermophoretic_star_.resize(grid_.Np() - 1);
				for (int i = 0; i < grid_.Np() - 1; i++)
					j_thermophoretic_star_[i].resize(polimi_soot_analyzer_->bin_indices().size());
			}
		}

		jc_star_.resize(grid_.Np() - 1);


		dX_over_dx_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			dX_over_dx_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		dY_over_dx_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			dY_over_dx_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		dT_over_dx_.resize(grid_.Np());
		dT_over_dx_centered_.resize(grid_.Np());
		lambda_d2T_over_dx2_.resize(grid_.Np());

		// Time derivatives
		dY_over_dt_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			dY_over_dt_[i].resize(thermodynamicsMap_.NumberOfSpecies());

		dT_over_dt_.resize(grid_.Np());
		dm_over_dt_.resize(grid_.Np());

		// Radiative heat transfer
		Q_radiation_.resize(grid_.Np());
		planck_mean_absorption_gas_.resize(grid_.Np());
		planck_mean_absorption_soot_.resize(grid_.Np());
		Q_radiation_.setZero();
		planck_mean_absorption_gas_.setZero();
		planck_mean_absorption_soot_.setZero();

		// Heat wall exchange
		Q_heat_wall_.resize(grid_.Np());
		Q_heat_wall_.setZero();
	}

	void OpenSMOKE_PremixedLaminarFlame1D::EnableSensitivityAnalysis(OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options)
	{
		sensitivity_analysis_ = true;
		indices_of_sensitivity_species_.resize(sensitivity_options.list_of_species().size());
		for (int i = 0; i<indices_of_sensitivity_species_.size(); i++)
			indices_of_sensitivity_species_[i] = thermodynamicsMap_.IndexOfSpecies(sensitivity_options.list_of_species()[i]) - 1;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Properties()
	{
		for (int i = 0; i < grid_.Np(); i++)
		{
			// Thermodynamics
			{
				thermodynamicsMap_.SetPressure(P_);
				thermodynamicsMap_.SetTemperature(T_(i));

				aux_Y.CopyFrom(Y_[i].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(i), aux_Y.GetHandle());
				aux_X.CopyTo(X_[i].data());

				// Concentrations [kmol/m3]
				const double cTot = P_ / PhysicalConstants::R_J_kmol / T_(i); // [kmol/m3]
				Product(cTot, aux_X, &aux_C);

				// Mixture density
				rho_(i) = cTot*mw_(i);	// [kg/m3]

				// Species constant pressure specific heats
				thermodynamicsMap_.cpMolar_Species(aux_prov.GetHandle());

				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					cp_species_[i](j) = aux_prov[j + 1] / thermodynamicsMap_.MW(j); // [J/kg/K]

				// Mixture constant pressure specific heat
				cp_(i) = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(aux_X.GetHandle());
				cp_(i) /= mw_(i); // [J/kg/K]
			}

			// Transport properties
			{
				transportMap_.SetTemperature(T_(i));
				transportMap_.SetPressure(P_);

				// Mixture thermal conductivity
				lambda_(i) = transportMap_.ThermalConductivity(aux_X.GetHandle());

				// Mixture diffusion coefficients
				if (mass_diffusion_coefficients_type_ == MASS_DIFFUSION_COEFFICIENTS_TYPE_MOLECULAR_THEORY_GASES)
				{
					transportMap_.MassDiffusionCoefficients(aux_prov.GetHandle(), aux_X.GetHandle(), transportMap_.is_species_bundling());
					aux_prov.CopyTo(gamma_fick_[i].data());

					// Correct diffusion coefficients for soot particles
					if (is_polimi_soot_ == true)
					{
						if (polimi_soot_analyzer_->physical_diffusion() == true)
						{
							const double Dref = gamma_fick_[i](polimi_soot_analyzer_->bin_diffusivity_reference_species());
							for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
							{
								const unsigned int index = polimi_soot_analyzer_->bin_indices()[jj];
								gamma_fick_[i](index) = Dref * polimi_soot_analyzer_->bin_diffusivity_correction_factors()[jj];
							}
						}

						const double reduction_coefficient = polimi_soot_analyzer_->physical_diffusion_reduction_coefficient();
						if (reduction_coefficient != 1.)
						{
							for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
							{
								const unsigned int index = polimi_soot_analyzer_->bin_indices()[jj];
								gamma_fick_[i](index) *= reduction_coefficient;
							}
						}
					}

					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						gamma_fick_star_[i](j) = gamma_fick_[i](j)*thermodynamicsMap_.MW(j) / mw_(i);
				}
				else if (mass_diffusion_coefficients_type_ == MASS_DIFFUSION_COEFFICIENTS_TYPE_LEWIS_NUMBERS)
				{
					const double alpha = lambda_(i) / rho_(i) / cp_(i);
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						gamma_fick_[i](j) = alpha / lewis_numbers_[j];

					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						gamma_fick_star_[i](j) = gamma_fick_[i](j)*thermodynamicsMap_.MW(j) / mw_(i);
				}

				// Thermal diffusion ratios
				if (soret_effect_ == true)
				{
					transportMap_.ThermalDiffusionRatios(aux_prov.GetHandle(), aux_X.GetHandle());
					for (unsigned int ii = 0; ii < transportMap_.iThermalDiffusionRatios().size(); ii++)
					{
						unsigned int index = transportMap_.iThermalDiffusionRatios()[ii];
						gamma_soret_star_[i](ii) = gamma_fick_[i](index - 1) * aux_prov[index] * thermodynamicsMap_.MW(index-1) / mw_(i) / T_(i);
					}
				}
			}

			// Kinetics
			{
				kineticsMap_.SetTemperature(T_(i));
				kineticsMap_.SetPressure(P_);

				kineticsMap_.ReactionRates(aux_C.GetHandle());
				kineticsMap_.FormationRates(aux_R.GetHandle());
				Q_(i) = kineticsMap_.HeatRelease(aux_R.GetHandle()); // [J/s/m3]
				ElementByElementProduct(aux_R.Size(), aux_R.GetHandle(), thermodynamicsMap_.MWs().data(), aux_R.GetHandle()); // [kg/m3/s]
				aux_R.CopyTo(omega_[i].data());
			}

			// Radiative heat transfer
			if (radiative_heat_transfer_ == true)
			{
				// Gaseous phase
				{
					planck_mean_absorption_gas_(i) = transportMap_.kPlanckMix(aux_X.GetHandle());
				}

				// Soot particles
				{
					if (is_hmom_soot_ == true)
					{
						if (hmom_->radiative_heat_transfer() == true)
						{
							hmom_->SetNormalizedMoments(hmom_M_[i](0), hmom_M_[i](1), hmom_M_[i](2), hmom_M_[i](3));
							planck_mean_absorption_soot_(i) = hmom_->planck_coefficient(T_(i), hmom_->SootVolumeFraction());
						}
					}
					else if (is_polimi_soot_ == true)
					{
						if (polimi_soot_analyzer_->radiative_heat_transfer() == true)
						{
							Eigen::VectorXd mass_fractions(thermodynamicsMap_.NumberOfSpecies());
							for (unsigned int j = 1; j <= thermodynamicsMap_.NumberOfSpecies(); j++)
								mass_fractions(j - 1) = aux_Y[j];

							planck_mean_absorption_soot_(i) = polimi_soot_analyzer_->planck_coefficient(rho_(i), T_(i), mass_fractions);
						}
					}
				}

				// Total Planck mean absorption coefficient [1/m]
				const double kPlanck = planck_mean_absorption_gas_(i) + planck_mean_absorption_soot_(i);
				
				// Total radiative heat transfer [W/m3]
				Q_radiation_(i) = 4.*PhysicalConstants::sigmaStefanBoltzmann*kPlanck * (boost::math::pow<4>(T_(i)) - boost::math::pow<4>(environment_temperature_));
			}

			// Heat wall exchange [W/m3]
			if (is_wall_heat_exchange_ == true)
			{
				if (wall_heat_nusselt_number_ > 0.)
					wall_heat_exchange_coefficient_ = wall_heat_nusselt_number_ * lambda_(i) / wall_heat_internal_diameter_;	// [W/m2/K]
				
				const double wall_temperature = wall_heat_temperature_profile_->Interpolate(grid_.x()(i));

				Q_heat_wall_(i) = -4.*wall_heat_exchange_coefficient_ / wall_heat_internal_diameter_ * (wall_temperature - T_(i));	// [W/m3]
			}
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::DiffusionFluxes()
	{
		// Fick diffusion velocity
		{
			grid_.Derivative(DERIVATIVE_1ST_BACKWARD, m_, X_, &dX_over_dx_);
			for (int i = 0; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					j_fick_star_[i](j) = -0.50 * (rho_(i)*gamma_fick_star_[i](j) + rho_(i+1)*gamma_fick_star_[i + 1](j)) *dX_over_dx_[i + 1](j);
			}
		}

		// Soret diffusion velocity
		if (soret_effect_ == true)
		{
			grid_.Derivative(DERIVATIVE_1ST_BACKWARD, m_, T_, &dT_over_dx_);
			for (int i = 0; i < grid_.Ni(); i++)
			{
				for (unsigned int jj = 0; jj < transportMap_.iThermalDiffusionRatios().size(); jj++)
					j_soret_star_[i](jj) = -0.50 * (rho_(i)*gamma_soret_star_[i](jj) + rho_(i+1)*gamma_soret_star_[i + 1](jj)) *dT_over_dx_[i + 1];
			}
		}

		// Thermophoretic diffusion velocity
		if (is_polimi_soot_ == true)
		{
			if (polimi_soot_analyzer_->thermophoretic_effect() == true)
			{
				// Dinamyc viscosity [kg/m/s]
				Eigen::VectorXd mu(grid_.Np());
				for (int i = 0; i < grid_.Np(); i++)
					mu(i) = SutherlandViscosity(T_(i));

				grid_.Derivative(DERIVATIVE_1ST_BACKWARD, m_, T_, &dT_over_dx_);
				for (int i = 0; i < grid_.Ni(); i++)
				{
					for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
					{
						unsigned int index = polimi_soot_analyzer_->bin_indices()[jj];
						j_thermophoretic_star_[i](jj) = -0.50 * (0.538*mu(i) / T_(i)*Y_[i](index) + 0.538*mu(i + 1) / T_(i + 1)*Y_[i + 1](index)) * dT_over_dx_[i + 1];
					}
				}
			}
		}

		// Correction diffusion velocity
		{
			jc_star_.setConstant(0.);

			// Fick contribution
			for (int i = 0; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					jc_star_(i) -= j_fick_star_[i](j);
			}

			// Soret contribution
			if (soret_effect_ == true)
			{
				for (int i = 0; i < grid_.Ni(); i++)
				{
					for (unsigned int jj = 0; jj < transportMap_.iThermalDiffusionRatios().size(); jj++)
						jc_star_(i) -= j_soret_star_[i](jj);
				}
			}

			// Thermophoretic contribution
			if (is_polimi_soot_ == true)
			{
				if (polimi_soot_analyzer_->thermophoretic_effect() == true)
				{
					for (int i = 0; i < grid_.Ni(); i++)
					{
						for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
							jc_star_(i) -= j_thermophoretic_star_[i](jj);
					}
				}
			}

			if (is_polimi_soot_ == true)
			{
				for (int i = 0; i < grid_.Ni(); i++)
				{
					double omega_soot = 0.;
					for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
					{
						unsigned int index = polimi_soot_analyzer_->bin_indices()[jj];
						omega_soot += Y_[i](index);
					}
					const double omega_gas = 1. - omega_soot;
					jc_star_(i) /= omega_gas;
				}
			}
		}

		// Total diffusion velocity (TOSIMPLIFY)
		{
			// Fick + Correction velocity
			if (is_polimi_soot_ == false)
			{
				for (int i = 0; i < grid_.Ni(); i++)
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						j_star_[i](j) = j_fick_star_[i](j) + 0.50*(Y_[i](j) + Y_[i + 1](j))*jc_star_(i);
				}
			}
			else
			{
				for (int i = 0; i < grid_.Ni(); i++)
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					{
						if (!OpenSMOKE::IsValuePresent(j, polimi_soot_analyzer_->bin_indices()))
							j_star_[i](j) = j_fick_star_[i](j) + 0.50*(Y_[i](j) + Y_[i + 1](j))*jc_star_(i);
						else if (OpenSMOKE::IsValuePresent(j, polimi_soot_analyzer_->bin_indices()))
							j_star_[i](j) = j_fick_star_[i](j);
					}
				}
			}

			// Soret contribution
			if (soret_effect_ == true)
			{
				for (int i = 0; i < grid_.Ni(); i++)
				{
					for (unsigned int jj = 0; jj < transportMap_.iThermalDiffusionRatios().size(); jj++)
					{
						unsigned int index = transportMap_.iThermalDiffusionRatios()[jj];
						j_star_[i](index - 1) += j_soret_star_[i](jj);
					}
				}
			}

			// Thermophoretic contribution
			if (is_polimi_soot_ == true)
			{
				if (polimi_soot_analyzer_->thermophoretic_effect() == true)
				{
					for (int i = 0; i < grid_.Ni(); i++)
					{
						for (unsigned int jj = 0; jj < polimi_soot_analyzer_->bin_indices().size(); jj++)
						{
							unsigned int index = polimi_soot_analyzer_->bin_indices()[jj];
							j_star_[i](index) += j_thermophoretic_star_[i](jj);
						}
					}
				}
			}
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::ResidenceTime(Eigen::VectorXd& tau)
	{
		tau.resize(grid_.Np());
		tau(0) = 0.;
		for (int i = 1; i < grid_.Np(); i++)
			tau(i) = tau(i - 1) + (grid_.x()(i) - grid_.x()(i - 1)) / (0.50*(m_(i) / rho_(i) + m_(i - 1) / rho_(i - 1)));
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Equations(const double t, const double* y, double* dy)
	{
		if (type_ == SIMULATION_TYPE_YTM)
			Equations_MassFractions_Temperature_MassFlowRate(t, y, dy);
		else if (type_ == SIMULATION_TYPE_YT)
			Equations_MassFractions_Temperature(t, y, dy);
		else if (type_ == SIMULATION_TYPE_Y)
			Equations_MassFractions(t, y, dy);
		else if (type_ == SIMULATION_TYPE_TM)
			Equations_Temperature_MassFlowRate(t, y, dy);
		else if (type_ == SIMULATION_TYPE_HMOM)
			Equations_HMOM(t, y, dy);
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
			Equations_MassFractions_HMOM(t, y, dy);
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
			Equations_MassFractions_Temperature_HMOM(t, y, dy);
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SubEquations_MassFractions()
	{
		// Inlet boundary
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			dY_over_dt_[0](j) = m_(0)*Y_[0](j) + j_star_[0](j) - m_(0)*Y_inlet_(j);

		// Internal points
		for (int i = 1; i < grid_.Ni(); i++)
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			{
				dY_over_dt_[i](j) = -m_(i)*dY_over_dx_[i](j)
					- (j_star_[i](j) - j_star_[i - 1](j)) / grid_.dxc_over_2()(i)
					+ omega_[i](j);
				dY_over_dt_[i](j) /= rho_(i);
			}

		// Outlet boundary
		for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
			dY_over_dt_[grid_.Ni()](j) = Y_[grid_.Ni()](j) - Y_[grid_.Ni() - 1](j);
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SubEquations_Temperature()
	{
		// Inlet boundary
		dT_over_dt_(0) = T_(0) - T_inlet_;

		// Internal points
		for (int i = 1; i < grid_.Ni(); i++)
		{
			double sumCp = 0.;
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				sumCp += cp_species_[i](j)* (j_star_[i - 1](j) + j_star_[i](j)) / 2.;
			sumCp *= dT_over_dx_centered_[i];

			dT_over_dt_(i) = -m_(i)*cp_(i)*dT_over_dx_(i)
				+ lambda_d2T_over_dx2_(i)
				- sumCp
				+ Q_(i);
			
			if (radiative_heat_transfer_ == true)
				dT_over_dt_(i) += -Q_radiation_(i);

			if (is_wall_heat_exchange_ == true)
				dT_over_dt_(i) += -Q_heat_wall_(i);

			dT_over_dt_(i) /= (rho_(i) * cp_(i));
		}

		// Outlet boundary
		dT_over_dt_(grid_.Ni()) = T_(grid_.Ni()) - T_(grid_.Ni() - 1);

		// In case of fixed outlet temperature
		if (is_fixed_outlet_temperature_ == true)
			dT_over_dt_(grid_.Ni()) = T_(grid_.Ni()) - fixed_outlet_temperature_;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SubEquations_MassFlowRate()
	{
		// Inlet boundary
		dm_over_dt_(0) = m_(1) - m_(0);

		// Internal points
		for (int i = 1; i < grid_.Ni(); i++)
		{
			// If turned on, inconsistent algebraic equations are found
			bool include_derivative_of_density = false;
			if (include_derivative_of_density == true)
			{
				double drho_over_dt = -rho_(i) / T_(i)*dT_over_dt_(i);
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					drho_over_dt -= rho_(i)*mw_(i) / thermodynamicsMap_.MW(j) * dY_over_dt_[i](j);

				if (i < static_cast<int>(grid_.fixed_point()))
					dm_over_dt_(i) = (m_(i + 1) - m_(i)) / (grid_.x()(i + 1) - grid_.x()(i)) + drho_over_dt;
				else
					dm_over_dt_(i) = (m_(i) - m_(i - 1)) / (grid_.x()(i) - grid_.x()(i - 1)) + drho_over_dt;
			}
			else
			{
				if (i < static_cast<int>(grid_.fixed_point()))
					dm_over_dt_(i) = m_(i) - m_(i + 1);
				else
					dm_over_dt_(i) = m_(i - 1) - m_(i);
			}

			if (i == grid_.fixed_point())
				dm_over_dt_(grid_.fixed_point()) = T_(grid_.fixed_point()) - fixed_T_;
		}

		// Outlet boundary
		dm_over_dt_(grid_.Ni()) = m_(grid_.Ni()) - m_(grid_.Ni() - 1);
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SubEquations_HMOM()
	{
		const int jN2 = thermodynamicsMap_.IndexOfSpecies("N2") - 1;
		const int jOH = thermodynamicsMap_.IndexOfSpecies("OH")-1;
		const int jH = thermodynamicsMap_.IndexOfSpecies("H") - 1;
		const int jH2O = thermodynamicsMap_.IndexOfSpecies("H2O") - 1;
		const int jH2 = thermodynamicsMap_.IndexOfSpecies("H2") - 1;
		const int jC2H2 = thermodynamicsMap_.IndexOfSpecies("C2H2") - 1;
		const int jO2 = thermodynamicsMap_.IndexOfSpecies("O2") - 1;
		std::vector<int> jPAH(hmom_->pah_species().size());
		for (unsigned int i = 0; i<hmom_->pah_species().size(); i++)
			jPAH[i] = thermodynamicsMap_.IndexOfSpecies(hmom_->pah_species()[i]) - 1;

		// Inlet boundary
		for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			dhmom_M_over_dt_[0](j) = hmom_M_[0](j) - 1e-12;

		// Internal points
		for (int i = 1; i <= grid_.Ni(); i++)
		{
			const double threshold = 1.e-16;

			if (hmom_M_[i](3) < threshold)
			{
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					dhmom_M_over_dt_[i](j) = 0.;
			}
			else
			{
				const double v = m_(i) / rho_(i);														// gas velocity [m/s]
				const double gamma = SutherlandViscosity(T_(i)) / rho_(i) / hmom_->schmidt_number();	// diffusion coefficient [m2/s]

				// Concentrations of relevant species [kmol/m3]
				const double ctot = P_ / PhysicalConstants::R_J_kmol / T_[i];
				const double conc_OH = ctot*X_[i](jOH);
				const double conc_H = ctot*X_[i](jH);
				const double conc_H2O = ctot*X_[i](jH2O);
				const double conc_H2 = ctot*X_[i](jH2);
				const double conc_C2H2 = ctot*X_[i](jC2H2);
				const double conc_O2 = ctot*X_[i](jO2);
				double conc_PAH = 0.;
				for (unsigned int j = 0; j<hmom_->pah_species().size(); j++)
					conc_PAH += ctot*X_[i](jPAH[j]);

				// Mass fractions of rlelevant species
				const double mass_fraction_OH = Y_[i](jOH);
				const double mass_fraction_H = Y_[i](jH);

				// Prepares the HMOM
				hmom_->SetNormalizedMoments(hmom_M_[i](0), hmom_M_[i](1), hmom_M_[i](2), hmom_M_[i](3));
				hmom_->SetTemperatureAndPressure(T_[i], P_);
				hmom_->SetMassFractions(mass_fraction_OH, mass_fraction_H);
				hmom_->SetConcentrations("kmol/m3", conc_OH, conc_H, conc_H2O, conc_H2, conc_C2H2, conc_O2, conc_PAH);
				hmom_->SetViscosity(SutherlandViscosity(T_[i]));

				// HMOM source terms [mol/m3/s]
				hmom_->CalculateSourceMoments();

				// Transport equations
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					dhmom_M_over_dt_[i](j) = -v*dhmom_M_over_dx_[i](j) + hmom_->sources()(j);

				// Source terms gas phase
				if (hmom_->PAHConsumption() == true)
				{
					if (conc_PAH > 1.e-64)
					{
						const double R_PAH = hmom_->PAHConsumptionRate() / 1000.;	// [kmol/m3/s]

						double Omega_PAH = 0.;	// [kg/m3]
						for (unsigned int j = 0; j < hmom_->pah_species().size(); j++)
						{
							const double fraction = ctot*X_[i](jPAH[j]) / conc_PAH;
							const double omega = thermodynamicsMap_.MW(jPAH[j]) * R_PAH * fraction;
							dY_over_dt_[i](jPAH[j]) -= omega / rho_[i];
							Omega_PAH += omega;
						}

						dY_over_dt_[i](jN2) += Omega_PAH / rho_[i];
					}
				}
			}
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Recover_Unknowns(const double* y)
	{
		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = y[count++];

				// Temperature
				T_(i) = y[count++];

				// Mass flow rate
				m_(i) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = y[count++];

				// Temperature
				T_(i) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Temperature
				T_(i) = y[count++];

				// Mass flow rate
				m_(i) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = y[count++];

				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = y[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = y[count++];
				T_(i) = y[count++];
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = y[count++];
			}
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Recover_Residuals(double* dy)
	{
		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					dy[count++] = dY_over_dt_[i](j);
				dy[count++] = dT_over_dt_(i);
				dy[count++] = dm_over_dt_(i);
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					dy[count++] = dY_over_dt_[i](j);
				dy[count++] = dT_over_dt_(i);
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				dy[count++] = dT_over_dt_(i);
				dy[count++] = dm_over_dt_(i);
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					dy[count++] = dY_over_dt_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					dy[count++] = dhmom_M_over_dt_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					dy[count++] = dY_over_dt_[i](j);

				// Moments
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					dy[count++] = dhmom_M_over_dt_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				// Species
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					dy[count++] = dY_over_dt_[i](j);

				// Temperature
				dy[count++] = dT_over_dt_(i);

				// Moments
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					dy[count++] = dhmom_M_over_dt_[i](j);
			}
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Equations_MassFractions_Temperature_MassFlowRate(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives
		{
			// Species
			grid_.Derivative(gas_mass_fractions_1st_derivative_type_, m_, Y_, &dY_over_dx_);

			// Temperature
			grid_.Derivative(gas_temperature_1st_derivative_type_, m_, T_, &dT_over_dx_);
			grid_.Derivative(DERIVATIVE_1ST_CENTERED, m_, T_, &dT_over_dx_centered_);
			grid_.SecondDerivative(lambda_, T_, &lambda_d2T_over_dx2_);
		}
		
		// Equations
		SubEquations_MassFractions();
		SubEquations_Temperature();
		SubEquations_MassFlowRate();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Equations_MassFractions_Temperature(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives
		{
			// Species
			grid_.Derivative(gas_mass_fractions_1st_derivative_type_, m_, Y_, &dY_over_dx_);

			// Temperature
			grid_.Derivative(gas_temperature_1st_derivative_type_, m_, T_, &dT_over_dx_);
			grid_.Derivative(DERIVATIVE_1ST_CENTERED, m_, T_, &dT_over_dx_centered_);
			grid_.SecondDerivative(lambda_, T_, &lambda_d2T_over_dx2_);
		}

		// Equations
		SubEquations_MassFractions();
		SubEquations_Temperature();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Equations_Temperature_MassFlowRate(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives
		{
			grid_.Derivative(gas_temperature_1st_derivative_type_, m_, T_, &dT_over_dx_);
			grid_.Derivative(DERIVATIVE_1ST_CENTERED, m_, T_, &dT_over_dx_centered_);
			grid_.SecondDerivative(lambda_, T_, &lambda_d2T_over_dx2_);
		}

		// Equations
		SubEquations_Temperature();
		SubEquations_MassFlowRate();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Equations_MassFractions(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives of Species
		grid_.Derivative(gas_mass_fractions_1st_derivative_type_, m_, Y_, &dY_over_dx_);

		// Equations
		SubEquations_MassFractions();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Equations_MassFractions_Temperature_HMOM(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives of Species
		grid_.Derivative(gas_mass_fractions_1st_derivative_type_, m_, Y_, &dY_over_dx_);
		grid_.Derivative(DERIVATIVE_1ST_BACKWARD, m_, hmom_M_, &dhmom_M_over_dx_);
		grid_.Derivative(gas_temperature_1st_derivative_type_, m_, T_, &dT_over_dx_);
		grid_.Derivative(DERIVATIVE_1ST_CENTERED, m_, T_, &dT_over_dx_centered_);
		grid_.SecondDerivative(lambda_, T_, &lambda_d2T_over_dx2_);

		// Equations
		SubEquations_MassFractions();
		SubEquations_Temperature();
		SubEquations_HMOM();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Equations_MassFractions_HMOM(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Fluxes
		DiffusionFluxes();

		// Derivatives of Species
		grid_.Derivative(gas_mass_fractions_1st_derivative_type_, m_, Y_, &dY_over_dx_);
		grid_.Derivative(DERIVATIVE_1ST_BACKWARD, m_, hmom_M_, &dhmom_M_over_dx_);

		// Equations
		SubEquations_MassFractions();
		SubEquations_HMOM();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Equations_HMOM(const double t, const double* y, double* dy)
	{
		// Recover unknowns
		Recover_Unknowns(y);

		// Properties
		Properties();

		// Derivatives of Species
		grid_.Derivative(DERIVATIVE_1ST_BACKWARD, m_, hmom_M_, &dhmom_M_over_dx_);

		// Equations
		SubEquations_HMOM();

		// Recover residuals
		Recover_Residuals(dy);
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Print(const double t, const Eigen::VectorXd& y, std::ofstream& fOutput)
	{
		// Recover unknowns
		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				OpenSMOKE::OpenSMOKEVectorDouble yy(thermodynamicsMap_.NumberOfSpecies());
				OpenSMOKE::OpenSMOKEVectorDouble xx(thermodynamicsMap_.NumberOfSpecies());

				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					yy(j + 1) = y(count++);
				const double T = y(count++);
				const double m = y(count++);

				double MW;
				thermodynamicsMap_.MoleFractions_From_MassFractions(xx.GetHandle(), MW, yy.GetHandle());

				fOutput << std::setprecision(9) << std::setw(18) << t;
				fOutput << std::setprecision(9) << std::setw(18) << grid_.x()(i);
				fOutput << std::setprecision(9) << std::setw(18) << T;
				fOutput << std::setprecision(9) << std::setw(18) << m;

				// Elements
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
					{
						double sum = 0.;
						for (unsigned int k = 0; k < thermodynamicsMap_.NumberOfSpecies(); k++)
							sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
						fOutput << std::setprecision(9) << std::setw(18) << sum*m/MW;
					}
				}

				// Species
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						fOutput << std::setprecision(9) << std::setw(18) << xx(j + 1);
				}

				fOutput << std::endl;
			}
			fOutput << std::endl;
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{

				OpenSMOKE::OpenSMOKEVectorDouble yy(thermodynamicsMap_.NumberOfSpecies());
				OpenSMOKE::OpenSMOKEVectorDouble xx(thermodynamicsMap_.NumberOfSpecies());

				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					yy(j + 1) = y(count++);
				const double T = y(count++);

				double MW;
				thermodynamicsMap_.MoleFractions_From_MassFractions(xx.GetHandle(), MW, yy.GetHandle());

				fOutput << std::setprecision(9) << std::setw(18) << t;
				fOutput << std::setprecision(9) << std::setw(18) << grid_.x()[i];
				fOutput << std::setprecision(9) << std::setw(18) << T;
				fOutput << std::setprecision(9) << std::setw(18) << m_(i);

				// Elements
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
					{
						double sum = 0.;
						for (unsigned int k = 0; k < thermodynamicsMap_.NumberOfSpecies(); k++)
							sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
						fOutput << std::setprecision(9) << std::setw(18) << sum * m_(i) / MW;
					}
				}
				
				// Species
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						fOutput << std::setprecision(9) << std::setw(18) << xx(j + 1);
				}

				fOutput << std::endl;
			}
			fOutput << std::endl;
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{

				OpenSMOKE::OpenSMOKEVectorDouble yy(thermodynamicsMap_.NumberOfSpecies());
				OpenSMOKE::OpenSMOKEVectorDouble xx(thermodynamicsMap_.NumberOfSpecies());

				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					yy(j + 1) = y(count++);

				double MW;
				thermodynamicsMap_.MoleFractions_From_MassFractions(xx.GetHandle(), MW, yy.GetHandle());

				fOutput << std::setprecision(9) << std::setw(18) << t;
				fOutput << std::setprecision(9) << std::setw(18) << grid_.x()(i);
				fOutput << std::setprecision(9) << std::setw(18) << T_(i);
				fOutput << std::setprecision(9) << std::setw(18) << m_(i);

				// Elements
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
					{
						double sum = 0.;
						for (unsigned int k = 0; k < thermodynamicsMap_.NumberOfSpecies(); k++)
							sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
						fOutput << std::setprecision(9) << std::setw(18) << sum * m_(i) / MW;
					}
				}

				// Species 
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						fOutput << std::setprecision(9) << std::setw(18) << xx(j + 1);
				}
				fOutput << std::endl;
			}
			fOutput << std::endl;
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{

				OpenSMOKE::OpenSMOKEVectorDouble yy(thermodynamicsMap_.NumberOfSpecies());
				OpenSMOKE::OpenSMOKEVectorDouble xx(thermodynamicsMap_.NumberOfSpecies());

				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					yy(j + 1) = Y_[i](j);
				const double T = y(count++);
				const double m = y(count++);

				double MW;
				thermodynamicsMap_.MoleFractions_From_MassFractions(xx.GetHandle(), MW, yy.GetHandle());

				fOutput << std::setprecision(9) << std::setw(18) << t;
				fOutput << std::setprecision(9) << std::setw(18) << grid_.x()(i);
				fOutput << std::setprecision(9) << std::setw(18) << T;
				fOutput << std::setprecision(9) << std::setw(18) << m;

				// Elements
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
					{
						double sum = 0.;
						for (unsigned int k = 0; k < thermodynamicsMap_.NumberOfSpecies(); k++)
							sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
						fOutput << std::setprecision(9) << std::setw(18) << sum*m/MW;
					}
				}

				// Species
				{
					for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
						fOutput << std::setprecision(9) << std::setw(18) << xx(j + 1);
				}
				fOutput << std::endl;
			}
			fOutput << std::endl;
		}
		else if (type_ == SIMULATION_TYPE_HMOM || type_ == SIMULATION_TYPE_Y_HMOM || type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				fOutput << std::setprecision(9) << std::setw(18) << t;
				fOutput << std::setprecision(9) << std::setw(18) << grid_.x()(i);

				// Moments
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				{
					fOutput << std::setprecision(9) << std::setw(18) << hmom_M_[i](j);;
				}

				fOutput << std::endl;
			}
			fOutput << std::endl;
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Print(const double* y, const double norm_residuals)
	{
		const double mass_flow_rate = m_(0);
		const double flame_speed = mass_flow_rate / rho_(0);
		const double max_T = T_.maxCoeff();

		std::cout << std::left << std::setw(14) << std::scientific << std::setprecision(3) << norm_residuals;
		std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(3) << flame_speed*100.;
		std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(3) << max_T;
		
		std::cout << std::endl;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Print(const double t, const double* y)
	{
		if (count_video_ == n_steps_video_)
		{
			const double mass_flow_rate = m_(0);
			const double flame_speed = mass_flow_rate / rho_(0);
			const double max_T = T_.maxCoeff();

			std::cout << std::left << std::setw(16) << std::scientific << t;
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(3) << flame_speed*100.;
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(3) << max_T;
			std::cout << std::endl;

			count_video_ = 0;
		}
		count_video_++;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Print(const std::string name_file)
	{
		std::ofstream fOutput(name_file.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		// Labels
		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "t[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "x[cm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "v[cm/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "m[kg/m2/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "rho[kg/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "MW[kg/kmol]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "phi[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "k[W/m/K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Cp[J/kg/K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Q[W/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Qrad[W/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Qwall[W/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "kPlaGas[1/m]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "kPlaSoot[1/m]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "vTherm[cm/s]", count);

			for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.elements()[j] + "[kmol/m2/s]", count);
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_x", count);
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, thermodynamicsMap_.NamesOfSpecies()[j] + "_w", count);
			fOutput << std::endl;
		}

		// Local residence time
		Eigen::VectorXd tau;
		ResidenceTime(tau);

		// Loop over all the points
		for (int i = 0; i < grid_.Np(); i++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble yy(thermodynamicsMap_.NumberOfSpecies());
			OpenSMOKE::OpenSMOKEVectorDouble xx(thermodynamicsMap_.NumberOfSpecies());

			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				yy(j + 1) = Y_[i](j);

			double MW;
			thermodynamicsMap_.MoleFractions_From_MassFractions(xx.GetHandle(), MW, yy.GetHandle());

			// Calculate thermophoretic velocity [m/s]
			const double vThermophoretic = -0.55*SutherlandViscosity(T_(i)) / rho_(i)*dT_over_dx_(i) / T_(i);

			fOutput << std::setprecision(9) << std::setw(20) << tau(i);
			fOutput << std::setprecision(9) << std::setw(20) << grid_.x()(i)*100.;
			fOutput << std::setprecision(9) << std::setw(20) << T_(i);
			fOutput << std::setprecision(9) << std::setw(20) << m_(i)/rho_(i)*100.;
			fOutput << std::setprecision(9) << std::setw(20) << m_(i);
			fOutput << std::setprecision(9) << std::setw(20) << rho_(i);
			fOutput << std::setprecision(9) << std::setw(20) << mw_(i);
			fOutput << std::setprecision(9) << std::setw(20) << thermodynamicsMap_.GetLocalEquivalenceRatioFromMoleFractions(xx.GetHandle());
			fOutput << std::setprecision(9) << std::setw(20) << lambda_(i);
			fOutput << std::setprecision(9) << std::setw(20) << cp_(i);
			fOutput << std::setprecision(9) << std::setw(20) << Q_(i);
			fOutput << std::setprecision(9) << std::setw(20) << Q_radiation_(i);
			fOutput << std::setprecision(9) << std::setw(20) << Q_heat_wall_(i);
			fOutput << std::setprecision(9) << std::setw(20) << planck_mean_absorption_gas_(i);
			fOutput << std::setprecision(9) << std::setw(20) << planck_mean_absorption_soot_(i);
			fOutput << std::setprecision(9) << std::setw(20) << vThermophoretic*100.;

			// Elements
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
				{
					double sum = 0.;
					for (unsigned int k = 0; k < thermodynamicsMap_.NumberOfSpecies(); k++)
						sum += thermodynamicsMap_.atomic_composition()(k, j) * xx(k + 1);
					fOutput << std::setprecision(9) << std::setw(20) << sum*m_(i)/MW;
				}
			}

			// Sum
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					fOutput << std::setprecision(9) << std::setw(20) << xx(j + 1);
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					fOutput << std::setprecision(9) << std::setw(20) << yy(j + 1);
			}
			fOutput << std::endl;
		}

		fOutput.close();
	}
	

	void OpenSMOKE_PremixedLaminarFlame1D::PrintSoot(const boost::filesystem::path output_folder)
	{
		// Prepare soot file
		std::ofstream fOutputSoot( (output_folder / "Solution.soot.out").c_str(), std::ios::out);
		fOutputSoot.setf(std::ios::scientific);

		// Labels
		{
			unsigned int count = 1;
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "t[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "x[cm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "v[cm/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "m[kg/m2/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "rho[kg/m3]", count);

			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "fv(L)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "x(L)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "w(L)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "rho(L)[kg/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "N(L)[#/cm3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "H/C(L)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "O/C(L)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "O/H(L)[-]", count);

			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "fv(S)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "x(S)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "w(S)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "rho(S)[kg/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "N(S)[#/cm3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "H/C(S)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "O/C(S)[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputSoot, "O/H(S)[-]", count);

			fOutputSoot << std::endl;
		}

		// Prepare soot distribution file
		std::ofstream fOutputSootDistribution((output_folder / "Solution.soot_distribution.out").c_str(), std::ios::out);
		fOutputSootDistribution.setf(std::ios::scientific);
		polimi_soot_analyzer_->WriteDistributionLabel(fOutputSootDistribution);

		// Local residence time
		Eigen::VectorXd tau;
		ResidenceTime(tau);

		// Loop over all the cells
		for (int i = 0; i < grid_.Np(); i++)
		{
			// Gas-phase properties
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << tau(i);
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << grid_.x()(i) * 100.;
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << T_(i);
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << m_(i) / rho_(i)*100.;
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << m_(i);
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << rho_(i);

			// Analysis of soot
			polimi_soot_analyzer_->Analysis(T_(i), P_, rho_(i), Y_[i], X_[i]);
			polimi_soot_analyzer_->Distribution();

			// Large sections (soot)
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->fv_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->x_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->omega_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->rho_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->N_large() / 1.e6;
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->h_over_c_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->o_over_c_large();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->o_over_h_large();

			// Small sections (PAH)
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->fv_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->x_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->omega_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->rho_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->N_small() / 1.e6;
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->h_over_c_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->o_over_c_small();
			fOutputSoot << std::scientific << std::setprecision(9) << std::setw(20) << polimi_soot_analyzer_->o_over_h_small();
			fOutputSoot << std::endl;

			// Distributions
			polimi_soot_analyzer_->WriteDistribution(fOutputSootDistribution, tau(i), grid_.x()(i), 0., 0., T_(i));
			fOutputSootDistribution << std::endl;
		}

		fOutputSoot.close();
		fOutputSootDistribution.close();
	}

	void OpenSMOKE_PremixedLaminarFlame1D::PrintHMOM(const boost::filesystem::path output_folder)
	{
		// Prepare soot file
		std::ofstream fOutputHMOM((output_folder / "Solution.hmom.out").c_str(), std::ios::out);
		fOutputHMOM.setf(std::ios::scientific);

		// Labels
		{
			unsigned int count = 1;

			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "t[s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "x[cm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "T[K]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "v[cm/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "m[kg/m2/s]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "rho[kg/m3]", count);

			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "fv[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "n[#/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "dp[nm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "dc[nm]", count);
			OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, "np[nm]", count);

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "M(" + label.str() + ")[mol/m3]";
				OpenSMOKE::PrintTagOnASCIILabel(20, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "Sall(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "Snuc(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "Sgro(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "Soxi(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "Scon(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(25, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaTot(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaDis(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaDisSS(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaDisSL(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaDisLL(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaCon(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaConSS(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaConSL(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
			{
				std::stringstream label;	label << j;
				std::string title = "ScoaConLL(" + label.str() + ")[mol/m3/s]";
				OpenSMOKE::PrintTagOnASCIILabel(30, fOutputHMOM, title, count);
			}

			fOutputHMOM << std::endl;
		}
		
		// Local residence time
		Eigen::VectorXd tau;
		ResidenceTime(tau);

		// Indices of relevant species jOH
		const int jOH = thermodynamicsMap_.IndexOfSpecies("OH") - 1;
		const int jH = thermodynamicsMap_.IndexOfSpecies("H") - 1;
		const int jH2O = thermodynamicsMap_.IndexOfSpecies("H2O") - 1;
		const int jH2 = thermodynamicsMap_.IndexOfSpecies("H2") - 1;
		const int jC2H2 = thermodynamicsMap_.IndexOfSpecies("C2H2") - 1;
		const int jO2 = thermodynamicsMap_.IndexOfSpecies("O2") - 1;
		std::vector<int> jPAH(hmom_->pah_species().size());
		for (unsigned int i = 0; i<hmom_->pah_species().size(); i++)
			jPAH[i] = thermodynamicsMap_.IndexOfSpecies(hmom_->pah_species()[i]) - 1;

		// Loop over all the cells
		for (int i = 0; i < grid_.Np(); i++)
		{
			// Gas-phase properties
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << tau(i);
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << grid_.x()(i) * 100.;
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << T_(i);
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << m_(i) / rho_(i)*100.;
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << m_(i);
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << rho_(i);

			// Analysis of soot
			const double ctot = P_ / PhysicalConstants::R_J_kmol / T_[i];	// [kmol/m3]
			const double mass_fraction_OH = Y_[i](jOH);
			const double mass_fraction_H = Y_[i](jH);
			const double conc_OH = ctot*X_[i](jOH);
			const double conc_H = ctot*X_[i](jH);
			const double conc_H2O = ctot*X_[i](jH2O);
			const double conc_H2 = ctot*X_[i](jH2);
			const double conc_C2H2 = ctot*X_[i](jC2H2);
			const double conc_O2 = ctot*X_[i](jO2);
			double conc_PAH = 0.;
			for (unsigned int j = 0; j<hmom_->pah_species().size(); j++)
				conc_PAH += ctot*X_[i](jPAH[j]);

			hmom_->SetNormalizedMoments(hmom_M_[i](0), hmom_M_[i](1), hmom_M_[i](2), hmom_M_[i](3));
			hmom_->SetTemperatureAndPressure(T_[i], P_);
			hmom_->SetMassFractions(mass_fraction_OH, mass_fraction_H);
			hmom_->SetConcentrations("kmol/m3", conc_OH, conc_H, conc_H2O, conc_H2, conc_C2H2, conc_O2, conc_PAH);
			hmom_->SetViscosity(SutherlandViscosity(T_[i]));
			hmom_->CalculateSourceMoments();
			
			// Soot properties
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_->SootVolumeFraction();
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_->SootParticleNumberDensity();
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_->SootParticleDiameter()*1e9;
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_->SootCollisionParticleDiameter()*1e9;
			fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_->SootNumberOfPrimaryParticles();
			
			// Soot moments
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(20) << hmom_M_[i](j);

			// Source terms (overall)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(25) << hmom_->sources()(j);

			// Source terms (nucleation)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(25) << hmom_->sources_nucleation()(j);

			// Source terms (surface growth)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(25) << hmom_->sources_growth()(j);

			// Source terms (oxidation)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(25) << hmom_->sources_oxidation()(j);

			// Source terms (condensation)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(25) << hmom_->sources_condensation()(j);

			// Source terms (coagulation)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_overall()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_discrete()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_discrete_ss()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_discrete_sl()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_discrete_ll()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_continous()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_continous_ss()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_continous_sl()(j);

			// Source terms (coagulation discrete)
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				fOutputHMOM << std::scientific << std::setprecision(9) << std::setw(30) << hmom_->sources_coagulation_continous_ll()(j);
			
			fOutputHMOM << std::endl;
		}

		fOutputHMOM.close();
	}

	void OpenSMOKE_PremixedLaminarFlame1D::PrintOnTheFlyPostProcessing()
	{
		// Local residence time
		Eigen::VectorXd tau;
		ResidenceTime(tau);

		// Output files
		on_the_fly_post_processing_->PrepareOutputFiles();
		for (int i = 0; i < grid_.Np(); i++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble omega(thermodynamicsMap_.NumberOfSpecies());
			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				omega[j+1] = Y_[i](j);
			on_the_fly_post_processing_->WriteOnFile(tau(i), grid_.x()[i], 0., 0., T_(i), P_, omega);
		}
		on_the_fly_post_processing_->CloseOutputFiles();
	}

	void OpenSMOKE_PremixedLaminarFlame1D::MinimumUnknownsVector(double* v)
	{
		const double minimum_temperature = 100.;		// [K]
		const double minimum_mass_flow_rate = 1.e-9;	// [kg/s]
		const double zero = 0.;

		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = zero;

				v[count++] = minimum_temperature;
				v[count++] = minimum_mass_flow_rate;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = zero;

				v[count++] = minimum_temperature;
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = zero;
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				v[count++] = minimum_temperature;
				v[count++] = minimum_mass_flow_rate;
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = zero;
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = zero;

				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = zero;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = zero;
				v[count++] = minimum_temperature;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = zero;
			}
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::MaximumUnknownsVector(double* v)
	{
		const double maximum_temperature = 10000.;		// [K]
		const double maximum_mass_flow_rate = 10000.;	// [kg/s]
		const double one = 1.;
		const double big = 1.e16;

		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = one;

				v[count++] = maximum_temperature;
				v[count++] = maximum_mass_flow_rate;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = one;

				v[count++] = maximum_temperature;
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = one;
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				v[count++] = maximum_temperature;
				v[count++] = maximum_mass_flow_rate;
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = big;
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = one;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = big;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = one;
				v[count++] = maximum_temperature;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = big;
			}
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::UnknownsVector(double* v)
	{
		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = Y_[i](j);

				v[count++] = T_(i);
				v[count++] = m_(i);
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = Y_[i](j);

				v[count++] = T_(i);
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = Y_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				v[count++] = T_(i);
				v[count++] = m_(i);
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = hmom_M_[i](j);
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = Y_[i](j);
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = hmom_M_[i](j);
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					v[count++] = Y_[i](j);
				v[count++] = T_(i);
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					v[count++] = hmom_M_[i](j);
			}
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::CorrectedUnknownsVector(const double* v)
	{
		if (type_ == SIMULATION_TYPE_YTM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = v[count++];

				T_(i) = v[count++];
				m_(i) = v[count++];
			}

			// Boundary
			m_inlet_ = m_(0);
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = v[count++];

				T_(i) = v[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = v[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				T_(i) = v[count++];
				m_(i) = v[count++];
			}

			// Boundary
			m_inlet_ = m_(0);
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = v[count++];
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = v[count++];
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = v[count++];
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			unsigned int count = 0;
			for (int i = 0; i < grid_.Np(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					Y_[i](j) = v[count++];
				T_(i) = v[count++];
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					hmom_M_[i](j) = v[count++];
			}
		}

		// Properties and fluxes
		Properties();
		DiffusionFluxes();
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SetAlgebraicDifferentialEquations()
	{
		if (type_ == SIMULATION_TYPE_YTM)
		{
			id_equations_.resize((thermodynamicsMap_.NumberOfSpecies() + 2)*grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = true;

				id_equations_[count++] = true;
				id_equations_[count++] = false;
			}

			// Outlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT)
		{
			id_equations_.resize((thermodynamicsMap_.NumberOfSpecies() + 1)*grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = true;
				id_equations_[count++] = true;
			}

			// Outlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_Y)
		{
			id_equations_.resize(thermodynamicsMap_.NumberOfSpecies()*grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = true;
			}

			// Outlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_TM)
		{
			id_equations_.resize(2*grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				id_equations_[count++] = true;
				id_equations_[count++] = false;
			}

			// Outlet boundary
			{
				id_equations_[count++] = false;
				id_equations_[count++] = false;
			}
		}
		else if (type_ == SIMULATION_TYPE_HMOM)
		{
			id_equations_.resize(hmom_->n_moments() * grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				id_equations_[count++] = false;

			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = true;

			// Outlet boundary
			for (unsigned int j = 0; j < hmom_->n_moments(); j++)
				id_equations_[count++] = true;
		}
		else if (type_ == SIMULATION_TYPE_Y_HMOM)
		{
			id_equations_.resize( (hmom_->n_moments()+ thermodynamicsMap_.NumberOfSpecies()) * grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = false;
			}
			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = true;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = true;
			}

			// Outlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = true;
			}
		}
		else if (type_ == SIMULATION_TYPE_YT_HMOM)
		{
			id_equations_.resize((hmom_->n_moments() + 1 + thermodynamicsMap_.NumberOfSpecies()) * grid_.Np());

			unsigned int count = 0;

			// Inlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = false;
			}
			// Internal points
			for (int i = 1; i < grid_.Ni(); i++)
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = true;
				id_equations_[count++] = true;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = true;
			}

			// Outlet boundary
			{
				for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
					id_equations_[count++] = false;
				id_equations_[count++] = false;
				for (unsigned int j = 0; j < hmom_->n_moments(); j++)
					id_equations_[count++] = true;
			}
		}

		unsigned int n_algebraic = std::count_if(id_equations_.begin(), id_equations_.end(), std::bind2nd(std::equal_to<bool>(), false));
		unsigned int n_differential = std::count_if(id_equations_.begin(), id_equations_.end(), std::bind2nd(std::equal_to<bool>(), true));

		algebraic_equations_.resize(n_algebraic);
		differential_equations_.resize(n_differential);

		int count_differential = 0;
		int count_algebraic = 0;
		for (unsigned int i = 0; i < id_equations_.size(); i++)
		{
			if (id_equations_[i] == true)  differential_equations_(count_differential++) = i;
			if (id_equations_[i] == false) algebraic_equations_(count_algebraic++) = i;
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::AlgebraicDifferentialVector(double* v)
	{
		int count = 0;
		for (unsigned int i = 0; i < id_equations_.size(); i++)
		{
			if (id_equations_[i] == true)  v[count++] = 1.;
			if (id_equations_[i] == false) v[count++] = 0.;
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::CorrectDifferentialEquations(double* upv, double* resv)
	{
		for (int i = 0; i < differential_equations_.size(); i++)
		{
			const int k = differential_equations_[i];
			resv[k] -= upv[k];
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::CorrectAlgebraicEquations(double* yp)
	{
		for (int i = 0; i < algebraic_equations_.size(); i++)
		{
			const int k = algebraic_equations_[i];
			yp[k] = 0.;
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::Update(const std::vector<Eigen::VectorXd>& phi)
	{
		MemoryAllocation();

		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
			for (int j = 0; j < grid_.Np(); j++)
					Y_[j](i) = phi[i](j);
		
		T_ = phi[thermodynamicsMap_.NumberOfSpecies()];
		m_ = phi[thermodynamicsMap_.NumberOfSpecies() + 1];

		m_inlet_ = m_(0);

		Properties();
		DiffusionFluxes();
	}



	// Print XML Files
	void OpenSMOKE_PremixedLaminarFlame1D::PrintXMLFile(const std::string file_name)
	{
		const unsigned int n_additional = 8;
		
		std::ofstream fXML;
		fXML.open(file_name.c_str(), std::ios::out);
		fXML.setf(std::ios::scientific);

		fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
		fXML << "<opensmoke version=\"0.1a\">" << std::endl;

		fXML << "<Type> Flame1D </Type>" << std::endl;

		fXML << "<additional>" << std::endl;
		fXML << n_additional << std::endl;
		fXML << "axial-coordinate [cm] 2" << std::endl;
		fXML << "temperature [K] 3" << std::endl;
		fXML << "pressure [Pa] 4" << std::endl;
		fXML << "mol-weight [kg/kmol] 5" << std::endl;
		fXML << "density [kg/m3] 6" << std::endl;
		fXML << "heat-release [W/m3] 7" << std::endl;
		fXML << "axial-velocity [m/s] 8" << std::endl;
		fXML << "mass-flow-rate [kg/m2/s] 9" << std::endl;
		fXML << "</additional>" << std::endl;

		fXML << "<t-p-mw>" << std::endl;
		fXML << "1 2 3" << std::endl;
		fXML << "</t-p-mw>" << std::endl;

		fXML << "<mass-fractions>" << std::endl;
		fXML << thermodynamicsMap_.NumberOfSpecies() << std::endl;
		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
			fXML << thermodynamicsMap_.NamesOfSpecies()[i] << " " << thermodynamicsMap_.MW(i) << " " << n_additional + (i+1) << std::endl;
		fXML << "</mass-fractions>" << std::endl;

		fXML << "<profiles>" << std::endl;
		for (int i = 0; i < grid_.Np(); i++)
		{
			fXML << 1.e2*grid_.x()(i) << " ";
			fXML << T_(i) << " ";
			fXML << P_ << " ";
			fXML << mw_(i) << " ";
			fXML << rho_(i) << " ";
			fXML << Q_(i) << " ";
			fXML << m_(i)/rho_(i) << " ";
			fXML << m_(i) << " ";

			for (unsigned int j = 0; j < thermodynamicsMap_.NumberOfSpecies(); j++)
				fXML << Y_[i](j) << " ";
			fXML << std::endl;
		}
		fXML << "</profiles>" << std::endl;

		fXML << "<profiles-size>" << std::endl;
		fXML << grid_.Np() << " " << thermodynamicsMap_.NumberOfSpecies() + n_additional << std::endl;
		fXML << "</profiles-size>" << std::endl;
		fXML << "</opensmoke>" << std::endl;
		fXML.close();
	}



	int OpenSMOKE_PremixedLaminarFlame1D::SolveInitialDAE(DaeSMOKE::DaeSolver_Parameters& dae_parameters, const double tEnd)
	{
		int flag = -1;

		// Solving the DAE system
		if (dae_parameters.type() != DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_BZZDAE)
		{
			if (dae_parameters.sparse_linear_algebra() == false)
				flag = DaeSMOKE::Solve_Band_OpenSMOKEppDae<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_OpenSMOKEpp_Premixed_FLAMESPEED>(this, dae_parameters, 0., tEnd);
			else
				flag = DaeSMOKE::Solve_Sparse_OpenSMOKEppDae<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_OpenSMOKEpp_Premixed_FLAMESPEED>(this, dae_parameters, 0., tEnd);
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else // if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_BZZDAE)
			flag = DaeSMOKE::Solve_TridiagonalBlock_BzzDae<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_Premixed_FLAMESPEED>(this, dae_object_, dae_parameters, 0., tEnd);
		#endif

		return flag;
	}

	int OpenSMOKE_PremixedLaminarFlame1D::SolveDAE(DaeSMOKE::DaeSolver_Parameters& dae_parameters, const double tEnd)
	{
		int flag = -1;

		if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_OPENSMOKEPP)
		{
			if (dae_parameters.sparse_linear_algebra() == false)
				flag = DaeSMOKE::Solve_Band_OpenSMOKEppDae<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_OpenSMOKEpp_Premixed_FLAMESPEED>(this, dae_parameters, 0., tEnd);
			else
				flag = DaeSMOKE::Solve_Sparse_OpenSMOKEppDae<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_OpenSMOKEpp_Premixed_FLAMESPEED>(this, dae_parameters, 0., tEnd);
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_BZZDAE)
			flag = DaeSMOKE::Solve_TridiagonalBlock_BzzDae<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyDaeSystem_Premixed_FLAMESPEED>(this, dae_object_, dae_parameters, 0., tEnd);
		#endif
		#if OPENSMOKE_USE_SUNDIALS == 1
		else if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_IDA)
			flag = DaeSMOKE::Solve_Band_Ida<OpenSMOKE_PremixedLaminarFlame1D>(this, dae_parameters, 0., tEnd);
		#endif
		#if OPENSMOKE_USE_DASPK == 1
		else if (dae_parameters.type() == DaeSMOKE::DaeSolver_Parameters::DAE_INTEGRATOR_DASPK)
			flag = DaeSMOKE::Solve_Band_Daspk<OpenSMOKE_PremixedLaminarFlame1D>(this, dae_parameters, 0., tEnd);
		#endif

		return flag;
	}

	int OpenSMOKE_PremixedLaminarFlame1D::SolveNLS(NlsSMOKE::NonLinearSolver_Parameters& nls_parameters)
	{
		int flag = -1;

		if (nls_parameters.type() == NlsSMOKE::NonLinearSolver_Parameters::NLS_SOLVER_OPENSMOKEPP)
		{
			if (nls_parameters.sparse_linear_algebra() == false)
				flag = NlsSMOKE::Solve_Band_OpenSMOKEppNls<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyNlsSystem_OpenSMOKEpp_Premixed_FLAMESPEED>(this, nls_parameters);
			else
				flag = NlsSMOKE::Solve_Sparse_OpenSMOKEppNls<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyNlsSystem_OpenSMOKEpp_Premixed_FLAMESPEED>(this, nls_parameters);
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else if (nls_parameters.type() == NlsSMOKE::NonLinearSolver_Parameters::NLS_SOLVER_BZZNLS)
			flag = NlsSMOKE::Solve_TridiagonalBlock_BzzNls<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyNlsSystem_Premixed_FLAMESPEED>(this, nls_object_, nls_parameters);
		#endif
		#if OPENSMOKE_USE_SUNDIALS == 1
		else if (nls_parameters.type() == NlsSMOKE::NonLinearSolver_Parameters::NLS_SOLVER_KINSOL)
			flag = NlsSMOKE::Solve_Band_KinSol<OpenSMOKE_PremixedLaminarFlame1D>(this, nls_parameters);
		#endif

		return flag;
	}

	int OpenSMOKE_PremixedLaminarFlame1D::SolveFalseTransient(NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		int flag = -1;

		if (false_transient_parameters.type() == NlsSMOKE::FalseTransientSolver_Parameters::FALSETRANSIENT_SOLVER_OPENSMOKEPP)
		{
			if (false_transient_parameters.sparse_linear_algebra() == false)
				flag = NlsSMOKE::Solve_Band_OpenSMOKEppFalseTransient<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyFalseTransientSystem_OpenSMOKEpp_Premixed_FLAMESPEED>(this, false_transient_parameters);
			else
				flag = NlsSMOKE::Solve_Sparse_OpenSMOKEppFalseTransient<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyFalseTransientSystem_OpenSMOKEpp_Premixed_FLAMESPEED>(this, false_transient_parameters);
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else if (false_transient_parameters.type() == NlsSMOKE::FalseTransientSolver_Parameters::FALSETRANSIENT_SOLVER_BZZNLS)
			flag = NlsSMOKE::Solve_Band_BzzNlsFalseTransient<OpenSMOKE_PremixedLaminarFlame1D, OpenSMOKE_Flame1D_MyFalseTransientSystem_Premixed_FLAMESPEED>(this, nls_object_, false_transient_parameters);
		#endif
		#if OPENSMOKE_USE_SUNDIALS == 1
		else if (false_transient_parameters.type() == NlsSMOKE::FalseTransientSolver_Parameters::FALSETRANSIENT_SOLVER_KINSOL)
			flag = NlsSMOKE::Solve_Band_KinSolFalseTransient<OpenSMOKE_PremixedLaminarFlame1D>(this, false_transient_parameters);
		#endif

		return flag;
	}

	int OpenSMOKE_PremixedLaminarFlame1D::InitialSolutionFlameSpeed(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
		NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
		NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		// Solution from scratch:
		// 1. Only mass fractions
		//    a. DAE: OpenSMOKE++ (default) || BzzDae
		//    b. NLS: User choice
		// 2. Only temperature and mass flow rate
		//    a. DAE: OpenSMOKE++ (default) || BzzDae
		//    b. NLS: User choice
		// 3. Complete
		//    a. DAE: OpenSMOKE++ (default) || BzzDae
		//    b. NLS: User choice

		// Step 1
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "                Initial Solution: Step 1 (Y)              " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_Y);

			int flag = -1;

			// Solve the DAE system
			SolveInitialDAE(dae_parameters, timeFixedTemperature_);

			// Then solve the NLS
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);
			
			// Write solution on ASCII file
			Print((output_folder_ / "Solution.initial.Y.out").string().c_str());
		}

		// Step 2
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "              Initial Solution: Step 2 (T+M)              " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_TM);

			int flag = -1;

			// Solve the DAE system
			SolveInitialDAE(dae_parameters, timeFixedTemperature_);

			// Then solve the NLS
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);

			// Write solution on ASCII file
			Print((output_folder_ / "Solution.initial.TM.out").string().c_str());
		}

		// Step 3
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "             Initial Solution: Step 3 (Y+T+M)             " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_YTM);

			int flag = -1;

			// Solve the DAE system
			SolveInitialDAE(dae_parameters, timeComplete_);

			// Then solve the NLS
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);

			// Write solution on ASCII file	
			Print((output_folder_ / "Solution.initial.YTM.out").string().c_str() );

			// Polimi soot
			if (is_polimi_soot_ == true)	
				PrintSoot(output_folder_);

			// On the fly post-processing
			if (is_on_the_fly_post_processing_ == true)	
				PrintOnTheFlyPostProcessing();
		}

		return 1;
	}

	int OpenSMOKE_PremixedLaminarFlame1D::InitialSolutionBurnerStabilized(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
		NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
		NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		// Solution from scratch:
		// 1. Only mass fractions
		//    a. DAE: OpenSMOKE++ (default) || BzzDae
		//    b. NLS: User choice
		// 2. Mass fractions and temperature
		//    a. DAE: OpenSMOKE++ (default) || BzzDae
		//    b. NLS: User choice

		// Step 1
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "                Initial Solution: Step 1 (Y)              " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_Y);

			int flag = -1;

			// Solve the DAE system
			SolveInitialDAE(dae_parameters, timeFixedTemperature_);

			// Then solve the NLS
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);

			// Write solution on ASCII file
			Print((output_folder_ / "Solution.initial.Y.out").string().c_str());
		}

		// Step 2
		if (is_fixed_temperature_profile_ == false)
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "             Initial Solution: Step 3 (Y+T+M)             " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_YT);

			int flag = -1;

			// Solve the DAE system
			SolveInitialDAE(dae_parameters, timeComplete_);

			// Then solve the NLS
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);

			// Write solution on ASCII file	
			Print((output_folder_ / "Solution.initial.YT.out").string().c_str());

			// Polimi soot
			if (is_polimi_soot_ == true)
				PrintSoot(output_folder_);

			// On the fly post-processing
			if (is_on_the_fly_post_processing_ == true)
				PrintOnTheFlyPostProcessing();
		}

		return 1;
	}

	int OpenSMOKE_PremixedLaminarFlame1D::FixedTemperatureSolution(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
		NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
		NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		int flag = -1;

		std::cout << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "              Solution at fixed temperature               " << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		SetType(SIMULATION_TYPE_Y);
		
		if (use_dae_solver_ == true)
		{
			flag = SolveDAE(dae_parameters, timeFixedTemperature_);

			if (flag < 0)
			{
				// Then solve the non linear system
				flag = SolveNLS(nls_parameters);

				// In case of failure
				if (flag < 0)
				{
					// Solve the false transient
					SolveFalseTransient(false_transient_parameters);

					// Then solves the non linear system
					flag = SolveNLS(nls_parameters);
				}
			}
		}
		else  // The NLS is solved by definition
		{
			// Then solve the non linear system
			flag = SolveNLS(nls_parameters);

			// In case of failure
			if (flag < 0)
			{
				// Solve the false transient
				SolveFalseTransient(false_transient_parameters);

				// Then solves the non linear system
				flag = SolveNLS(nls_parameters);
			}
		}

		return 1;
	}

	int OpenSMOKE_PremixedLaminarFlame1D::CompleteSolution(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
		NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
		NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		int flag = -1;

		// Step 2
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "               Preparation: Step 2 (T+M)                  " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_TM);

			// Solve first the DAE system
			flag = SolveDAE(dae_parameters, timeFixedComposition_);

			// Solving the NLS
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);
		}

		std::cout << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "                   Solution (complete)                    " << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		SetType(SIMULATION_TYPE_YTM);

		if (use_dae_solver_ == true)
		{
			// Solve first the DAE system
			flag = SolveDAE(dae_parameters, timeComplete_);

			// Then solve the non linear system
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);
		}
		else  // The NLS is solved by definition
		{
			// Then solve the non linear system
			flag = SolveNLS(nls_parameters);

			// In case of failure
			if (flag < 0)
			{
				// Solve the false transient
				SolveFalseTransient(false_transient_parameters);

				// Then solves the non linear system
				flag = SolveNLS(nls_parameters);
			}
		}

		CheckForAdiabaticity();
		AtomicAnalysis();
		CheckForInlet();

		return 1;
	}

	int OpenSMOKE_PremixedLaminarFlame1D::SolveBurnerStabilized(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
		NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
		NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		int flag = -1;

		// Step 1
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "               Preparation: Step 1 (Y)                    " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_Y);

			const double tEnd = 1.;

			// Solve first the DAE system
			flag = SolveDAE(dae_parameters, tEnd);

			// Then solve the non linear system
			if (use_nls_solver_ == true)
				flag = SolveNLS(nls_parameters);
		}

		// Step 2
		if (is_fixed_temperature_profile_ == false)
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "                   Solution (complete)                    " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			SetType(SIMULATION_TYPE_YT);

			if (use_dae_solver_ == true)
			{
				// Solve first the DAE system
				flag = SolveDAE(dae_parameters, timeComplete_);

				// Then solve the non linear system
				if (use_nls_solver_ == true)
					flag = SolveNLS(nls_parameters);
			}
			else  // The NLS is solved by definition
			{
				// Solve directly the non linear system
				flag = SolveNLS(nls_parameters);

				// In case of failure
				if (flag < 0)
				{
					// Solve the false transient
					SolveFalseTransient(false_transient_parameters);

					// Then solves the non linear system
					flag = SolveNLS(nls_parameters);
				}
			}
		}

		CheckForAdiabaticity();
		AtomicAnalysis();
		CheckForInlet();

		return 1;
	}

	int OpenSMOKE_PremixedLaminarFlame1D::SolveHMOMFromExistingSolution(	
						OpenSMOKE::HMOM& hmom,
						DaeSMOKE::DaeSolver_Parameters& dae_parameters,
						NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
						NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters )
	{
		// HMOM pointer
		is_hmom_soot_ = true;
		hmom_ = &hmom;

		// -----------------------------------------------------------------------------------
		//				                   Memory allocation
		// -----------------------------------------------------------------------------------
		
		dhmom_M_over_dx_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			dhmom_M_over_dx_[i].resize(hmom_->n_moments());

		dhmom_M_over_dt_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			dhmom_M_over_dt_[i].resize(hmom_->n_moments());

		hmom_M_.resize(grid_.Np());
		for (int i = 0; i < grid_.Np(); i++)
			hmom_M_[i].resize(hmom_->n_moments());

		// -----------------------------------------------------------------------------------
		//				                   Initial values
		// -----------------------------------------------------------------------------------
		for (int i = 0; i < grid_.Np(); i++)
			hmom_M_[i].setConstant(1e-12);

		// -----------------------------------------------------------------------------------

		std::cout << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "                       Solving HMOM                       " << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		if (hmom_->PAHConsumption() == false)
		{
			SetType(SIMULATION_TYPE_HMOM);
		}
		else
		{
			if (is_fixed_temperature_profile_ == true)
				SetType(SIMULATION_TYPE_Y_HMOM);
			else
				SetType(SIMULATION_TYPE_YT_HMOM);
		}

		const double tEnd = 10.;

		// Solve first the DAE system
		int flag = SolveDAE(dae_parameters, tEnd);

		// Print on file
		PrintHMOM(output_folder_);

		// Print main solution on XML file
		PrintXMLFile((output_folder_ / "Output.xml").string().c_str());

		// Write current solution
		Print((output_folder_ / "Solution.final.out").string().c_str());

		// On the fly post-processing
		if (is_on_the_fly_post_processing_ == true)
			PrintOnTheFlyPostProcessing();

		return flag;
	}

	OpenSMOKE::Adapter_Grid1D_Status OpenSMOKE_PremixedLaminarFlame1D::RefineGrid(const unsigned int count)
	{
		OpenSMOKE::Adapter_Grid1D_Status refinement_status;

		std::cout << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "                  Grid refinement: " << count               << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;

		std::vector<Eigen::VectorXd> phi(thermodynamicsMap_.NumberOfSpecies() + 2);
		for (unsigned int i = 0; i < phi.size(); i++)
			phi[i].resize(grid_.Np());

		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
		for (int j = 0; j < grid_.Np(); j++)
			phi[i](j) = Y_[j](i);
		phi[thermodynamicsMap_.NumberOfSpecies()] = T_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 1] = m_;

		std::vector<Eigen::VectorXd> phi_new;
		refinement_status = grid_.Refine(phi, phi_new);

		// If the temperature profile has to be kept fixed, we need to perform
		// interpolation from the user-defined profile, not from the previous profile
		if (is_fixed_temperature_profile_ == true && refinement_status != OpenSMOKE::NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED)
		{
			std::cout << "Interpolating temperature from user-defined profile" << std::endl;
			fixed_temperature_profile_->Interpolate(grid_.x(), phi_new[thermodynamicsMap_.NumberOfSpecies()]);
		}

		// If the specific mass flow rate profile profile has to be kept fixed, we need to perform
		// interpolation from the user-defined profile, not from the previous profile
		if (is_fixed_specific_mass_flow_rate_profile_ == true && refinement_status != OpenSMOKE::NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED)
		{
			std::cout << "Interpolating specific mass flow rate profile from user-defined profile" << std::endl;
			fixed_specific_mass_flow_rate_profile_->Interpolate(grid_.x(), phi_new[thermodynamicsMap_.NumberOfSpecies() + 1]);
		}

		// In case new points have been added
		if (refinement_status != OpenSMOKE::NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED)
		{
			Update(phi_new);
		}
		
		return refinement_status;
	}

	OpenSMOKE::Adapter_Grid1D_Status OpenSMOKE_PremixedLaminarFlame1D::RefineGrid(const double xA, const double xB, unsigned int count)
	{
		std::cout << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "                  Grid refinement (local): " << count       << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;

		std::vector<Eigen::VectorXd> phi(thermodynamicsMap_.NumberOfSpecies() + 2);
		for (unsigned int i = 0; i < phi.size(); i++)
			phi[i].resize(grid_.Np());

		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
		for (int j = 0; j < grid_.Np(); j++)
			phi[i](j) = Y_[j](i);
		phi[thermodynamicsMap_.NumberOfSpecies()] = T_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 1] = m_;

		std::vector<Eigen::VectorXd> phi_new;
		grid_.Refine(xA, xB, phi, phi_new);
		Update(phi_new);

		return OpenSMOKE::REGRID_SUCCESS;
	}

	OpenSMOKE::Adapter_Grid1D_Status OpenSMOKE_PremixedLaminarFlame1D::Doubling(unsigned int count)
	{
		std::cout << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "                  Grid doubling: " << count << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;

		std::vector<Eigen::VectorXd> phi(thermodynamicsMap_.NumberOfSpecies() + 2);
		for (unsigned int i = 0; i < phi.size(); i++)
			phi[i].resize(grid_.Np());

		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
		for (int j = 0; j < grid_.Np(); j++)
			phi[i](j) = Y_[j](i);
		phi[thermodynamicsMap_.NumberOfSpecies()] = T_;
		phi[thermodynamicsMap_.NumberOfSpecies() + 1] = m_;

		std::vector<Eigen::VectorXd> phi_new;
		grid_.Double(phi, phi_new);
		Update(phi_new);

		return OpenSMOKE::REGRID_SUCCESS;
	}


	int OpenSMOKE_PremixedLaminarFlame1D::SolveFlameSpeedFromScratch(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
		NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
		NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		std::ofstream fMonitoring;
		fMonitoring.open((output_folder_ / "monitoring.out").string().c_str(), std::ios::out);
		fMonitoring.setf(std::ios::scientific);
		
		// Initial solution

		InitialSolutionFlameSpeed(dae_parameters, nls_parameters, false_transient_parameters);
				
		// Monitoring
		fMonitoring << grid_.Np() << "\t" << m_(0) / rho_(0)*100. << std::endl;

		// Loop
		unsigned int max_refinement_attempts_ = 1000;
		for (unsigned int k = 1; k <= max_refinement_attempts_; k++)
		{	
			OpenSMOKE::Adapter_Grid1D_Status refinement_status;
			refinement_status = RefineGrid(k);


			if (refinement_status == OpenSMOKE::NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED)
			{
				std::cout << "\n criteria satisified \n";
				break;
			}
			else
			{
			
				// Without energy and mass flow rate
				FixedTemperatureSolution(dae_parameters, nls_parameters, false_transient_parameters);

				// With energy and mass flow rates
				CompleteSolution(dae_parameters, nls_parameters, false_transient_parameters);
			}

			// Update statistics
			//norm();

			// Write current solution
			{
				Print((output_folder_ / "Solution.current.out").string().c_str());

				// Polimi soot
				if (is_polimi_soot_ == true)
					PrintSoot(output_folder_);

				// On the fly post-processing
				if (is_on_the_fly_post_processing_ == true)
					PrintOnTheFlyPostProcessing();

				PrintXMLFile((output_folder_ / "Output.xml").string().c_str());
			}

			// Monitoring
			fMonitoring << grid_.Np() << "\t" << m_(0) / rho_(0)*100. << std::endl;

			if (refinement_status == OpenSMOKE::MAXIMUM_NUMBER_POINTS)
				break;
		}
		
		fMonitoring.close();
		
		// Print final solution
		{
			Print((output_folder_ / "Solution.final.out").string().c_str());
			
			// Polimi soot
			if (is_polimi_soot_ == true)
				PrintSoot(output_folder_);

			// On the fly post-processing
			if (is_on_the_fly_post_processing_ == true)
				PrintOnTheFlyPostProcessing();
			
			PrintXMLFile((output_folder_ / "Output.xml").string().c_str());
			
		}
		
		// Sensitivity Analysis
		//if (sensitivity_analysis() == true)
			// Sarebbe da mettere un avviso o un errore per modificare l'input file!!! TODO
			//SensitivityAnalysis();
		
		return 1;
		
	}

	int OpenSMOKE_PremixedLaminarFlame1D::SolveFlameSpeedFromScratchForOptimization(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
		NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
		NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters)
	{
		// Initial solution
		InitialSolutionFlameSpeed(dae_parameters, nls_parameters, false_transient_parameters);

		// Loop
		unsigned int max_refinement_attempts_ = 1000;
		for (unsigned int k = 1; k <= max_refinement_attempts_; k++)
		{
			OpenSMOKE::Adapter_Grid1D_Status refinement_status;
			refinement_status = RefineGrid(k);

			if (refinement_status == OpenSMOKE::NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED)
			{
				break;
			}
			else
			{
				// Without energy and mass flow rate
				FixedTemperatureSolution(dae_parameters, nls_parameters, false_transient_parameters);

				// With energy and mass flow rates
				CompleteSolution(dae_parameters, nls_parameters, false_transient_parameters);
			}

			// Update statistics
			norm();

			if (refinement_status == OpenSMOKE::MAXIMUM_NUMBER_POINTS)
				break;
		}

		return 1;
	}

	int OpenSMOKE_PremixedLaminarFlame1D::SolveBurnerStabilizedFromScratch(DaeSMOKE::DaeSolver_Parameters& dae_parameters,
		NlsSMOKE::NonLinearSolver_Parameters& nls_parameters,
		NlsSMOKE::FalseTransientSolver_Parameters& false_transient_parameters
		)
	{
		// Initial solution
		InitialSolutionBurnerStabilized(dae_parameters, nls_parameters, false_transient_parameters);

		// Loop
		unsigned int max_refinement_attempts_ = 1000;
		for (unsigned int k = 1; k <= max_refinement_attempts_; k++)
		{
			OpenSMOKE::Adapter_Grid1D_Status refinement_status;
			refinement_status = RefineGrid(k);

			if (refinement_status == OpenSMOKE::NO_ADDED_POINTS_BECAUSE_CRITERIA_SATISFIED)
			{
				break;
			}
			else
			{
				SolveBurnerStabilized(dae_parameters, nls_parameters, false_transient_parameters);
			}

			// Update statistics
			norm();

			// Write current solution
			{
				Print((output_folder_ / "Solution.current.out").string().c_str());

				// Polimi soot
				if (is_polimi_soot_ == true)
					PrintSoot(output_folder_);

				// On the fly post-processing
				if (is_on_the_fly_post_processing_ == true)
					PrintOnTheFlyPostProcessing();

				PrintXMLFile((output_folder_ / "Output.xml").string().c_str());
			}

			if (refinement_status == OpenSMOKE::MAXIMUM_NUMBER_POINTS)
				break;
		}

		// Print final solution
		{
			Print((output_folder_ / "Solution.final.out").string().c_str());
			
			// Polimi soot
			if (is_polimi_soot_ == true)
				PrintSoot(output_folder_);

			// On the fly post-processing
			if (is_on_the_fly_post_processing_ == true)
				PrintOnTheFlyPostProcessing();

			PrintXMLFile((output_folder_ / "Output.xml").string().c_str());
		}

		// Sensitivity Analysis (TODO)
		if (sensitivity_analysis() == true)
			// Sarebbe da mettere warning su errore o input come sopra TODO
			//SensitivityAnalysis();
		
		std::cout<<"The address of wall T function "<< wall_heat_temperature_profile_ << std::endl;
		return 1;
	}

	// Atomic analysis
	void OpenSMOKE_PremixedLaminarFlame1D::AtomicAnalysis()
	{
		// Inlet stream
		std::vector<double> sum_inlet(thermodynamicsMap_.elements().size());
		{
			double MWmix;
			aux_Y.CopyFrom(Y_inlet_.data());
			thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), MWmix, aux_Y.GetHandle());

			for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
			{
				sum_inlet[j] = 0.;
				for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
					sum_inlet[j] += thermodynamicsMap_.atomic_composition()(i, j) * aux_X[i + 1];
			}
		}
		
		// Outlet stream
		std::vector<double> sum_final(thermodynamicsMap_.elements().size());
		for (unsigned int j = 0; j < thermodynamicsMap_.elements().size(); j++)
		{
			sum_final[j] = 0.;
			for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
				sum_final[j] += thermodynamicsMap_.atomic_composition()(i, j) * X_[grid_.Ni()](i);
		}

		const double moles_inlet = m_(0) / mw_(0);
		const double moles_final = m_(grid_.Ni()) / mw_(grid_.Ni());

		// Write on the screen
		{

			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout	<< std::setw(16) << std::left << "Atomic balances"
						<< std::setw(15) << std::left << "inlet"
						<< std::setw(15) << std::left << "outlet"
						<< std::setw(15) << std::left << "error(%)"
						<< std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;

			for (unsigned int j = 0; j<thermodynamicsMap_.elements().size(); j++)
			if (sum_inlet[j] > 0.)
			{
				std::string label = thermodynamicsMap_.elements()[j] + "[kmol/s]";
				std::cout << std::setw(16) << std::left << label
					<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
					<< sum_inlet[j] * moles_inlet
					<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
					<< sum_final[j] * moles_final
					<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
					<< std::fabs((sum_final[j] * moles_final) - (sum_inlet[j] * moles_inlet)) / (sum_inlet[j] * moles_inlet) * 100.
					<< std::endl;
			}

			std::cout << std::endl;
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::CheckForAdiabaticity()
	{
		// Velocities and kinetic energy
		const double v_inlet = m_(0) / rho_(0);
		const double v_outlet = m_(grid_.Ni()) / rho_(grid_.Ni());
		const double Ek_inlet = 0.50*m_(0)*v_inlet*v_inlet;
		const double Ek_outlet = 0.50*m_(grid_.Ni())*v_outlet*v_outlet;

		// Enthalpy: inlet 
		double H_inlet;
		{
			// Set thermodynamic map
			thermodynamicsMap_.SetTemperature(T_(0));
			thermodynamicsMap_.SetPressure(P_);

			// Mole fractions and molecular weight
			double MWmix;
			aux_Y.CopyFrom(Y_inlet_.data());
			thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), MWmix, aux_Y.GetHandle());

			// Enthalpy [W]
			H_inlet = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(aux_X.GetHandle());
			H_inlet *= m_(0)/MWmix;
		}

		// Enthalpy: outlet
		double H_outlet;
		{
			// Set thermodynamic map
			thermodynamicsMap_.SetTemperature(T_(grid_.Ni()));
			thermodynamicsMap_.SetPressure(P_);

			// Mole fractions
			aux_X.CopyFrom(X_[grid_.Ni()].data());

			// Enthalpy [W]
			H_outlet = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(aux_X.GetHandle());
			H_outlet *= m_(grid_.Ni()) / mw_(grid_.Ni());
		}

		// Conductive fluxes [W]
		double Q_inlet  = -lambda_(0)*			(T_(1)-T_(0)) / 
												(grid_.x()(1)-grid_.x()(0));						// implicitly assume area equal to 1 m2
		double Q_outlet = -lambda_(grid_.Ni())*	(T_(grid_.Ni()) - T_(grid_.Ni() - 1)) /
												(grid_.x()(grid_.Ni()) - grid_.x()(grid_.Ni()-1));	// implicitly assume area equal to 1 m2

		// Write on the screen
		{
			std::cout	<< "----------------------------------------------------------" << std::endl;
			std::cout	<< std::setw(16) << std::left << "Energy balances"
						<< std::setw(15) << std::left << "inlet"
						<< std::setw(15) << std::left << "outlet"
						<< std::setw(15) << std::left << "error(%)"
						<< std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;

			std::cout	<< std::setw(16) << std::left << "Ek[W]"
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< Ek_inlet
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< Ek_outlet
						<< std::endl;

			std::cout	<< std::setw(16) << std::left << "Q[W]"
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< Q_inlet
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< Q_outlet
						<< std::endl;


			std::cout	<< std::setw(16) << std::left << "H[W]"
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< H_inlet
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< H_outlet
						<< std::endl;

			std::cout	<< std::setw(16) << std::left << "Etot[W]"
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< H_inlet + Ek_inlet + Q_inlet
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< H_outlet + Ek_outlet
						<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
						<< std::fabs((H_inlet + Ek_inlet + Q_inlet) - (H_outlet + Ek_outlet+Q_outlet)) / (H_inlet + Ek_inlet + Q_inlet) * 100.
						<< std::endl;

			std::cout << std::endl;
		}

	}

	void OpenSMOKE_PremixedLaminarFlame1D::CheckForInlet()
	{
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << std::setw(16) << std::left << "Checking inlet"
			<< std::setw(16) << std::left << "nominal"
			<< std::setw(16) << std::left << "current"
			<< std::setw(16) << std::left << "error(%)"
			<< std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;

		// Check for species
		for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
		{
			if (Y_inlet_(i) > 0.)
				std::cout	<< std::setw(16) << std::left << thermodynamicsMap_.NamesOfSpecies()[i]
							<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
							<< Y_inlet_(i)
							<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
							<< Y_[0](i)
							<< std::setw(15) << std::left << std::scientific << std::setprecision(3)
							<< std::fabs(Y_inlet_(i) - Y_[0](i)) / Y_inlet_(i) * 100.
							<< std::endl;
		}

		std::cout << std::endl;
	}

	OpenSMOKE::Adapter_Grid1D_Status OpenSMOKE_PremixedLaminarFlame1D::Regrid()
	{
		if (grid_.Np() > static_cast<int>(grid_.grid_adapter().regrid_points()))
		{
			std::cout << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "                 Regridding...                            " << std::endl;
			std::cout << "----------------------------------------------------------" << std::endl;

			std::vector<Eigen::VectorXd> phi(thermodynamicsMap_.NumberOfSpecies() + 2);
			for (unsigned int i = 0; i < phi.size(); i++)
				phi[i].resize(grid_.Np());

			for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
			for (int j = 0; j < grid_.Np(); j++)
				phi[i](j) = Y_[j](i);
			phi[thermodynamicsMap_.NumberOfSpecies()] = T_;
			phi[thermodynamicsMap_.NumberOfSpecies() + 1] = m_;

			std::vector<Eigen::VectorXd> phi_new;
			grid_.Regrid(thermodynamicsMap_.NumberOfSpecies(), phi, phi_new);
			Update(phi_new);
		}
		else
		{
			OpenSMOKE::FatalErrorMessage("Regrid operation cannot be performed because the original number of points is too small.");
		}

		return OpenSMOKE::REGRID_SUCCESS;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::norm()
	{
		double* y_ = new double[NumberOfEquations()];
		double* f_ = new double[NumberOfEquations()];

		UnknownsVector(y_);
		Equations(0., y_, f_);

		double norm2 = 0.;
		for (int i = 0; i < NumberOfEquations(); i++)
			norm2 += f_[i] * f_[i];
		norm2 = std::sqrt(norm2);

		double yp_mean = 0.;
		for (int i = 0; i < NumberOfEquations(); i++)
			yp_mean += std::fabs(f_[i]);
		yp_mean /= NumberOfEquations();

		std::cout << " n=" << NumberOfEquations() / BlockDimensions() << std::scientific << "  ||f||=" << norm2 << "  |y'|=" << yp_mean << std::endl;
	}

	void OpenSMOKE_PremixedLaminarFlame1D::InitializeFromBackupFile(const boost::filesystem::path path_file, const bool use_userdefined_grid)
	{
		std::vector<double> x_old;
		std::vector<double> T_old;
		std::vector<double> P_old;
		std::vector<double> m_old;
		std::vector < std::vector<double> > omega_old;
		ReadFromBackupFile(path_file, thermodynamicsMap_, x_old, T_old, P_old, m_old, omega_old);

		// Interpolation (if needed)
		if (use_userdefined_grid == true)
		{
			// Grid: use the user-defined grid
			{
				std::cout << " * Building the new mesh from user-defined mesh..." << std::endl;

				std::cout << "   Position of fixed point (mm): " << grid_.x_fixed_point()*1e3 << std::endl;
				std::cout << "   Index of fixed point:         " << grid_.fixed_point() << "/" << grid_.Np() << std::endl;

				// Only in case of a burner stabilized flame
				if (solver_type_ == SOLVER_TYPE_BURNERSTABILIZED)
					grid_.ResetFixedPoint();
			}

			// Setup 
			{
				std::cout << " * Building the first-guess solution from backup data and user-defined mesh..." << std::endl;

				MemoryAllocation();

				// Pressure
				// The user is free to choose a different pressure (useful to investigate the role of pressure)

				// Profiles
				{
					// Memory allocation
					std::vector<Eigen::VectorXd> v_old(thermodynamicsMap_.NumberOfSpecies() + 1);
					std::vector<Eigen::VectorXd> v_new(thermodynamicsMap_.NumberOfSpecies() + 1);
					for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies() + 1; i++)
					{
						v_old[i].resize(x_old.size());
						v_new[i].resize(grid_.Np());
					}

					// Construct the input matrix
					Eigen::VectorXd x_old_(x_old.size());
					for (unsigned int j = 0; j < x_old.size(); j++)
					{
						x_old_(j) = x_old[j];
						v_old[thermodynamicsMap_.NumberOfSpecies()](j) = T_old[j];
						for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
							v_old[i](j) = omega_old[i][j];
					}

					// Linear interpolation
					std::cout << "Linear interpolation from backup mesh to user-defined mesh" << std::endl;
					linear_interpolation(x_old_, v_old, grid_.x(), v_new);

					// Construct the input matrix
					for (unsigned int j = 0; j < grid_.Np(); j++)
					{
						T_(j) = v_new[thermodynamicsMap_.NumberOfSpecies()](j);
						for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
							Y_[j](i) = v_new[i](j);
					}

					std::cout << "Interpolated temperature profile on the user-defined mesh" << std::endl;
					for (unsigned int j = 0; j < grid_.Np(); j++)
						std::cout << grid_.x()(j) << "\t" << T_(j) << std::endl;
					std::cout << std::endl;
				}
			}
		}
		else
		{
			// Grid: use the last grid available from the backup file
			{
				std::cout << " * Building the new mesh from backup data..." << std::endl;

				Eigen::VectorXd x_eigen(x_old.size());
				std::cout << "   Position of fixed point (mm): " << grid_.x_fixed_point()*1e3 << std::endl;
				std::cout << "   Index of fixed point:         " << grid_.fixed_point() << "/" << grid_.Np() << std::endl;

				// The grid used is the old one
				for (unsigned int i = 0; i < x_old.size(); i++)
					x_eigen(i) = x_old[i];

				// Only in case of a burner stabilized flame
				if (solver_type_ == SOLVER_TYPE_BURNERSTABILIZED)
					grid_.ResetFixedPoint();

				// Update the grid
				grid_.Update(x_eigen);
			}

			// Setup 
			{
				std::cout << " * Building the first-guess solution from backup data..." << std::endl;

				MemoryAllocation();

				// Pressure
				// The user is free to choose a different pressure (useful to investigate the role of pressure)

				// Temperature
				for (int i = 0; i < grid_.Np(); i++)
					T_(i) = T_old[i];

				// Composition
				for (unsigned int i = 0; i < thermodynamicsMap_.NumberOfSpecies(); i++)
				{
					for (int j = 0; j < grid_.Np(); j++)
						Y_[j](i) = omega_old[i][j];
				}
			}
		}
		
		// Complete setup
		{
			// First guess mass flow rate (i.e. velocity)
			{
				// Molecular weight
				aux_Y.CopyFrom(Y_[0].data());
				thermodynamicsMap_.MoleFractions_From_MassFractions(aux_X.GetHandle(), mw_(0), aux_Y.GetHandle());

				// Density
				const double rho = P_*mw_(0) / PhysicalConstants::R_J_kmol / T_inlet_;

				// Mass flow rate 
				if (solver_type_ == SOLVER_TYPE_BURNERSTABILIZED)
				{
					// Velocity is the one provided through the input file
					if (is_v_inlet_ == true)	m_inlet_ = rho * v_inlet_;
					if (is_v_inlet_ == false)	v_inlet_ = m_inlet_ / rho;
					m_.setConstant(m_inlet_);

					std::cout << "   Old mass flow rate [kg/m2/s]: " << m_old[0] << std::endl;
					std::cout << "   New mass flow rate [kg/m2/s]: " << m_inlet_ << std::endl;
					std::cout << "   Old velocity [cm/s]:          " << (m_old[0] / rho)*100. << std::endl;
					std::cout << "   New velocity [cm/s]:          " << v_inlet_*100. << std::endl;
				}
				else
				{
					m_inlet_ = m_old[0];
					m_.setConstant(m_inlet_);
					v_inlet_ = m_inlet_ / rho;
				}
			}

			Properties();
			DiffusionFluxes();

			fixed_T_ = T_(grid_.fixed_point());
		}

		// Messages on the screen
		if (P_ != P_old[0])
		{
			std::cout << "* WARNING: The backup pressure was equal to: " << P_old[0] << " Pa, while the current pressure is: " << P_ << " Pa" << std::endl;
			std::cout << "           Changing the pressure is allowed, but be sure that this is really what you want to do!" << std::endl;
		}

		if (T_inlet_ != T_old[0])
		{
			std::cout << "* WARNING: The backup temperature was equal to: " << T_old[0] << " K, while the current temperature is: " << T_ << " K" << std::endl;
			std::cout << "           Changing the temperature is allowed, but be sure that this is really what you want to do!" << std::endl;
		}
	}

	void OpenSMOKE_PremixedLaminarFlame1D::SparsityPattern(std::vector<unsigned int>& rows, std::vector<unsigned int>& cols)
	{
		std::vector<unsigned int> rows_single;
		std::vector<unsigned int> cols_single;
		OpenSMOKE::SparsityPatternTridiagonal(NumberOfEquations() / BlockDimensions(), rows_single, cols_single);
		OpenSMOKE::SparsityPatternBlock(NumberOfEquations() / BlockDimensions(), BlockDimensions(), rows_single, cols_single, rows, cols);
	}

	double OpenSMOKE_PremixedLaminarFlame1D::SutherlandViscosity(const double T)
	{
		return 1.716e-5*std::pow(T / 273.15, 1.5)*(273.15 + 110.4) / (T + 110.4);	// [kg/m/s]
	}
}
