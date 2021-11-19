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

#include "PerfectlyStirredReactor_OdeInterfaces.h"
#include "math/PhysicalConstants.h"
#include "math/native-ode-solvers/MultiValueSolver"

namespace OpenSMOKE
{
	#if OPENSMOKE_USE_BZZMATH == 1
	ODESystem_BzzOde_PerfectlyStirredReactor* ptOde_PerfectlyStirredReactor_NonIsothermal_ConstantPressure;
	void ODE_Print_PerfectlyStirredReactor_NonIsothermal_ConstantPressure(BzzVector &Y, double t)
	{
		ptOde_PerfectlyStirredReactor_NonIsothermal_ConstantPressure->MyPrint(Y,t);
	}
	#endif

	PerfectlyStirredReactor_NonIsothermal_ConstantPressure::PerfectlyStirredReactor_NonIsothermal_ConstantPressure(	
								OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								OpenSMOKE::ODE_Parameters& ode_parameters,
								OpenSMOKE::PerfectlyStirredReactor_Options& psr_options,
								OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa,
								OpenSMOKE::OnTheFlyPostProcessing& on_the_fly_post_processing,
								OpenSMOKE::PolimiSoot_Analyzer& polimi_soot_analyzer,
								const double T0, const double P0, const OpenSMOKE::OpenSMOKEVectorDouble& omega0,
								const double TInlet, const double PInlet, const OpenSMOKE::OpenSMOKEVectorDouble& omegaI,
								const double residence_time, const double volume, const double mass_flow_rate,
								const double global_thermal_exchange_coefficient,
								const double exchange_area, const double T_environment) :

	PerfectlyStirredReactor(thermodynamicsMap, kineticsMap, ode_parameters, psr_options, on_the_fly_ropa, on_the_fly_post_processing, polimi_soot_analyzer)

	{
		type_ = PERFECTLYSTIRRED_REACTOR_NONISOTHERMAL_CONSTANTP;

		T0_ = T0;
		P0_ = P0;
		omega0_ = omega0;
		
		TInlet_ = TInlet;
		PInlet_ = PInlet;
		omegaInlet_ = omegaI;

		tau0_ = residence_time;
		V0_ = volume;
		mass_flow_rate0_ = mass_flow_rate;

		global_thermal_exchange_coefficient_ = global_thermal_exchange_coefficient;
		exchange_area_ = exchange_area;
		T_environment_ = T_environment;

		iteration_ = 0;
		counter_file_video_ = 0;
		counter_file_ASCII_ = 0;
		counter_file_ASCII_ReactionRates_ = 0;
		counter_file_XML_ = 0;
		counter_sensitivity_XML_ = 0;

		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();
		NE_ = NC_+1;

		MemoryAllocation();

		// Constraints Analysis
		if (V0_<0.)
			constraint_type_ = MASSFLOWRATE_RESIDENCETIME;
		else if (tau0_ < 0.)
			constraint_type_ = MASSFLOWRATE_VOLUME;
		else if (mass_flow_rate0_ < 0.)
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
			V0_ = mass_flow_rate0_*tau0_/rho0_;
		else if (constraint_type_ == MASSFLOWRATE_VOLUME)
			tau0_ = rho0_*V0_/mass_flow_rate0_;
		else if (constraint_type_ == VOLUME_RESIDENCETIME)
			mass_flow_rate0_ = rho0_*V0_/tau0_;

		// Current values
		mass_flow_rate_ = mass_flow_rate0_;
		tau_ = tau0_;
		rho_ = rho0_;
		mass_ = rho_ * V0_;
		T_		= T0_;
		P_		= P0_;
		V_      = V0_;
		omega_	= omega0_;
		MW_     = MW0_;
		QR_     = 0.;

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

	int PerfectlyStirredReactor_NonIsothermal_ConstantPressure::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
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

		// Calculates the volume and the concentrations of species
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
		cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
		Product(cTot_, x_, &c_);
		rho_ = cTot_*MW_;
		
		// Constraints analysis
		if (constraint_type_ == MASSFLOWRATE_RESIDENCETIME)
			V_ = mass_flow_rate_*tau_/rho_;
		else if (constraint_type_ == MASSFLOWRATE_VOLUME)
			tau_ = rho_*V_/mass_flow_rate_;
		else if (constraint_type_ == VOLUME_RESIDENCETIME)
			mass_flow_rate_ = rho_*V_/tau_;

		// Calculates thermodynamic properties
		thermodynamicsMap_.SetTemperature(T_);
		thermodynamicsMap_.SetPressure(P_);
		
		const double CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
		CvMixMass_ = (CpMixMolar - PhysicalConstants::R_J_kmol) / MW_;
		CpMixMass_ = CpMixMolar / MW_;

		// Calculates kinetics
		kineticsMap_.SetTemperature(T_);
		kineticsMap_.SetPressure(P_);

		kineticsMap_.ReactionRates(c_.GetHandle());
		kineticsMap_.FormationRates(R_.GetHandle());
		QR_ = kineticsMap_.HeatRelease(R_.GetHandle());

		// Recovering residuals (mass fractions)
		for (unsigned int i=1;i<=NC_;++i)	
			dy[i] = (omegaInlet_[i]-omega_[i])/tau_ + thermodynamicsMap_.MW(i-1)*R_[i]/rho_;

		// Recovering residuals (temperature)
		const double Qexchange = global_thermal_exchange_coefficient_*exchange_area_*(T_environment_ - T_);		// [W]
		const double HStarI = Dot(xInlet_.Size(), xInlet_.GetHandle(), thermodynamicsMap_.Species_H_over_RT().data())*(PhysicalConstants::R_J_kmol*T_);	// [J/kmol]
		dy[NC_+1]  = ( (HInlet_-HStarI)/MWInlet_/tau_ + QR_/rho_ + Qexchange/(rho_*V_)) / CpMixMass_;

		return 0;
	}

	int PerfectlyStirredReactor_NonIsothermal_ConstantPressure::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
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
					std::cout << std::setw(14) << std::left << "m[g/s]";
					std::cout << std::endl;
				}
				std::cout << std::fixed << std::setw(7) << std::left << iteration_;
				std::cout << std::scientific << std::setw(15) << std::setprecision(6) << std::left << t;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(2) << T_;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(4) << P_/101325.;
				std::cout << std::setw(14) << std::left << std::scientific << std::setprecision(3) << V_*1.e6;
				std::cout << std::setw(14) << std::left << std::scientific << std::setprecision(3) << mass_flow_rate_*1000.;
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

					// On the fly post-processing
					if (on_the_fly_post_processing_.is_active() == true)
						on_the_fly_post_processing_.WriteOnFile(t, 0., 0., 0., T_, P_, omega_);
				}

					if(psr_options_.verbose_reactionrates_file() == true)
					{
						counter_file_ASCII_ReactionRates_++;
						rr_ = kineticsMap_.GiveMeReactionRates();
						fReactionRatesTime_.setf(std::ios::scientific);
						fReactionRatesTime_ << std::setw(20) << std::left << t;
						fReactionRatesTime_ << std::setw(20) << std::left << T0_;
						fReactionRatesTime_ << std::setw(20) << std::left << P0_;
						fReactionRatesTime_ << std::setw(20) << std::left << V0_;
						fReactionRatesTime_ << std::setw(20) << std::left << T_;
						fReactionRatesTime_ << std::setw(20) << std::left << P_;
						fReactionRatesTime_ << std::setw(20) << std::left << V_;
						fReactionRatesTime_ << std::setw(20) << std::left << rho_;
						fReactionRatesTime_ << std::setw(20) << std::left << MW_;
						for (unsigned int j = 0; j < kineticsMap_.NumberOfReactions(); j++)							
							fReactionRatesTime_ << std::setw(20) << std::left << rr_[j];		

						fReactionRatesTime_ << std::endl;
					}
			}
		}

		return 0;
	}

	void PerfectlyStirredReactor_NonIsothermal_ConstantPressure::PrintFinalStatus(std::ostream& fOutput)
	{
		fOutput.setf(std::ios::scientific);
		fOutput << std::setw(20) << std::left << tau_;
		fOutput << std::setw(20) << std::left << T0_;
		fOutput << std::setw(20) << std::left << P0_;
		fOutput << std::setw(20) << std::left << V0_;
		fOutput << std::setw(20) << std::left << T_;
		fOutput << std::setw(20) << std::left << P_;
		fOutput << std::setw(20) << std::left << V_;
		fOutput << std::setw(20) << std::left << rho_;
		fOutput << std::setw(20) << std::left << MW_;

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
	}

	void PerfectlyStirredReactor_NonIsothermal_ConstantPressure::PrintParametricFinalStatus(std::ostream& fOutput)
	{
		fOutput.setf(std::ios::scientific);
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

	void PerfectlyStirredReactor_NonIsothermal_ConstantPressure::PrintSolution()
	{
		if (psr_options_.verbose_output() == true)
		{
			// ASCII file output
			if (psr_options_.verbose_ascii_file() == true)
			{
				counter_file_ASCII_++;
				fASCIIFinal_.setf(std::ios::scientific);
				PrintFinalStatus(fASCIIFinal_);
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
					fXML_ << QR_ << " ";
					for (unsigned int i = 1; i <= NC_; i++)
						fXML_ << std::setprecision(12) << omega_[i] << " ";
					fXML_ << std::endl;
				}
			}

			// Rate of Production Analysis (on the fly)
			if (on_the_fly_ropa_.is_active() == true)
				on_the_fly_ropa_.Analyze(fROPA_, 0, tau_, T_, P_, c_, omega_, omegaInlet_);

			// Polimi Soot Analyzer (on the fly)
			if (polimi_soot_analyzer_.is_active() == true)
				PolimiSootAnalysis(tau_);
		}
	}

	void PerfectlyStirredReactor_NonIsothermal_ConstantPressure::Solve(const double tf)
	{
		if (psr_options_.verbose_video() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << " Solving the steady-state, adiabatic perfectly stirred reactor...            " << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
		}

		// Initial conditions
		for(unsigned int i=1;i<=NC_;++i)
			y0_[i] = omega0_[i];
		y0_[NC_+1] = T0_;

		// Print intial conditions
		{
			OpenSMOKE::OpenSMOKEVectorDouble dy0(y0_.Size());
			Equations(0., y0_, dy0);
			Print(0., y0_);
		}

		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			// Min and max values
			Eigen::VectorXd yMin(NE_); for (unsigned int i = 0; i<NE_; i++) yMin(i) = 0.;
			Eigen::VectorXd yMax(NE_); for (unsigned int i = 0; i<NE_; i++) yMax(i) = 1.;  yMax(NC_) = 50000.;

			// Initial conditions
			Eigen::VectorXd y0_eigen(y0_.Size());
			y0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Create the solver
			typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE_PerfectlyStirredReactor> denseOde;
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
			// Min and max values
			BzzVector yMin(NE_); yMin=0.; yMin[NC_+1] = 0.;
			BzzVector yMax(NE_); yMax=1.; yMax[NC_+1] = 50000.;
			
			// Initial conditions
			BzzVector y0_bzz(y0_.Size());
			y0_.CopyTo(y0_bzz.GetHandle());
			
			// Final solution
			BzzVector yf_bzz(y0_bzz.Size());

			ODESystem_BzzOde_PerfectlyStirredReactor odepsr(*this);
			BzzOdeStiffObject o(y0_bzz, 0., &odepsr);
			
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			o.SetTolAbs(ode_parameters_.absolute_tolerance());
			o.SetTolRel(ode_parameters_.relative_tolerance());

			ptOde_PerfectlyStirredReactor_NonIsothermal_ConstantPressure = &odepsr;
			o.StepPrint(ODE_Print_PerfectlyStirredReactor_NonIsothermal_ConstantPressure);

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

	void PerfectlyStirredReactor_NonIsothermal_ConstantPressure::SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		// Scaling factors
		for(unsigned int j=1;j<=NC_;j++)
			scaling_Jp[j] = thermodynamicsMap_.MW(j-1)/rho_;
		scaling_Jp[NC_+1] = 1./(rho_*CpMixMass_);

		// Calculates the current Jacobian
		NumericalJacobian(t, y, J);

		// Recover variables
		{
			#ifdef CHECK_MASSFRACTIONS
			for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = std::max(y[i], 0.);
			#else
			for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = y[i];
			#endif

			// Recover temperature
			T_ = y[NC_+1];

			// Calculates the pressure and the concentrations of species
			thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
			cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
			Product(cTot_, x_, &c_);
			rho_ = cTot_*MW_;

			// Constraints analysis
			if (constraint_type_ == MASSFLOWRATE_RESIDENCETIME)
				V_ = mass_flow_rate_*tau_/rho_;
			else if (constraint_type_ == MASSFLOWRATE_VOLUME)
				tau_ = rho_*V_/mass_flow_rate_;
			else if (constraint_type_ == VOLUME_RESIDENCETIME)
				mass_flow_rate_ = rho_*V_/tau_;
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

			for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
				fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << sensitivityMap_->sensitivity_coefficients()(NC_,j-1) << " ";
			fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << std::endl;
		}
	}
}
