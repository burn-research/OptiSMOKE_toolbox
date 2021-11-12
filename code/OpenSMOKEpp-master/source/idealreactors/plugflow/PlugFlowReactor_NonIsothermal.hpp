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

#include "PlugFlowReactor_OdeInterfaces.h"
#include "math/PhysicalConstants.h"
#include "math/native-ode-solvers/MultiValueSolver"

namespace OpenSMOKE
{
	#if OPENSMOKE_USE_BZZMATH == 1
	ODESystem_BzzOde_PlugFlowReactor* ptOde_PlugFlowReactor_NonIsothermal;
	void ODE_Print_PlugFlowReactor_NonIsothermal(BzzVector &Y, double t)
	{
		ptOde_PlugFlowReactor_NonIsothermal->MyPrint(Y,t);
	}
	#endif


	PlugFlowReactor_NonIsothermal::PlugFlowReactor_NonIsothermal(	
								OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								OpenSMOKE::ODE_Parameters& ode_parameters,
								OpenSMOKE::PlugFlowReactor_Options& plugflow_options,
								OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa,
								OpenSMOKE::OnTheFlyPostProcessing& on_the_fly_post_processing,
								OpenSMOKE::IgnitionDelayTimes_Analyzer& idts_analyzer,
								OpenSMOKE::PolimiSoot_Analyzer& polimi_soot_analyzer,
								const bool time_independent, 
								const bool constant_pressure, 
								const double v0, const double T0, const double P0, 
								const OpenSMOKE::OpenSMOKEVectorDouble& omega0,
                                const double global_thermal_exchange_coefficient,
                                const double cross_section_over_perimeter,
                                const double T_environment) :

		PlugFlowReactor(thermodynamicsMap, kineticsMap, ode_parameters, plugflow_options, on_the_fly_ropa, on_the_fly_post_processing, idts_analyzer, polimi_soot_analyzer)
	{
		type_ = PLUGFLOW_REACTOR_NONISOTHERMAL;

		time_independent_ = time_independent;
		constant_pressure_ = constant_pressure;
		iteration_ = 0;
		counter_file_video_ = 0;
		counter_file_ASCII_ = 0;
		counter_file_XML_ = 0;
		counter_sensitivity_XML_ = 0;

		v0_ = v0;
		T0_ = T0 ;
		P0_ = P0;
		omega0_ = omega0;
                
                global_thermal_exchange_coefficient_ = global_thermal_exchange_coefficient;
                cross_section_over_perimeter_ = cross_section_over_perimeter;
                T_environment_ = T_environment;

		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();
		NE_ = NC_+2;

		MemoryAllocation();

		MW0_ = thermodynamicsMap_.MolecularWeight_From_MassFractions(omega0_.GetHandle());
		rho0_ = P0_ * MW0_ / (PhysicalConstants::R_J_kmol*T0_);
		specificmassflowrate_ = rho0_ * v0_;

		T_		= T0_;
		P_		= P0_;
		v_      = v0_;
		omega_	= omega0_;
		MW_     = MW0_;
		QR_     = 0.;
		rho_    = rho0_;

		thermodynamicsMap_.MoleFractions_From_MassFractions(x0_.GetHandle(), MW0_, omega0_.GetHandle());
		thermodynamicsMap_.SetPressure(P0_);
		thermodynamicsMap_.SetTemperature(T0_);
		H0_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x0_.GetHandle());
		U0_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(x0_.GetHandle());

		OpenAllFiles();
	}

	int PlugFlowReactor_NonIsothermal::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
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

		// Calculates the volume and the concentrations of species
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
		if (constant_pressure_ == true)
		{
			cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
			Product(cTot_, x_, &c_);
			rho_ = cTot_*MW_;
			v_ = specificmassflowrate_/rho_;
		}
		else
		{
			P_ = specificmassflowrate_/v_ * PhysicalConstants::R_J_kmol * T_ / MW_;
			cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
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
		kineticsMap_.SetTemperature(T_);
		kineticsMap_.SetPressure(P_);
		
		//kineticsMap_.KineticConstants();
		kineticsMap_.ReactionRates(c_.GetHandle());
		kineticsMap_.FormationRates(R_.GetHandle());
		QR_ = kineticsMap_.HeatRelease(R_.GetHandle());

		// Recovering residuals
		for (unsigned int i=1;i<=NC_;++i)	
			dy[i] = thermodynamicsMap_.MW(i-1)*R_[i]/rho_/v_;

		// Exchange of heat with external environment
		double Qexchange = 0.;
		if (global_thermal_exchange_coefficient_ != 0.)
			Qexchange = global_thermal_exchange_coefficient_ / cross_section_over_perimeter_ * (T_environment_ - T_);		// [W/m3]

		// Energy equation
		const double v2 = v_*v_;
		const double Beta = 1.;
		const double delta = R_.SumElements()*MW_/rho_/v_;						// [1/m]
		const double T_coefficient = CpMixMass_ + v2/Beta/T_;					// [m2/s2/K]
		const double T_extra_terms = -v2/Beta*delta;							// [m/s2]
        dy[NC_+1] = ( (QR_+Qexchange)/rho_/v_ + T_extra_terms) / T_coefficient;	// [K/m]
	
		// Time-axial coordinate
		if (time_independent_ == true)	dy[NC_+2] = 1.;
		else                            dy[NC_+2] = 1./v_;

		if (time_independent_ == true)
			dy *= v_;

		return 0;
	}

	int PlugFlowReactor_NonIsothermal::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		// AB
		if (iteration_ == 0)
			time_vector.resize(0);

		iteration_++;

		time_vector.push_back(t);

		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = y[i];

		thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
		species_matrix.push_back(x_);

		if (plugflow_options_.verbose_video() == true)
		{
			// Video output
			if (iteration_%plugflow_options_.n_step_video() == 1 || plugflow_options_.n_step_video() == 1)
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
					std::cout << std::endl;
				}
				std::cout << std::fixed << std::setw(8) << std::left << iteration_;
				std::cout << std::scientific << std::setw(16) << std::setprecision(6) << std::left << t;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << T_;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << P_ / 101325.;
				std::cout << std::endl;
			}
		}
		
		if (plugflow_options_.verbose_output() == true)
		{
			// ASCII file output
			if (plugflow_options_.verbose_ascii_file() == true)
			{
				if (iteration_%plugflow_options_.n_step_file() == 1 || plugflow_options_.n_step_file() == 1)
				{
					counter_file_ASCII_++;
					PrintFinalStatus(fASCII_);

					if (on_the_fly_post_processing_.is_active() == true)
						on_the_fly_post_processing_.WriteOnFile(t, csi_, 0., 0., T_, P_, omega_);
				}
			}

			// XML file output
			if (plugflow_options_.verbose_xml_file() == true)
			{
				if (iteration_%plugflow_options_.n_step_file() == 1 || plugflow_options_.n_step_file() == 1)		
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
					fXML_ << QR_ << " ";
					fXML_ << v_ << " ";

					for (unsigned int i=1;i<=NC_;i++)
						fXML_ << std::setprecision(12) << omega_[i] << " ";
					fXML_ << std::endl;
				}
			}

			// Rate of Production Analysis (on the fly)
			if (on_the_fly_ropa_.is_active() == true)
				on_the_fly_ropa_.Analyze(fROPA_, iteration_, t, T_, P_, c_, omega_, omega0_);

			// Polimi Soot Analyzer (on the fly)
			if (polimi_soot_analyzer_.is_active() == true)
				PolimiSootAnalysis(t);
		}
		
		// Ignition delay times (on the fly)
		if (idts_analyzer_.is_active() == true)
			idts_analyzer_.Analyze(tau_, T_, P_, x_.GetHandle());

		// Sensitivity analysis (on the fly)
		if (plugflow_options_.sensitivity_analysis() == true)
			SensitivityAnalysis(t, y);

		return 0;
	}

	void PlugFlowReactor_NonIsothermal::PrintFinalStatus(std::ostream& fOutput)
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


	void PlugFlowReactor_NonIsothermal::PrintParametricFinalStatus(std::ostream& fOutput)
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

	void PlugFlowReactor_NonIsothermal::Solve(const double tf)
	{
		if (plugflow_options_.verbose_video() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << " Solving the plug flow reactor...                                            " << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
		}

		// Initial conditions
		for (unsigned int i=1;i<=NC_;i++)
			y0_[i] = omega0_[i];
		y0_[NC_+1] = T0_;
		y0_[NC_+2] = 0.;

		// Print intial conditions
		{
			OpenSMOKE::OpenSMOKEVectorDouble dy0(y0_.Size());
			Equations(0., y0_, dy0);
			Print(0., y0_);
		}

		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			// Min and max values
			Eigen::VectorXd yMin(NE_); for (unsigned int i = 0; i<NE_; i++) yMin(i) = 0.;  yMin(NC_) = 0.;
			Eigen::VectorXd yMax(NE_); for (unsigned int i = 0; i<NE_; i++) yMax(i) = 1.;  yMax(NC_) = 50000.; yMax(NC_ + 1) = 1.e16;

			// Initial conditions
			Eigen::VectorXd y0_eigen(y0_.Size());
			y0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Create the solver
			typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE_PlugFlowReactor> denseOde;
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
			BzzVector yMin(NE_); yMin=0.; 
			BzzVector yMax(NE_); yMax=1.; yMax[NC_+1] = 50000; yMax[NC_+2] = 1e16;
			
			// Initial conditions
			BzzVector y0_bzz(y0_.Size());
			y0_.CopyTo(y0_bzz.GetHandle());
			
			// Final solution
			BzzVector yf_bzz(y0_bzz.Size());

			ODESystem_BzzOde_PlugFlowReactor odeplugflow(*this);
			BzzOdeStiffObject o(y0_bzz, 0., &odeplugflow);
			
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			o.SetTolAbs(ode_parameters_.absolute_tolerance());
			o.SetTolRel(ode_parameters_.relative_tolerance());

			ptOde_PlugFlowReactor_NonIsothermal = &odeplugflow;
			o.StepPrint(ODE_Print_PlugFlowReactor_NonIsothermal);

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

		if (plugflow_options_.verbose_video() == true)
			ode_parameters_.Status(std::cout);

		FinalStatus(tf);
		FinalSummary(plugflow_options_.output_path() / "FinalSummary.out", tf);

		if (idts_analyzer_.is_active() == true)
			idts_analyzer_.PrintOnFile(plugflow_options_.output_path() / "IDT.out");

		CloseAllFiles();
	}

	void PlugFlowReactor_NonIsothermal::SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		if (iteration_ == 1)
		{
			// Writes the coefficients on file (only on request)
			if (iteration_%plugflow_options_.n_step_file() == 1 || plugflow_options_.n_step_file() == 1)		
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
					v_ = specificmassflowrate_/rho_;
				}
				else
				{
					P_ = specificmassflowrate_/v_ * PhysicalConstants::R_J_kmol * T_ / MW_;
					cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
					Product(cTot_, x_, &c_);
					rho_ = cTot_*MW_;
				}
			}

			// Calculates the current sensitivity coefficients
			sensitivityMap_->Calculate(t, T_, P_, c_, J, scaling_Jp);

			// Writes the coefficients on file (only on request)
			if (iteration_%plugflow_options_.n_step_file() == 1 || plugflow_options_.n_step_file() == 1)		
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
