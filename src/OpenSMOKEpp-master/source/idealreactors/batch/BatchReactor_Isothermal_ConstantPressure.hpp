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

#include "BatchReactor_OdeInterfaces.h"
#include "math/PhysicalConstants.h"
#include "math/native-ode-solvers/MultiValueSolver"

namespace OpenSMOKE
{
	#if OPENSMOKE_USE_BZZMATH == 1
	ODESystem_BzzOde_BatchReactor* ptOde_BatchReactor_Isothermal_ConstantPressure;
	void ODE_Print_BatchReactor_Isothermal_ConstantPressure(BzzVector &Y, double t)
	{
		ptOde_BatchReactor_Isothermal_ConstantPressure->MyPrint(Y,t);
	}
	#endif

	BatchReactor_Isothermal_ConstantPressure::BatchReactor_Isothermal_ConstantPressure(	
								OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								OpenSMOKE::ODE_Parameters& ode_parameters,
								OpenSMOKE::BatchReactor_Options& batch_options,
								OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa,
								OpenSMOKE::OnTheFlyCEMA& on_the_fly_cema,
								OpenSMOKE::OnTheFlyPostProcessing& on_the_fly_post_processing,
								OpenSMOKE::IgnitionDelayTimes_Analyzer& idts_analyzer,
								OpenSMOKE::PolimiSoot_Analyzer& polimi_soot_analyzer,
								const double V0, const double T0, const double P0, 
								const OpenSMOKE::OpenSMOKEVectorDouble& omega0) :

		BatchReactor(thermodynamicsMap, kineticsMap, ode_parameters, batch_options, on_the_fly_ropa, on_the_fly_cema, on_the_fly_post_processing, idts_analyzer, polimi_soot_analyzer)

	{
		type_ = BATCH_REACTOR_ISOTHERMAL_CONSTANTP;

		iteration_ = 0;
		counter_file_video_ = 0;
		counter_file_ASCII_ = 0;
		counter_file_XML_ = 0;
		counter_sensitivity_XML_ = 0;

		V0_ = V0;
		T0_ =T0 ;
		P0_ = P0;
		omega0_ = omega0;

		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();
		NE_ = NC_;

		MemoryAllocation();

		MW0_ = thermodynamicsMap_.MolecularWeight_From_MassFractions(omega0_.GetHandle());
		rho0_ = P0_ * MW0_ / (PhysicalConstants::R_J_kmol*T0_);
		mass_ = rho0_ * V0_;

		T_		= T0_;
		P_		= P0_;
		V_      = V0_;
		omega_	= omega0_;
		MW_     = MW0_;
		rho_    = rho0_;
		QR_     = 0.;

		thermodynamicsMap_.MoleFractions_From_MassFractions(x0_.GetHandle(), MW0_, omega0_.GetHandle());
		thermodynamicsMap_.SetPressure(P0_);
		thermodynamicsMap_.SetTemperature(T0_);
		H0_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x0_.GetHandle());
		U0_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(x0_.GetHandle());

		OpenAllFiles();
	}

	int BatchReactor_Isothermal_ConstantPressure::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
	{
		// Recover mass fractions
		#ifdef CHECK_MASSFRACTIONS
		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = std::max(y[i], 0.);
		#else
			omega_ = y;
		#endif

		// Calculates the volume and the concentrations of species
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
		cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
		Product(cTot_, x_, &c_);
		rho_ = cTot_*MW_;
		V_ = mass_/rho_;

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

		// Recovering residuals
		for (unsigned int i=1;i<=NC_;++i)	
			dy[i] = thermodynamicsMap_.MW(i-1)*R_[i]/rho_;

		return 0;
	}

	int BatchReactor_Isothermal_ConstantPressure::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		iteration_++;

		//AB
		if (iteration_ == 0)
			time_vector.resize(0);
		
		time_vector.push_back(t);

		// Calculates the pressure
		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = y[i];

		thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
		species_matrix.push_back(x_);
		//AB

		if (batch_options_.verbose_video() == true)
		{
			// Video output
			if (iteration_%batch_options_.n_step_video() == 1 || batch_options_.n_step_video() == 1)
			{
				counter_file_video_++;
				if (counter_file_video_ % 100 == 1)
				{
					std::cout << std::endl;
					std::cout << std::setw(6) << std::left << "#Step";
					std::cout << std::setw(16) << std::left << "Time[s]";
					std::cout << std::setw(10) << std::left << "T[K]";
					std::cout << std::setw(10) << std::left << "P[atm]";
					std::cout << std::endl;
				}
				std::cout << std::fixed << std::setw(6) << std::left << iteration_;
				std::cout << std::scientific << std::setw(16) << std::setprecision(6) << std::left << t;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << T_;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << P_ / 101325.;
				std::cout << std::endl;
			}

			// ASCII file output
			if (batch_options_.verbose_output() == true)
			{
				if (batch_options_.verbose_ascii_file() == true)
				{
					if (iteration_%batch_options_.n_step_file() == 1 || batch_options_.n_step_file() == 1)
					{
						counter_file_ASCII_++;
						PrintFinalStatus(fASCII_, t);

						if (on_the_fly_post_processing_.is_active() == true)
							on_the_fly_post_processing_.WriteOnFile(t, 0., 0., 0., T_, P_, omega_);
					}
				}

				// XML file output
				if (batch_options_.verbose_xml_file() == true)
				{
					if (iteration_%batch_options_.n_step_file() == 1 || batch_options_.n_step_file() == 1)
					{
						counter_file_XML_++;
						fXML_ << t << " ";
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
			}

			// Rate of Production Analysis (on the fly)
			if (on_the_fly_ropa_.is_active() == true)
				on_the_fly_ropa_.Analyze(fROPA_, iteration_, t, T_, P_, c_, omega_, omega0_);

			// Polimi Soot Analyzer (on the fly)
			if (polimi_soot_analyzer_.is_active() == true)
				PolimiSootAnalysis(t);
		}

		if (batch_options_.sensitivity_analysis() == true)
			SensitivityAnalysis(t, y);

		return 0;
	}

	void BatchReactor_Isothermal_ConstantPressure::Solve(const double t0, const double tf)
	{
		if (batch_options_.verbose_video() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << " Solving the batch reactor...                                                " << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
		}

		// Initial conditions
		y0_ = omega0_;

		// Print intial conditions
		{
			OpenSMOKE::OpenSMOKEVectorDouble dy0(y0_.Size());
			Equations(t0, y0_, dy0);
			Print(t0, y0_);
		}

		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			// Min and max values
			Eigen::VectorXd yMin(NE_); for (unsigned int i = 0; i < NE_; i++) yMin(i) = 0.;
			Eigen::VectorXd yMax(NE_); for (unsigned int i = 0; i < NE_; i++) yMax(i) = 1.;

			// Initial conditions
			Eigen::VectorXd y0_eigen(y0_.Size());
			y0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Create the solver
			typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE_BatchReactor> denseOde;
			typedef OdeSMOKE::MethodGear<denseOde> methodGear;
			OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
			ode_solver.SetReactor(this);

			// Set initial conditions
			ode_solver.SetInitialConditions(t0, y0_eigen);

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
			BzzVector yMax(NE_); yMax=1.; 
			
			// Initial conditions
			BzzVector y0_bzz(y0_.Size());
			y0_.CopyTo(y0_bzz.GetHandle());
			
			// Final solution
			BzzVector yf_bzz(y0_bzz.Size());

			ODESystem_BzzOde_BatchReactor odebatch(*this);
			BzzOdeStiffObject o(y0_bzz, t0, &odebatch);
			
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			o.SetTolAbs(ode_parameters_.absolute_tolerance());
			o.SetTolRel(ode_parameters_.relative_tolerance());

			ptOde_BatchReactor_Isothermal_ConstantPressure = &odebatch;
			o.StepPrint(ODE_Print_BatchReactor_Isothermal_ConstantPressure);

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
			SolveOpenSourceSolvers(t0, tf);
		}

		if (batch_options_.verbose_video() == true)
			ode_parameters_.Status(std::cout);

		FinalStatus(t0, tf);
		FinalSummary(batch_options_.output_path() / "FinalSummary.out", t0, tf);

		CloseAllFiles();
	}

	void BatchReactor_Isothermal_ConstantPressure::SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		if (iteration_ == 1)
		{
			// Writes the coefficients on file (only on request)
			if (iteration_%batch_options_.n_step_file() == 1 || batch_options_.n_step_file() == 1)		
			{
				counter_sensitivity_XML_++;
				for (unsigned int k=0;k<indices_of_sensitivity_species_.size();k++)
				{
					for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
						fSensitivityChildXML_[k] << 0. << " ";		
					fSensitivityChildXML_[k] << std::endl;	
				}
			}
		}
		else
		{
			// Scaling factors
			for(unsigned int j=1;j<=NC_;j++)
				scaling_Jp_[j] = thermodynamicsMap_.MW(j-1)/rho_;

			// Calculates the current Jacobian
			NumericalJacobian(t, y, Jnum_);

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
				V_ = mass_/rho_;
			}

			// Calculates the current sensitivity coefficients
			sensitivityMap_->Calculate(t, T_, P_, c_, Jnum_, scaling_Jp_);

			// Writes the coefficients on file (only on request)
			if (iteration_%batch_options_.n_step_file() == 1 || batch_options_.n_step_file() == 1)		
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

	void BatchReactor_Isothermal_ConstantPressure::ChemicalExplosiveModeAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		OpenSMOKE::ErrorMessage("BatchReactor_Isothermal_ConstantPressure", "ChemicalExplosiveModeAnalysis is not yet available for Adiabtic User Defined Batch Reactor");
	}

	void BatchReactor_Isothermal_ConstantPressure::SparseAnalyticalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::SparseMatrix<double> &J)
	{
		OpenSMOKE::ErrorMessage("BatchReactor_Isothermal_ConstantPressure", "SparseAnalyticalJacobian is not yet available for Isothermal Constant Pressure Batch Reactor");
	}

	void BatchReactor_Isothermal_ConstantPressure::DenseAnalyticalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::MatrixXd &J)
	{
		OpenSMOKE::ErrorMessage("BatchReactor_Isothermal_ConstantPressure", "SparseAnalyticalJacobian is not yet available for Isothermal Constant Pressure Batch Reactor");
	}
}
