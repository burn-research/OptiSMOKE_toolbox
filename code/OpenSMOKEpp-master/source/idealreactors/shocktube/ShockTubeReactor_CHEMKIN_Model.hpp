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

#include "ShockTubeReactor_OdeInterfaces.h"
#include "math/PhysicalConstants.h"
#include "math/native-ode-solvers/MultiValueSolver"

namespace OpenSMOKE
{
	#if OPENSMOKE_USE_BZZMATH == 1
	ODESystem_BzzOde_ShockTubeReactor* ptOde_ShockTubeReactor_CHEMKIN_Model;
	void ODE_Print_ShockTubeReactor_CHEMKIN_Model(BzzVector &Y, double t)
	{
		ptOde_ShockTubeReactor_CHEMKIN_Model->MyPrint(Y,t);
	}
	#endif

	ShockTubeReactor_CHEMKIN_Model::ShockTubeReactor_CHEMKIN_Model(	
																	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
																	OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
																	OpenSMOKE::ODE_Parameters& ode_parameters,
																	OpenSMOKE::ShockTubeReactor_Options& shocktube_options,
																	OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa,
																	OpenSMOKE::OnTheFlyPostProcessing& on_the_fly_post_processing,
																	OpenSMOKE::PolimiSoot_Analyzer& polimi_soot_analyzer,
																	const double v0, const double T0, const double P0, 
																	const OpenSMOKE::OpenSMOKEVectorDouble& omega0, const double rho1,
																	const double d0, const double lm) :

	ShockTubeReactor(thermodynamicsMap, kineticsMap, ode_parameters, shocktube_options, on_the_fly_ropa, on_the_fly_post_processing, polimi_soot_analyzer)

	{
		lm_ = lm;
		rho1_ = rho1;

		iteration_ = 0;
		counter_file_video_ = 0;
		counter_file_ASCII_ = 0;
		counter_file_XML_ = 0;
		counter_sensitivity_XML_ = 0;

		d0_ = d0;
		if (d0>0.)
			A0_ = PhysicalConstants::pi_over_4*d0_*d0_;
		else
			A0_ = 1.;
		v0_ = v0;
		T0_ = T0 ;
		P0_ = P0;
		omega0_ = omega0;

		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();
		NE_ = NC_+5;

		MemoryAllocation();

		thermodynamicsMap_.SetPressure(P0_);
		thermodynamicsMap_.SetTemperature(T0_);
		thermodynamicsMap_.MoleFractions_From_MassFractions(x0_.GetHandle(), MW0_, omega0_.GetHandle());
		H0_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x0_.GetHandle());
		U0_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(x0_.GetHandle());
		rho0_ = P0_ * MW0_ / (PhysicalConstants::R_J_kmol*T0_);

		d_		= d0_;
		A_		= A0_;
		rho_    = rho0_;
		T_		= T0_;
		P_		= P0_;
		v_      = v0_;
		omega_	= omega0_;
		MW_     = MW0_;
		QR_     = 0.;

		OpenAllFiles();
	}

	int ShockTubeReactor_CHEMKIN_Model::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
	{
		// Recover mass fractions
		#ifdef CHECK_MASSFRACTIONS
		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = std::max(y[i], 0.);
		#else
		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = y[i];
		#endif

		// Recover additional variables
		T_   = y[NC_+1];
		rho_ = y[NC_+2];
		v_   = y[NC_+3];
		z_   = y[NC_+4];
		tau_ = y[NC_+5];

		// Calculates the volume and the concentrations of species
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
		P_ = rho_*PhysicalConstants::R_J_kmol*T_/MW_;
		cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
		Product(cTot_, x_, &c_);

		// Calculates thermodynamic properties
		thermodynamicsMap_.SetTemperature(T_);
		thermodynamicsMap_.SetPressure(P_);
		
		const double CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
		CpMixMass_ = CpMixMolar / MW_;

		// Calculates kinetics
		kineticsMap_.SetTemperature(T_);
		kineticsMap_.SetPressure(P_);

		kineticsMap_.ReactionRates(c_.GetHandle());
		kineticsMap_.FormationRates(R_.GetHandle());
		QR_ = kineticsMap_.HeatRelease(R_.GetHandle());

		// Species equations
		{
			for (unsigned int i=1;i<=NC_;++i)	
				dy[i] = thermodynamicsMap_.MW(i-1)*R_[i]/rho_;
		}

		// Additional equations
		{
			const double v2 = v_*v_;
			const double v3 = v2*v_;
			const double CpT = CpMixMass_*T_;

			// Correct area
			double Acorrection = BoundaryLayerCorrection();

			// Density equation
			{
				const double sumRtilde = R_.SumElements();
				const double T1 = -P_*MW_*sumRtilde - P_/CpT*QR_;
				const double T2 = rho_*rho_*v3*(1.-P_/rho_/CpT)*Acorrection;
				const double T3 = P_*(1.+v2/CpT)-rho_*v2;
			
				dy[NC_+2] = (T1+T2)/T3;
			}

			// Energy equation
			dy[NC_+1] = ( v2*dy[NC_+2] + QR_ + rho_*v3*Acorrection) / (rho_*CpMixMass_);

			// Velocity equation
			dy[NC_+3] = -v_/rho_*dy[NC_+2] - v2*Acorrection;

			// Axial coordinate equation
			dy[NC_+4] = v_;

			// Laboratory time equation
			dy[NC_+5] = rho1_*A0_/rho_/A_;	
		}

		return 0;
	}

	int ShockTubeReactor_CHEMKIN_Model::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		iteration_++;

		if (shocktube_options_.verbose_output() == true)
		{
			// Video output
			if (iteration_%shocktube_options_.n_step_video() == 1 || shocktube_options_.n_step_video() == 1)
			{
				counter_file_video_++;
				if (counter_file_video_%100 == 1)
				{
					std::cout << std::endl;
					std::cout << std::setw(6)  << std::left << "#Step";
					std::cout << std::setw(16) << std::left << "Time[s]";
					std::cout << std::setw(10) << std::left << "T[K]";
					std::cout << std::setw(10) << std::left << "P[atm]";
					std::cout << std::endl;
				}
				std::cout << std::fixed << std::setw(6) << std::left << iteration_;
				std::cout << std::scientific << std::setw(16) << std::setprecision(6) << std::left << t;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << T_;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(5) << P_/101325.;
				std::cout << std::endl;
			}

			// ASCII file output
			if (shocktube_options_.verbose_ascii_file() == true)
			{
				if (iteration_%shocktube_options_.n_step_file() == 1 || shocktube_options_.n_step_file() == 1)
				{
					counter_file_ASCII_++;
					fASCII_.setf(std::ios::scientific);
					fASCII_ << std::setw(20) << std::left << t;
					fASCII_ << std::setw(20) << std::left << T_;
					fASCII_ << std::setw(20) << std::left << P_;
					fASCII_ << std::setw(20) << std::left << v_;
					fASCII_ << std::setw(20) << std::left << rho_;
					fASCII_ << std::setw(20) << std::left << MW_;
					fASCII_ << std::setw(20) << std::left << v_;
					fASCII_ << std::setw(20) << std::left << z_;
					fASCII_ << std::setw(20) << std::left << tau_;
					
					if (indices_of_output_species_.size() != 0)
					{
						for(unsigned int i=0;i<indices_of_output_species_.size();i++)
							fASCII_ << std::setw(widths_of_output_species_[i]) << std::left << x_[indices_of_output_species_[i]];
						for(unsigned int i=0;i<indices_of_output_species_.size();i++)
							fASCII_ << std::setw(widths_of_output_species_[i]) << std::left << omega_[indices_of_output_species_[i]];
					}
					else
					{
						for(unsigned int i=1;i<=NC_;i++)
							fASCII_ << std::setw(widths_of_output_species_[i-1]) << std::left << x_[i];
						for(unsigned int i=1;i<=NC_;i++)
							fASCII_ << std::setw(widths_of_output_species_[i-1]) << std::left << omega_[i];
					}
					fASCII_ << std::endl;

					// On the fly post-processing
					if (on_the_fly_post_processing_.is_active() == true)
						on_the_fly_post_processing_.WriteOnFile(t, z_, 0., 0., T_, P_, omega_);
				}
			}

			// XML file output
			if (shocktube_options_.verbose_xml_file() == true)
			{
				if (iteration_%shocktube_options_.n_step_file() == 1 || shocktube_options_.n_step_file() == 1)		
				{
					counter_file_XML_++;
					fXML_ << t << " ";
					fXML_ << T_ << " ";
					fXML_ << P_ << " ";
					fXML_ << MW_ << " ";
					fXML_ << rho_ << " ";
					fXML_ << QR_ << " ";
					fXML_ << v_ << " ";
					fXML_ << z_ << " ";
					fXML_ << tau_ << " ";
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

		if (shocktube_options_.sensitivity_analysis() == true)
			SensitivityAnalysis(t, y);

		return 0;
	}

	void ShockTubeReactor_CHEMKIN_Model::Solve(const double tf)
	{
		if (shocktube_options_.verbose_output() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << " Solving the shock tube reactor...                                                " << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
		}

		for(unsigned int i=1;i<=NC_;++i)
			y0_[i] = omega0_[i];
		y0_[NC_+1] = T0_;
		y0_[NC_+2] = rho0_;
		y0_[NC_+3] = v0_;
		y0_[NC_+4] = 0.;
		y0_[NC_+5] = 0.;


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
			Eigen::VectorXd yMax(NE_); for (unsigned int i = 0; i<NE_; i++) yMax(i) = 1.;
			yMax(NC_ + 0) = 50000.;
			yMax(NC_ + 1) = 10000.;
			yMax(NC_ + 2) = 50000.;
			yMax(NC_ + 3) = 1.e16;
			yMax(NC_ + 4) = 1.e16;

			// Initial conditions
			Eigen::VectorXd y0_eigen(y0_.Size());
			y0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Create the solver
			typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE_ShockTubeReactor> denseOde;
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
			BzzVector yMax(NE_); yMax=1.; 
			yMax[NC_+1] = 50000.;
			yMax[NC_+2] = 10000.;
			yMax[NC_+3] = 50000.;
			yMax[NC_+4] = 1.e16;
			yMax[NC_+5] = 1.e16;
			
			// Initial conditions
			BzzVector y0_bzz(y0_.Size());
			y0_.CopyTo(y0_bzz.GetHandle());
			
			// Final solution
			BzzVector yf_bzz(y0_bzz.Size());

			ODESystem_BzzOde_ShockTubeReactor odeshocktube(*this);
			BzzOdeStiffObject o(y0_bzz, 0., &odeshocktube);
			
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			o.SetTolAbs(ode_parameters_.absolute_tolerance());
			o.SetTolRel(ode_parameters_.relative_tolerance());

			ptOde_ShockTubeReactor_CHEMKIN_Model = &odeshocktube;
			o.StepPrint(ODE_Print_ShockTubeReactor_CHEMKIN_Model);

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

		if (shocktube_options_.verbose_output() == true)
		{
			ode_parameters_.Status(std::cout);
			FinalStatus(tf);
			FinalSummary(shocktube_options_.output_path() / "FinalSummary.out", tf);
		}

		CloseAllFiles();
	}

	void ShockTubeReactor_CHEMKIN_Model::SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		if (iteration_ == 1)
		{
			// Writes the coefficients on file (only on request)
			if (iteration_%shocktube_options_.n_step_file() == 1 || shocktube_options_.n_step_file() == 1)		
			{
				counter_sensitivity_XML_++;
				for (unsigned int k=0;k<indices_of_sensitivity_species_.size();k++)
				{
					for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
						fSensitivityChildXML_[k] << 0. << " ";		
					fSensitivityChildXML_[k] << std::endl;	
				}

				for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
					fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << 0. << " ";
				fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << std::endl;
			}
		}
		else
		{
			// Scaling factors
			for(unsigned int j=1;j<=NC_;j++)
				scaling_Jp[j] = thermodynamicsMap_.MW(j-1)/rho_;
			scaling_Jp[NC_+1] = 1./(rho_*CpMixMass_);

			const double CpT = CpMixMass_*T_;
			const double v2 = v_*v_;
			const double T3 = P_*(1.+v2/CpT)-rho_*v2;
			scaling_Jp[NC_+2] = P_*MW_/T3/T_;
			scaling_Jp[NC_+3] = 0.;
			scaling_Jp[NC_+4] = 0.;
			scaling_Jp[NC_+5] = 0.;

			
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

				// Recover additional variables
				T_   = y[NC_+1];
				rho_ = y[NC_+2];
				v_   = y[NC_+3];
				z_   = y[NC_+4];
				tau_ = y[NC_+5];

				// Calculates the volume and the concentrations of species
				thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
				P_ = rho_*PhysicalConstants::R_J_kmol*T_/MW_;
				cTot_ = P_/(PhysicalConstants::R_J_kmol * T_);
				Product(cTot_, x_, &c_);
			}

			// Calculates the current sensitivity coefficients
			sensitivityMap_->Calculate(t, T_, P_, c_, J, scaling_Jp);

			// Writes the coefficients on file (only on request)
			if (iteration_%shocktube_options_.n_step_file() == 1 || shocktube_options_.n_step_file() == 1)		
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

				for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
					fSensitivityChildXML_[indices_of_sensitivity_species_.size()+1] << sensitivityMap_->sensitivity_coefficients()(NC_+1,j-1) << " ";
				fSensitivityChildXML_[indices_of_sensitivity_species_.size()+1] << std::endl;

				for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
					fSensitivityChildXML_[indices_of_sensitivity_species_.size()+2] << sensitivityMap_->sensitivity_coefficients()(NC_+2,j-1) << " ";
				fSensitivityChildXML_[indices_of_sensitivity_species_.size()+2] << std::endl;
			}
		}
	}
}
