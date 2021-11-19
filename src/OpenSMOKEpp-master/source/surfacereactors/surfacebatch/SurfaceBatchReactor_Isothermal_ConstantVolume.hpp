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

#include "SurfaceBatchReactor_OdeInterfaces.h"
#include "math/PhysicalConstants.h"
#include "math/native-ode-solvers/MultiValueSolver"

namespace OpenSMOKE
{
	#if OPENSMOKE_USE_BZZMATH == 1
	ODESystem_BzzOde_SurfaceBatchReactor* ptOde_SurfaceBatchReactor_Isothermal_ConstantVolume;
	void ODE_Print_SurfaceBatchReactor_Isothermal_ConstantVolume(BzzVector &Y, double t)
	{
		ptOde_SurfaceBatchReactor_Isothermal_ConstantVolume->MyPrint(Y,t);
	}
	#endif

	SurfaceBatchReactor_Isothermal_ConstantVolume::SurfaceBatchReactor_Isothermal_ConstantVolume(	
								OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap, 
								OpenSMOKE::KineticsMap_Surface_CHEMKIN& kineticsSurfaceMap,
								OpenSMOKE::ODE_Parameters& ode_parameters,
								OpenSMOKE::SurfaceBatchReactor_Options& batch_options,
								const double V0, const double T0, const double P0, const double A0, 
								const OpenSMOKE::OpenSMOKEVectorDouble& omega0,
								const OpenSMOKE::OpenSMOKEVectorDouble& Z0,
								const OpenSMOKE::OpenSMOKEVectorDouble& Gamma0,
								const std::vector<bool>& site_non_conservation) :
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap), 
	thermodynamicsSurfaceMap_(thermodynamicsSurfaceMap), 
	kineticsSurfaceMap_(kineticsSurfaceMap), 
	ode_parameters_(ode_parameters),
	batch_options_(batch_options)
	{
		type_ = SURFACEBATCH_REACTOR_ISOTHERMAL_CONSTANTV;

		iteration_ = 0;
		counter_file_video_ = 0;
		counter_file_ASCII_ = 0;
		counter_file_XML_ = 0;
		counter_sensitivity_XML_ = 0;

		V0_ = V0;
		T0_ =T0 ;
		P0_ = P0;
		A0_ =  A0;
		omega0_ = omega0;
		Z0_ = Z0;
		Gamma0_ = Gamma0;
		massBulk0_ = 0.;

		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();

		// Surface phases
		SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
		SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
		SURF_NR_ = kineticsSurfaceMap_.NumberOfReactions();

		// Bulk phases
		BULK_NP_ = thermodynamicsSurfaceMap_.number_of_bulk_phases(0);
		BULK_NC_ = thermodynamicsSurfaceMap_.number_of_bulk_species();

		NE_ = NC_ + SURF_NC_ + SURF_NP_ + 1 + BULK_NC_;

		MemoryAllocation();

		MW0_ = thermodynamicsMap_.MolecularWeight_From_MassFractions(omega0_.GetHandle());
		rho_ = P0_ * MW0_ / (PhysicalConstants::R_J_kmol*T0_);
		mass0_ = rho_ * V0_;

		site_non_conservation_ = site_non_conservation;

		T_		= T0_;
		P_		= P0_;
		V_      = V0_;
		A_      = A0_;
		omega_	= omega0_;
		Z_	    = Z0_;
		Gamma_  = Gamma0_;
		GammaFromEqn_ = Gamma0_;
		MW_     = MW0_;
		QRGas_		= 0.;
		QRSurface_	= 0.;
		mass_		= mass0_;
		massBulk_ = massBulk0_;

		thermodynamicsMap_.MoleFractions_From_MassFractions(x0_.GetHandle(), MW0_, omega0_.GetHandle());
		thermodynamicsMap_.SetPressure(P0_);
		thermodynamicsMap_.SetTemperature(T0_);
		H0_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x0_.GetHandle());
		U0_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(x0_.GetHandle());

		// Initial enthalpy and internal energy for surfaces
		thermodynamicsSurfaceMap_.SetPressure(P0_);
		thermodynamicsSurfaceMap_.SetTemperature(T0_);
		H0_Surface_ = 0.;
		U0_Surface_ = 0.;
		for(unsigned int i=1;i<=SURF_NC_;i++)
		{
			const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i-1]+1;
			H0_Surface_ += Z0_[i] * thermodynamicsSurfaceMap_.Species_H_over_RT()[i+NC_-1]*PhysicalConstants::R_J_kmol*T0_ *
				           Gamma0_[index_phase] * A0_;
			U0_Surface_ += Z0_[i] * (thermodynamicsSurfaceMap_.Species_H_over_RT()[i+NC_-1]-1.)*PhysicalConstants::R_J_kmol*T0_ *
				           Gamma0_[index_phase] * A0_;
		}

		OpenAllFiles(thermodynamicsSurfaceMap_, batch_options_);
	}

	int SurfaceBatchReactor_Isothermal_ConstantVolume::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
	{
		// Recover unknowns
		{
			unsigned int k=1;
			#ifdef CHECK_MASSFRACTIONS
			for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = std::max(y[k++], 0.);
			mass_ = std::max(y[k++], 0.);
			for(unsigned int i=1;i<=SURF_NP_;++i)
				GammaFromEqn_[i] = std::max(y[k++], 0.);
			for(unsigned int i=1;i<=SURF_NC_;++i)
				Z_[i] = std::max(y[k++], 0.);
			for (unsigned int i = 1; i <= BULK_NC_; ++i)
				massBulk_[i] = y[k++];
			#else
			for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = y[k++];
			mass_ = y[k++];
			for(unsigned int i=1;i<=SURF_NP_;++i)
				GammaFromEqn_[i] = y[k++];
			for(unsigned int i=1;i<=SURF_NC_;++i)
				Z_[i] = y[k++];
			for (unsigned int i = 1; i <= BULK_NC_; ++i)
				massBulk_[i] = y[k++];
			#endif

			for (unsigned int i = 1; i <= SURF_NP_; ++i)
				if (site_non_conservation_[i - 1] == true)
					Gamma_[i] = GammaFromEqn_[i];
		}

		// Calculates the pressure and the concentrations of species
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
		rho_ = mass_/V_;
		cTot_ = rho_/MW_;
		Product(cTot_, x_, &c_);
		P_ = cTot_ * PhysicalConstants::R_J_kmol * T_;

		
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

		double dm_over_dt = A_*Dot(RfromSurface_.Size(), RfromSurface_.GetHandle(), thermodynamicsMap_.MWs().data());

		// Recovering residuals
		{
			unsigned int k=1;

			// Gas species
			for (unsigned int i=1;i<=NC_;++i)	
				dy[k++] = thermodynamicsMap_.MW(i-1)*RfromGas_[i]/rho_ +
						 (-omega_[i]*dm_over_dt + A_*RfromSurface_[i]*thermodynamicsMap_.MW(i-1))/mass_;
			
			// Gas total mass
			dy[k++] = dm_over_dt;

			// Surface phases
			for (unsigned int i=1;i<=SURF_NP_;++i)	
				dy[k++] = RsurfacePhases_[i];

			// Surface species
			for (unsigned int i=1;i<=SURF_NC_;++i)	
			{
				const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i-1]+1;
				dy[k++] = (thermodynamicsSurfaceMap_.vector_occupancies_site_species()[i-1]*Rsurface_[i] - 
										Z_[i]*dy[NC_+1+index_phase])/Gamma_[index_phase];
			}

			// Bulk species
			for (unsigned int i = 1; i <= BULK_NC_; ++i)
			{
				const unsigned int index_bulk = NC_ + SURF_NC_ + i;
				dy[k++] = Rbulk_[i] * thermodynamicsSurfaceMap_.MW(index_bulk-1) * A_;
			}
		}
		
		return 0;
	}

	int SurfaceBatchReactor_Isothermal_ConstantVolume::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		iteration_++;

		if (batch_options_.verbose_output() == true)
		{
			// Video output
			if (iteration_%batch_options_.n_step_video() == 1 || batch_options_.n_step_video() == 1)
			{
				counter_file_video_++;
				if (counter_file_video_%100 == 1)
				{
					std::cout << std::endl;
					std::cout << std::setw(6)  << std::left << "#Step";
					std::cout << std::setw(16) << std::left << "Time[s]";
					std::cout << std::setw(10) << std::left << "T[K]";
					std::cout << std::setw(10) << std::left << "P[atm]";
					if (BULK_NC_ != 0)
						std::cout << std::setw(16) << std::left << "bulk[kg]";
					std::cout << std::endl;
				}
				std::cout << std::fixed << std::setw(6) << std::left << iteration_;
				std::cout << std::scientific << std::setw(16) << std::setprecision(6) << std::left << t;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << T_;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << P_/101325.;
				if (BULK_NC_ != 0)
					std::cout << std::scientific << std::setw(16) << std::setprecision(6) << std::left << massBulk_.SumElements();
				std::cout << std::endl;
			}

			// ASCII file output
			if (batch_options_.verbose_ascii_file() == true)
			{
				if (iteration_%batch_options_.n_step_file() == 1 || batch_options_.n_step_file() == 1)
				{
					counter_file_ASCII_++;
					fASCII_.setf(std::ios::scientific);
					fASCII_ << std::setw(20) << std::left << t;
					fASCII_ << std::setw(20) << std::left << T_;
					fASCII_ << std::setw(20) << std::left << P_;
					fASCII_ << std::setw(20) << std::left << V_;
					fASCII_ << std::setw(20) << std::left << rho_;
					fASCII_ << std::setw(20) << std::left << MW_;
					fASCII_ << std::setw(20) << std::left << QRGas_;
					fASCII_ << std::setw(20) << std::left << QRSurface_;

					if (BULK_NC_ != 0)
					{
						fASCII_ << std::setw(20) << std::left << massBulk_.SumElements();
						for (unsigned int i = 1; i <= BULK_NC_; i++)
							fASCII_ << std::setw(20) << std::left << massBulk_[i];
					}

					for(unsigned int i=1;i<=SURF_NP_;i++)
						fASCII_ << std::setw(20) << std::left << Gamma_[i];

					for(unsigned int i=1;i<=SURF_NC_;i++)
						fASCII_ << std::setw(widths_of_output_surface_species_[i-1]) << std::left << Z_[i];
					
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

					if (BULK_NC_ != 0)
					{
						fBulkASCII_.setf(std::ios::scientific);
						fBulkASCII_ << std::setw(20) << std::left << t;
						fBulkASCII_ << std::setw(20) << std::left << T_;
						fBulkASCII_ << std::setw(20) << std::left << P_;
						fBulkASCII_ << std::setw(20) << std::left << V_;
						fBulkASCII_ << std::setw(20) << std::left << rho_;
						fBulkASCII_ << std::setw(20) << std::left << MW_;
						fBulkASCII_ << std::setw(20) << std::left << QRGas_;
						fBulkASCII_ << std::setw(20) << std::left << QRSurface_;

						OpenSMOKE::OpenSMOKEVectorDouble omegaBulk(BULK_NC_);	// [kg/m2/s]
						for (unsigned int i = 1; i <= BULK_NC_; i++)
						{
							unsigned int index = NC_ + SURF_NC_ + i;
							omegaBulk[i] = Rbulk_[i] * thermodynamicsSurfaceMap_.MW(index-1);
						}

						fBulkASCII_ << std::setw(24) << std::left << massBulk_.SumElements();		// [kg]
						fBulkASCII_ << std::setw(24) << std::left << omegaBulk.SumElements();		// [kg/m2/s]
						fBulkASCII_ << std::setw(24) << std::left << omegaBulk.SumElements()*A_;	// [kg/s]

						for (unsigned int i = 1; i <= BULK_NC_; i++)
						{
							const double rho = thermodynamicsSurfaceMap_.vector_densities_bulk_species()[i - 1];	// [kg/m3]
							fBulkASCII_ << std::setw(24) << std::left << massBulk_[i];								// [kg]
							fBulkASCII_ << std::setw(24) << std::left << omegaBulk[i];								// [kg/m2/s]
							fBulkASCII_ << std::setw(24) << std::left << omegaBulk[i] * A_;							// [kg/s]
							fBulkASCII_ << std::setw(24) << std::left << massBulk_[i] / rho / A_;					// [m]
							fBulkASCII_ << std::setw(24) << std::left << massBulk_[i] / rho / A_ * 1e6;				// [micron]
							fBulkASCII_ << std::setw(24) << std::left << omegaBulk[i] / rho;						// [m/s]
							fBulkASCII_ << std::setw(24) << std::left << omegaBulk[i] / rho * 1e3 * 3600;			// [mm/hr]
						}

						fBulkASCII_ << std::endl;
					}
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
					fXML_ << QRGas_ << " ";
					fXML_ << QRSurface_ << " ";
					for (unsigned int i=1;i<=NC_;i++)
						fXML_ << std::setprecision(12) << omega_[i] << " ";
					fXML_ << std::endl;
				}
			}
		}

		if (batch_options_.sensitivity_analysis() == true)
			SensitivityAnalysis(t, y);

		return 0;
	}

	void SurfaceBatchReactor_Isothermal_ConstantVolume::Solve(const double tf)
	{
		if (batch_options_.verbose_output() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << " Solving the surface batch reactor...                                      " << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
		}

		unsigned int k = 1;
		for(unsigned int i=1;i<=NC_;++i)
			y0_[k++] = omega0_[i];
		y0_[k++] = mass0_;
		for(unsigned int i=1;i<=SURF_NP_;++i)
			y0_[k++] = Gamma0_[i];
		for(unsigned int i=1;i<=SURF_NC_;++i)
			y0_[k++] = Z0_[i];
		for (unsigned int i = 1; i <= BULK_NC_; ++i)
			y0_[k++] = massBulk0_;
		
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
			typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE_SurfaceBatchReactor> denseOde;
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

			ODESystem_BzzOde_SurfaceBatchReactor odebatch(*this);
			BzzOdeStiffObject o(y0_bzz, 0., &odebatch);
			
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			o.SetTolAbs(ode_parameters_.absolute_tolerance());
			o.SetTolRel(ode_parameters_.relative_tolerance());

			ptOde_SurfaceBatchReactor_Isothermal_ConstantVolume = &odebatch;
			o.StepPrint(ODE_Print_SurfaceBatchReactor_Isothermal_ConstantVolume);

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
			SolveOpenSourceSolvers(tf, ode_parameters_);
		}

		if (batch_options_.verbose_output() == true)
		{
			ode_parameters_.Status(std::cout);
			FinalStatus(tf, thermodynamicsMap_, thermodynamicsSurfaceMap_);
			FinalSummary(batch_options_.output_path() / "FinalSummary.out", tf, thermodynamicsMap_, thermodynamicsSurfaceMap_);
		}

		CloseAllFiles(batch_options_);
	}

	void SurfaceBatchReactor_Isothermal_ConstantVolume::SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		FatalErrorMessage("SurfaceBatchReactor_Isothermal_ConstantVolume::SensitivityAnalysis not yet implemented");

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
				scaling_Jp[j] = thermodynamicsMap_.MW(j-1)/rho_;

			// Calculates the current Jacobian
			NumericalJacobian(t, y, J);

			// Recover variables
			{
				#ifdef CHECK_MASSFRACTIONS
				for(unsigned int i=1;i<=NC_;++i)
					omega_[i] = std::max(y[i], 0.);
				#else
				for(unsigned int i=1;i<=NC_;++i)
					omega_ = y;
				#endif

				// Calculates the pressure and the concentrations of species
				thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
				cTot_ = rho_/MW_;
				Product(cTot_, x_, &c_);
				P_ = cTot_ * PhysicalConstants::R_J_kmol * T_;
			}

			// Calculates the current sensitivity coefficients
			sensitivityMap_->Calculate(t, T_, P_, c_, J, scaling_Jp);

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
}
