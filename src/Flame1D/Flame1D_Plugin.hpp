/*-----------------------------------------------------------------------*\
|     ____            _  ______ __  __  ____  _  ________                 |
|    / __ \       _  (_)/  ___ |  \/  |/ __ \| |/ /  ____|                |
|   | |  | |_ __ | |_ _|  (___ | \  / | |  | | ' /| |__    _     _        |
|   | |  | | '_ \|  _| |\___  \| |\/| | |  | |  < |  __| _| |_ _| |_      |
|   | |__| | |_) | |_| |____)  | |  | | |__| | . \| |___|_   _|_   _|     |
|    \____/| .__/\___|_|______/|_|  |_|\____/|_|\_\______||_|   |_|       |
|          | |                                                            |
|          |_|                                                            |
|                                                                         |
|            Authors: Magnus Fürst <magnus.furst@ulb.ac.be>               |
|                     Andrea Bertolino <andrea.bertolino@ulb.be>          |
|-------------------------------------------------------------------------|
|   License                                                               |
|                                                                         |
|   This file is part of OptiSMOKE.                                       |
|   Copyright (C) 2020 by Magnus Fürst and Andrea Bertolino               |
|                                                                         |
|   OptiSMOKE is free software: you can redistribute it and/or modify     |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OptiSMOKE is distributed in the hope that it will be useful,          |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OptiSMOKE. If not, see <http://www.gnu.org/licenses/>.     |
|                                                                         |
\*-----------------------------------------------------------------------*/
#include <cmath>
namespace OpenSMOKE
{
	std::vector<double> Flame1D_Plugin::Setup_and_Solve(const std::string input_file_name,
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML,
		OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML,
		OpenSMOKE::TransportPropertiesMap_CHEMKIN* 	transportMapXML,
		std::string Species_of_interest,
		std::vector<double> Abscissa_x_m)
	{
		// Pointers
		thermodynamicsMapXML_ = thermodynamicsMapXML;
		kineticsMapXML_ = kineticsMapXML;
		transportMapXML_ = transportMapXML;

		// Defines the grammar rules
		OpenSMOKE::Grammar_LaminarFlame grammar_laminarflame;

		// Define the dictionaries
		std::string main_dictionary_name_ = "PremixedLaminarFlame1D";
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(main_dictionary_name_).SetGrammar(grammar_laminarflame);

		// Input parameters
		std::vector<double> inlet_T;
		std::vector<double> P_Pa;
		double inlet_velocity = 0.;
		double inlet_mass_flux = 0.;
		std::vector<OpenSMOKE::OpenSMOKEVectorDouble> inlet_omega;
		std::vector<double> equivalence_ratios;

		double outlet_T;
		OpenSMOKE::OpenSMOKEVectorDouble outlet_omega;

		// Read inlet conditions
		{
		std::vector<std::string> list_of_strings;
		if (dictionaries(main_dictionary_name_).CheckOption("@InletStream") == true)
			dictionaries(main_dictionary_name_).ReadOption("@InletStream", list_of_strings);

		// If multiple inlet streams are specified
		if (list_of_strings.size() != 1)
		{
			inlet_T.resize(list_of_strings.size());
			P_Pa.resize(list_of_strings.size());
			inlet_omega.resize(list_of_strings.size());
			equivalence_ratios.resize(list_of_strings.size());
			for (unsigned int i = 0; i < list_of_strings.size(); i++)
				GetGasStatusFromDictionary(dictionaries(list_of_strings[i]), *thermodynamicsMapXML_, inlet_T[i], P_Pa[i], inlet_omega[i]);
		}
		// If a single inlet stream is defined
		else
		{
			GetGasStatusFromDictionary(dictionaries(list_of_strings[0]), *thermodynamicsMapXML_, inlet_T, P_Pa, inlet_omega, equivalence_ratios);
		}
	}
	
	// Read outlet conditions
	{
		double P_Pa_outlet;
		std::string name_of_gas_status_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@OutletStream") == true)
			dictionaries(main_dictionary_name_).ReadDictionary("@OutletStream", name_of_gas_status_subdictionary);

		GetGasStatusFromDictionary(dictionaries(name_of_gas_status_subdictionary), *thermodynamicsMapXML_, outlet_T, P_Pa_outlet, outlet_omega);

		if (P_Pa_outlet != P_Pa[0])
		{
			OpenSMOKE::FatalErrorMessage("The pressure of outlet stream does not match with the inlet stream");
		}
	}

	// Read inlet velocity
	{
		double value;
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@InletVelocity") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@InletVelocity", value, units);
			if (units == "m/s")			inlet_velocity = value;
			else if (units == "cm/s")   inlet_velocity = value/100.;
			else if (units == "mm/s")	inlet_velocity = value/1000.;
			else if (units == "km/h")   inlet_velocity = value*10./36.;
			else 
			{
				OpenSMOKE::FatalErrorMessage("Unknown velocity units");
			}
		}
	}

	Eigen::VectorXd w;
	
	// Adaptive grid
	//OpenSMOKE::Grid1D* grid;
	{
		std::string name_of_adaptive_grid_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@Grid") == true)
			dictionaries(main_dictionary_name_).ReadDictionary("@Grid", name_of_adaptive_grid_subdictionary);

		grid = new OpenSMOKE::Grid1D(dictionaries(name_of_adaptive_grid_subdictionary), w);
	}
	
	{
		std::shared_ptr<OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D> ciao_dario = std::make_shared<OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D>(*thermodynamicsMapXML, *kineticsMapXML, *transportMapXML, *grid);
		
		flame_premixed = ciao_dario;
	}
	// Output folder
	if (dictionaries(main_dictionary_name_).CheckOption("@Output") == true)
	{
		boost::filesystem::path output_folder;
		dictionaries(main_dictionary_name_).ReadPath("@Output", output_folder);
		flame_premixed->SetOutputFolder(output_folder);
	}

	// Solver type
	{
		std::string solver_type;
		if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true)
		{
			dictionaries(main_dictionary_name_).ReadString("@Type", solver_type);
			flame_premixed->SetSolverType(solver_type);
		}
	}

	// Soret effect
	{
		bool soret = true;
		if (dictionaries(main_dictionary_name_).CheckOption("@Soret") == true)
			dictionaries(main_dictionary_name_).ReadBool("@Soret",soret);
		flame_premixed->SetSoret(soret);
	}

	// Radiative heat transfer
	{
		bool radiative_heat_transfer = false;
		if (dictionaries(main_dictionary_name_).CheckOption("@Radiation") == true)
			dictionaries(main_dictionary_name_).ReadBool("@Radiation", radiative_heat_transfer);
		flame_premixed->SetRadiativeHeatTransfer(radiative_heat_transfer);
	}

	// Read environment temperature
	{
		double value;
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@EnvironmentTemperature") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@EnvironmentTemperature", value, units);
			if (units == "K")			value *= 1.;
			else if (units == "C")		value += 273.15;
			else 
			{
				OpenSMOKE::FatalErrorMessage("Unknown temperature units");
			}

			flame_premixed->SetEnvironmentTemperature(value);
		}
	}

	// Polimi soot
	polimi_soot_used = false;
	//OpenSMOKE::PolimiSoot_Analyzer* polimi_soot;
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@PolimiSoot") == true)
		{
			std::string name_of_polimisoot_analyzer_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@PolimiSoot", name_of_polimisoot_analyzer_subdictionary);
			polimi_soot_used = true;
			polimi_soot = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML, dictionaries(name_of_polimisoot_analyzer_subdictionary));
			
			if (polimi_soot->number_sections() != 0)
				flame_premixed->SetPolimiSoot(polimi_soot);
		}
	}

	// On the fly PostProcessing
	//OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing;
	{
		on_the_fly_post_processing = new OpenSMOKE::OnTheFlyPostProcessing(*thermodynamicsMapXML, *kineticsMapXML, flame_premixed->output_folder());

		if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyPostProcessing") == true)
		{
			std::string name_of_options_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyPostProcessing", name_of_options_subdictionary);
			on_the_fly_post_processing->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
			flame_premixed->SetOnTheFlyPostProcessing(on_the_fly_post_processing);
		}
	}

	// Fixed temperature profile
	{
		std::string name_of_fixed_temperature_profile_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@FixedTemperatureProfile") == true)
		{
			dictionaries(main_dictionary_name_).ReadDictionary("@FixedTemperatureProfile", name_of_fixed_temperature_profile_subdictionary);

			OpenSMOKE::OpenSMOKEVectorDouble x, y;
			std::string x_variable, y_variable;
			GetXYProfileFromDictionary(dictionaries(name_of_fixed_temperature_profile_subdictionary), x, y, x_variable, y_variable);

			if (x_variable != "length")
			{
				OpenSMOKE::FatalErrorMessage("The @FixedTemperatureProfile must be defined versus spacee");
			}
			if (y_variable != "temperature")
			{
				OpenSMOKE::FatalErrorMessage("The @FixedTemperatureProfile must be define the temperature profile");
			}

			flame_premixed->SetFixedTemperatureProfile(x, y);
		}
	}

	// Fixed specific (i.e. per unit area) mass flow rate profile
	{
		std::string name_of_fixed_specific_mass_flow_rate_profile_subdictionary;
		if (dictionaries(main_dictionary_name_).CheckOption("@FixedSpecificMassFlowRateProfile") == true)
		{
			dictionaries(main_dictionary_name_).ReadDictionary("@FixedSpecificMassFlowRateProfile", name_of_fixed_specific_mass_flow_rate_profile_subdictionary);

			OpenSMOKE::OpenSMOKEVectorDouble x, y;
			std::string x_variable, y_variable;
			GetXYProfileFromDictionary(dictionaries(name_of_fixed_specific_mass_flow_rate_profile_subdictionary), x, y, x_variable, y_variable);

			if (x_variable != "length")
			{
				OpenSMOKE::FatalErrorMessage("The @FixedSpecificMassFlowRateProfile must be defined versus spacee");
			}
			if (y_variable != "specific-mass-flow-rate")
			{
				OpenSMOKE::FatalErrorMessage("The @FixedSpecificMassFlowRateProfile must be define the specific (i.e. per unit area) mass flow rate profile");
			}

			flame_premixed->SetFixedSpecificMassFlowRateProfile(x, y);
		}
	}

	// Read inlet velocity
	{
		double value;
		std::string units;
		if (dictionaries(main_dictionary_name_).CheckOption("@FixedOutletTemperature") == true)
		{
			dictionaries(main_dictionary_name_).ReadMeasure("@FixedOutletTemperature", value, units);
			if (units == "K")		value = value;
			else if (units == "C")  value = value+273.15;
			else 
			{
				OpenSMOKE::FatalErrorMessage("Unknown temperature units");
			}

			flame_premixed->SetFixedOutletTemperature(value);
		}
	}

	// Wall heat exchange
	if ( dictionaries(main_dictionary_name_).CheckOption("@WallHeatExchangeCoefficient") == true ||
		 dictionaries(main_dictionary_name_).CheckOption("@WallHeatNusseltNumber") == true )
	{
		
		double wall_heat_exchange_coefficient = 0.;
		double wall_heat_nusselt_number = 0.;

		if (dictionaries(main_dictionary_name_).CheckOption("@WallHeatExchangeCoefficient") == true)
		{
			std::string units_wall_heat_exchange_coefficient;
			dictionaries(main_dictionary_name_).ReadMeasure("@WallHeatExchangeCoefficient", wall_heat_exchange_coefficient, units_wall_heat_exchange_coefficient);
			if (units_wall_heat_exchange_coefficient == "W/m2/K")			wall_heat_exchange_coefficient = wall_heat_exchange_coefficient;
			else if (units_wall_heat_exchange_coefficient == "W/m2/C")		wall_heat_exchange_coefficient = wall_heat_exchange_coefficient;
			else if (units_wall_heat_exchange_coefficient == "kcal/m2/K")	wall_heat_exchange_coefficient = wall_heat_exchange_coefficient * 4186.8;
			else if (units_wall_heat_exchange_coefficient == "kcal/m2/C")	wall_heat_exchange_coefficient = wall_heat_exchange_coefficient * 4186.8;

			else 
			{
				OpenSMOKE::FatalErrorMessage("Unknown heat exchange coefficient");
			}
		}
		else if (dictionaries(main_dictionary_name_).CheckOption("@WallHeatNusseltNumber") == true)
		{
			dictionaries(main_dictionary_name_).ReadDouble("@WallHeatNusseltNumber", wall_heat_nusselt_number);
		}

		OpenSMOKE::OpenSMOKEVectorDouble wall_temperature_profile_x;
		OpenSMOKE::OpenSMOKEVectorDouble wall_temperature_profile_y;
		{
			std::string name_of_wall_temperature_profile_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@WallTemperatureProfile", name_of_wall_temperature_profile_subdictionary);
			std::string x_variable, y_variable;
			GetXYProfileFromDictionary(dictionaries(name_of_wall_temperature_profile_subdictionary), wall_temperature_profile_x, wall_temperature_profile_y, x_variable, y_variable);

			if (x_variable != "length")
			{
				OpenSMOKE::FatalErrorMessage("The @WallTemperatureProfile must be defined versus space");
			}
			if (y_variable != "temperature")
			{
				OpenSMOKE::FatalErrorMessage("The @WallTemperatureProfile must be define the temperature profile");
			}
		}
		
		double internal_diameter = 0.;
		{
			std::string units_internal_diameter;
			dictionaries(main_dictionary_name_).ReadMeasure("@InternalDiameter", internal_diameter, units_internal_diameter);
			if (units_internal_diameter == "m")			internal_diameter = internal_diameter;
			else if (units_internal_diameter == "cm")	internal_diameter = internal_diameter * 0.01;
			else if (units_internal_diameter == "mm")	internal_diameter = internal_diameter * 0.001;
			else 
			{
				OpenSMOKE::FatalErrorMessage("Unknown length units");
			}
		}

		flame_premixed->SetWallHeatExchange(	wall_heat_exchange_coefficient, wall_heat_nusselt_number, internal_diameter,
												wall_temperature_profile_x, wall_temperature_profile_y);
	}
	
	// Lewis numbers
	if (dictionaries(main_dictionary_name_).CheckOption("@LewisNumbers") == true)
	{
		std::string name_of_lewis_numbers_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@LewisNumbers", name_of_lewis_numbers_subdictionary);

		std::vector<double> lewis_numbers;
		OpenSMOKE::GetLewisNumbersFromDictionary(dictionaries(name_of_lewis_numbers_subdictionary), thermodynamicsMapXML[0], lewis_numbers);
		flame_premixed->SetLewisNumbers(lewis_numbers);
	}

	// Initialize from backup
	std::cout.setstate(std::ios_base::failbit); // Disable video output
	bool use_userdefined_grid_for_backup = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@Backup") == true)
	{
		if (dictionaries(main_dictionary_name_).CheckOption("@DontUseBackupGrid") == true)
			dictionaries(main_dictionary_name_).ReadBool("@DontUseBackupGrid", use_userdefined_grid_for_backup);

		if (flame_premixed->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_FLAMESPEED)
		{
			boost::filesystem::path path_backup;
			dictionaries(main_dictionary_name_).ReadPath("@Backup", path_backup);

			// Set inlet and outlet values and first guess velocity according to user values
			flame_premixed->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
			flame_premixed->SetOutlet(outlet_T, outlet_omega);
			if (inlet_velocity>0.) flame_premixed->SetInletVelocity(inlet_velocity);
			else                   flame_premixed->SetInletMassFlux(inlet_mass_flux);

			// Setup the solution, accordingly to backup file
			flame_premixed->InitializeFromBackupFile(path_backup, use_userdefined_grid_for_backup);
		}
		else if (flame_premixed->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_BURNERSTABILIZED)
		{
			boost::filesystem::path path_backup;
			dictionaries(main_dictionary_name_).ReadPath("@Backup", path_backup);

			// Set inlet and outlet values and first guess velocity according to user values
			flame_premixed->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
			flame_premixed->SetOutlet(outlet_T, outlet_omega);
			if (inlet_velocity > 0.) flame_premixed->SetInletVelocity(inlet_velocity);
			else                     flame_premixed->SetInletMassFlux(inlet_mass_flux);

			// Setup the solution, accordingly to backup file
			flame_premixed->InitializeFromBackupFile(path_backup, use_userdefined_grid_for_backup);
		}
	}
	else
	{
		if (flame_premixed->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_FLAMESPEED)
		{
			flame_premixed->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
			flame_premixed->SetOutlet(outlet_T, outlet_omega);
			if (inlet_velocity > 0.) flame_premixed->SetInletVelocity(inlet_velocity);
			else                     flame_premixed->SetInletMassFlux(inlet_mass_flux);
			flame_premixed->SetupForFlameSpeed(w);
		}
		else if (flame_premixed->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_BURNERSTABILIZED)
		{
			flame_premixed->SetInlet(inlet_T[0], P_Pa[0], inlet_omega[0]);
			flame_premixed->SetOutlet(outlet_T, outlet_omega);
			if (inlet_velocity > 0.) flame_premixed->SetInletVelocity(inlet_velocity);
			else                     flame_premixed->SetInletMassFlux(inlet_mass_flux);
			flame_premixed->SetupForBurnerStabilized(w);
		}
	}
	std::cout.clear(); // Re-enable video output

	// Use NLS Solver
	{
		bool use_nls_solver = true;
		if (dictionaries(main_dictionary_name_).CheckOption("@UseNlsSolver") == true)
			dictionaries(main_dictionary_name_).ReadBool("@UseNlsSolver", use_nls_solver);
		flame_premixed->SetUseNlsSolver(use_nls_solver);
	}

	// Use DAE Solver
	{
		bool use_dae_solver = true;
		if (dictionaries(main_dictionary_name_).CheckOption("@UseDaeSolver") == true)
			dictionaries(main_dictionary_name_).ReadBool("@UseDaeSolver", use_dae_solver);
		flame_premixed->SetUseDaeSolver(use_dae_solver);
	}

	// Sensitivity Options
	sensitivity_options_used = false;
	//OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options;
	if (dictionaries(main_dictionary_name_).CheckOption("@SensitivityAnalysis") == true)
	{
		sensitivity_options_used = true;
		sensitivity_options = new OpenSMOKE::SensitivityAnalysis_Options();
		std::string name_of_sensitivity_options_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@SensitivityAnalysis", name_of_sensitivity_options_subdictionary);
		sensitivity_options->SetupFromDictionary(dictionaries(name_of_sensitivity_options_subdictionary));
		flame_premixed->EnableSensitivityAnalysis(*sensitivity_options);
	}

	// Dae Options
	//DaeSMOKE::DaeSolver_Parameters* dae_parameters;
	dae_parameters = new DaeSMOKE::DaeSolver_Parameters();
	if (dictionaries(main_dictionary_name_).CheckOption("@DaeParameters") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@DaeParameters", name_of_subdictionary);
		dae_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
	}

	// Nls Options
	//NlsSMOKE::NonLinearSolver_Parameters* nls_parameters;
	nls_parameters = new NlsSMOKE::NonLinearSolver_Parameters();
	if (dictionaries(main_dictionary_name_).CheckOption("@NlsParameters") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@NlsParameters", name_of_subdictionary);
		nls_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
	}
	
	// Pseudo Transient Options
	//NlsSMOKE::FalseTransientSolver_Parameters* false_transient_parameters;
	false_transient_parameters = new NlsSMOKE::FalseTransientSolver_Parameters();
	if (dictionaries(main_dictionary_name_).CheckOption("@FalseTransientParameters") == true)
	{
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@FalseTransientParameters", name_of_subdictionary);
		false_transient_parameters->SetupFromDictionary(dictionaries(name_of_subdictionary));
	}

	if (dictionaries(main_dictionary_name_).CheckOption("@DerivativeGasTemperature") == true)
	{
		std::string value;
		dictionaries(main_dictionary_name_).ReadString("@DerivativeGasTemperature", value);
		if (value == "upwind")					flame_premixed->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_UPWIND);
		else if (value == "backward")			flame_premixed->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_BACKWARD);
		else if (value == "forward")			flame_premixed->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_FORWARD);
		else if (value == "centered")			flame_premixed->SetDerivativeGasTemperature(OpenSMOKE::DERIVATIVE_1ST_CENTERED);
		else 
		{
			OpenSMOKE::FatalErrorMessage("Unknown derivative type for gas temperature");
		}
	}

	if (dictionaries(main_dictionary_name_).CheckOption("@DerivativeGasMassFractions") == true)
	{
		std::string value;
		dictionaries(main_dictionary_name_).ReadString("@DerivativeGasMassFractions", value);
		if (value == "upwind")					flame_premixed->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_UPWIND);
		else if (value == "backward")			flame_premixed->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_BACKWARD);
		else if (value == "forward")			flame_premixed->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_FORWARD);
		else if (value == "centered")			flame_premixed->SetDerivativeGasMassFractions(OpenSMOKE::DERIVATIVE_1ST_CENTERED);
		else 
		{
			OpenSMOKE::FatalErrorMessage("Unknown derivative type for gas mass fractions");
		}
	}

	// Hybrid Method of Moments
	hmom_used = false;
	//OpenSMOKE::HMOM* hmom;
	if (dictionaries(main_dictionary_name_).CheckOption("@HMOM") == true)
	{
		hmom_used = true;
		hmom = new OpenSMOKE::HMOM();
		std::string name_of_subdictionary;
		dictionaries(main_dictionary_name_).ReadDictionary("@HMOM", name_of_subdictionary);
		hmom->SetupFromDictionary(dictionaries(name_of_subdictionary));

		flame_premixed->SolveHMOMFromExistingSolution(*hmom, *dae_parameters, *nls_parameters, *false_transient_parameters);
		std::vector<double> empty_vector;
		return(empty_vector);
	}
	
	if (flame_premixed->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_FLAMESPEED)
	{       
		//std::cout << "\n sono dentro quello stupido solver \n";
		// Solve only for a single flame
		if (inlet_omega.size() == 1)
		{	
			time_t timerStart;
			time_t timerEnd;
			
			time(&timerStart);
			
			std::cout.setstate(std::ios_base::failbit); // Disable video output
			
			flame_premixed->SolveFlameSpeedFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);
			//std::cout << "Ho finito figlio di puttana \n";	
			time(&timerEnd);

			std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
			std::cout.clear(); // Re-enable video output
		} else
		{
			OpenSMOKE::FatalErrorMessage("Only one flame can be solved at time!");		
		}
		//double flame_speed_return = 7.0;
		std::vector<double> flame_speed_return(2,0); 
		flame_speed_return[0] = flame_premixed->flame_speed()*100.;
		clean_up();
		//std::cout << "Debug \n ";
        	//std::cout << flame_speed_return[0]<< std::endl;
		return flame_speed_return;
	}


	if (flame_premixed->solver_type() == OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D::SOLVER_TYPE_BURNERSTABILIZED)
	{
		// Solve only for a single flame
		if (inlet_omega.size() == 1)
		{
			time_t timerStart;
			time_t timerEnd;
			
			time(&timerStart);
			std::cout.setstate(std::ios_base::failbit); // Disable video output
			// solve the burner stabilized flame from backup!
			flame_premixed->SolveBurnerStabilizedFromScratch(*dae_parameters, *nls_parameters, *false_transient_parameters);
			time(&timerEnd);

			std::cout << "Total time: " << difftime(timerEnd, timerStart) << " s" << std::endl;
			std::cout.clear(); // Re-enable video output
		} else
		{
			OpenSMOKE::FatalErrorMessage("Only one flame can be solved at a time!");		
		}


	// Get the coordinates list from the flame object
	Eigen::VectorXd coordinates_list = flame_premixed->grid().x();


	// Get the mole fractions of all species
	std::vector<Eigen::VectorXd> solution_tot = flame_premixed->X();
	// Get the species and the water index from the thermo map
	int species_index = thermodynamicsMapXML_->IndexOfSpecies(Species_of_interest)-1;
	int species_water_index = thermodynamicsMapXML_->IndexOfSpecies("H2O")-1;


	// Make the transformation to dry composition 
	std::vector<double> dry_species_profile(solution_tot.size());
	for (int i=0; i< solution_tot.size(); i++){
		dry_species_profile[i] =std::pow(10,6)*solution_tot[i](species_index)/(1-solution_tot[i](species_water_index));
	}

	// initialize vector for the solution and fill it up
	std::vector<double> solution_vector;
	for (int i=0; i< Abscissa_x_m.size(); i++)
	{
		double temp_solution = Interpolate(coordinates_list, dry_species_profile, Abscissa_x_m[i]);
		solution_vector.push_back(temp_solution);
	}

	clean_up();
	std::cout << "Debug \n ";
	std::cout << solution_vector[0]<< std::endl;
	return solution_vector;

	}

}

double Flame1D_Plugin::Interpolate(Eigen::VectorXd grid_coordinates, std::vector<double> species_profile, double Abscissa_temp)
{
	double interpolated_point;

	if (Abscissa_temp > grid_coordinates[grid_coordinates.size()-1])
	{
		std::cout<<"WARNING_OptiSMOKE: The time (" + std::to_string(Abscissa_temp) + ") that you specified into the EXP_DATA is outside the runtime of this simulation (" + std::to_string(grid_coordinates[grid_coordinates.size()-1]) + ")! Extrapolation was not possible! Please be sure that this is what you want."<<std::endl;
		interpolated_point=10000000;

	} else {

		for (int j=0; j<grid_coordinates.size(); j++)
		{
			if(Abscissa_temp <= grid_coordinates(j))
			{
				interpolated_point = species_profile[j-1] + (species_profile[j]-species_profile[j-1])/(grid_coordinates(j)-grid_coordinates(j-1))*(Abscissa_temp-grid_coordinates(j-1));
				//std::cout<<"Index of the time vector where the condition is satisfied is: "<< j << std::endl;
				//std::cout<<"The molar fraction of the species "<< Species_name <<" is " << interpolated_point <<" at abscissa "<<Abscissa_temp<< std::endl;
				break;
			}
		}
		
	}

	return interpolated_point;
}

void Flame1D_Plugin::clean_up()
{
	delete grid;
	
	//delete flame_premixed;
	
	if (sensitivity_options_used)
	{
		delete sensitivity_options;
	}
	delete on_the_fly_post_processing;
	if (polimi_soot_used)
	{
		delete polimi_soot;
	}
	delete dae_parameters;
	delete nls_parameters;
	delete false_transient_parameters;
	if (hmom_used)
	{
		delete hmom;
	}
	//std::cout<< " Address flame_premixed " << flame_premixed <<std::endl;
	//delete [] flame_premixed;
	//std::cout<< " Address flame_premixed " << flame_premixed <<std::endl;
}

}
