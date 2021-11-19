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

namespace OpenSMOKE
{
	void BatchReactor_Plugin::Setup(const std::string input_file_name,
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML,
		OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML)
	{
		// Pointers
		thermodynamicsMapXML_ = thermodynamicsMapXML;
		kineticsMapXML_ = kineticsMapXML;

		// Defines the grammar rules
		OpenSMOKE::Grammar_BatchReactor grammar_BatchReactor;

		// Define the dictionaries
		std::string main_dictionary_name_ = "BatchReactor";
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(main_dictionary_name_).SetGrammar(grammar_BatchReactor);

		tEnd_ = 0.;
		tStart_ = 0.;                       // default 0
		volume = 1.;					// default value [1 m3]

		// Read initial conditions
		{
			std::string name_of_gas_status_subdictionary;
			if (dictionaries(main_dictionary_name_).CheckOption("@InitialStatus") == true)
				dictionaries(main_dictionary_name_).ReadDictionary("@InitialStatus", name_of_gas_status_subdictionary);

			GetGasStatusFromDictionary(dictionaries(name_of_gas_status_subdictionary), *thermodynamicsMapXML_, T, P_Pa, omega);
		}

		// Read end time
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@EndTime") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@EndTime", value, units);
				if (units == "s")		  tEnd_ = value;
				else if (units == "ms")   tEnd_ = value / 1000.;
				else if (units == "min")  tEnd_ = value * 60.;
				else if (units == "h")    tEnd_ = value * 3600.;
				else OpenSMOKE::FatalErrorMessage("Unknown time units");
			}
		}

		// Read start time
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@StartTime") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@StartTime", value, units);
				if (units == "s")         tStart_ = value;
				else if (units == "ms")   tStart_ = value / 1000.;
				else if (units == "min")  tStart_ = value * 60.;
				else if (units == "h")    tStart_ = value * 3600.;
				else OpenSMOKE::FatalErrorMessage("Unknown time units");
			}
		}

		// Read volume
		{
			if (dictionaries(main_dictionary_name_).CheckOption("@Volume") == true)
			{
				double value;
				std::string units;

				dictionaries(main_dictionary_name_).ReadMeasure("@Volume", value, units);
				if (units == "m3")		  volume = value;
				else if (units == "dm3")  volume = value / 1.e3;
				else if (units == "cm3")  volume = value / 1.e6;
				else if (units == "mm3")  volume = value / 1.e9;
				else if (units == "l")    volume = value / 1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown volume units");
			}
		}

		// Read exchange area
		exchange_area = 0.;
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@ExchangeArea") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@ExchangeArea", value, units);
				if (units == "m2")        exchange_area = value;
				else if (units == "dm2")  exchange_area = value / 1.e2;
				else if (units == "cm2")  exchange_area = value / 1.e4;
				else if (units == "mm2")  exchange_area = value / 1.e6;
				else OpenSMOKE::FatalErrorMessage("Unknown area units");
			}
		}

		// Read global thermal exchange coefficient
		global_thermal_exchange_coefficient = 0.;
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@GlobalThermalExchangeCoefficient") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@GlobalThermalExchangeCoefficient", value, units);
				if (units == "W/m2/K")			global_thermal_exchange_coefficient = value;
				else if (units == "W/m2/C")		global_thermal_exchange_coefficient = value;
				else if (units == "kcal/m2/K")		global_thermal_exchange_coefficient = value * 4186.8;
				else if (units == "kcal/m2/C")		global_thermal_exchange_coefficient = value * 4186.8;
				else OpenSMOKE::FatalErrorMessage("Unknown global thermal exchange coefficient units");
			}
		}

		// Environment temperature
		T_environment = T;
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@EnvironmentTemperature") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@EnvironmentTemperature", value, units);
				if (units == "K")		T_environment = value;
				else if (units == "C")	T_environment = value + 273.15;
				else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
			}
		}


		//Type
		{
			std::string value;
			if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true)
			{
				dictionaries(main_dictionary_name_).ReadString("@Type", value);
				if (value == "Isothermal-ConstantVolume")				type_ = OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV;
				else if (value == "Isothermal-ConstantPressure")		type_ = OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP;
				else if (value == "NonIsothermal-ConstantVolume")		type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV;
				else if (value == "NonIsothermal-ConstantPressure")		type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP;
				else if (value == "NonIsothermal-UserDefinedVolume")	type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME;
				else OpenSMOKE::FatalErrorMessage("Unknown batch reactor type: " + value);

			}
		}

		// Options
		{
			batch_options_ = new OpenSMOKE::BatchReactor_Options();
			if (dictionaries(main_dictionary_name_).CheckOption("@Options") == true)
			{
				std::string name_of_options_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@Options", name_of_options_subdictionary);
				batch_options_->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
			}
			batch_options_->SetVerboseVideo(false);
			batch_options_->SetVerboseOutput(false);
			batch_options_->SetSensitivityAnalysis(false);
		}

		// ODE Parameters
		{
			ode_parameters_ = new OpenSMOKE::ODE_Parameters();
			if (dictionaries(main_dictionary_name_).CheckOption("@OdeParameters") == true)
			{
				std::string name_of_ode_parameters_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@OdeParameters", name_of_ode_parameters_subdictionary);
				ode_parameters_->SetupFromDictionary(dictionaries(name_of_ode_parameters_subdictionary));
			}
		}

		// Sensitivity Options
		{
			sensitivity_options_ = new OpenSMOKE::SensitivityAnalysis_Options();
			if (dictionaries(main_dictionary_name_).CheckOption("@SensitivityAnalysis") == true)
			{
				//std::string name_of_sensitivity_options_subdictionary;
				//dictionaries(main_dictionary_name_).ReadDictionary("@SensitivityAnalysis", name_of_sensitivity_options_subdictionary);

				//batch_options_->SetSensitivityAnalysis(true);
				//sensitivity_options_->SetupFromDictionary(dictionaries(name_of_sensitivity_options_subdictionary));
			}
		}

		// On the fly ROPA
		onTheFlyROPA_ = new OpenSMOKE::OnTheFlyROPA(*thermodynamicsMapXML_, *kineticsMapXML_);
		if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyROPA") == true)
		{
			std::string name_of_options_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyROPA", name_of_options_subdictionary);
			// No ROPA (disabled)
			//onTheFlyROPA_->SetupFromDictionary(dictionaries(name_of_options_subdictionary), path_kinetics_output);
		}

		// On the fly CEMA
		onTheFlyCEMA_ = new OpenSMOKE::OnTheFlyCEMA(*thermodynamicsMapXML_, *kineticsMapXML_, batch_options_->output_path());
		if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyCEMA") == true)
		{
			std::string name_of_options_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyCEMA", name_of_options_subdictionary);
			onTheFlyCEMA_->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
		}

		// On the fly PostProcessing
		{
			on_the_fly_post_processing_ = new OpenSMOKE::OnTheFlyPostProcessing(*thermodynamicsMapXML_, *kineticsMapXML_, batch_options_->output_path());

			if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyPostProcessing") == true)
			{
				std::string name_of_options_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyPostProcessing", name_of_options_subdictionary);
				on_the_fly_post_processing_->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
			}
		}

		// Ignition Delay Times
		{
			idt = new OpenSMOKE::IgnitionDelayTimes_Analyzer();
			if (dictionaries(main_dictionary_name_).CheckOption("@IgnitionDelayTimes") == true)
			{
				if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV ||
					type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP)
					OpenSMOKE::FatalErrorMessage("The @IgnitionDelayTimes can be used only for NonIsothermal reactors");

				std::string name_of_idt_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@IgnitionDelayTimes", name_of_idt_subdictionary);
				idt->SetupFromDictionary(dictionaries(name_of_idt_subdictionary), *thermodynamicsMapXML_);
			}
		}

		{
		volume_profile = false; 
			// Read pressure coefficient
			if (dictionaries(main_dictionary_name_).CheckOption("@PressureCoefficient") == true)
			{
				batchreactor_volumeprofile = new OpenSMOKE::BatchReactor_VolumeProfile();

				if (type_ != OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
					OpenSMOKE::FatalErrorMessage("The @VolumeProfile can be used only for NonIsothermal-UserDefinedVolume reactors");

				double value;
				std::string units;

				dictionaries(main_dictionary_name_).ReadMeasure("@PressureCoefficient", value, units);
				if (units == "Pa/s")		 batchreactor_volumeprofile->SetPressureCoefficient(value);
				else if (units == "bar/s")   batchreactor_volumeprofile->SetPressureCoefficient(value*1.e5);
				else if (units == "atm/s")   batchreactor_volumeprofile->SetPressureCoefficient(value*101325.);
				else if (units == "Pa/ms")	 batchreactor_volumeprofile->SetPressureCoefficient(value*1000.);
				else if (units == "bar/ms")  batchreactor_volumeprofile->SetPressureCoefficient(value*1.e5*1000.);
				else if (units == "atm/ms")  batchreactor_volumeprofile->SetPressureCoefficient(value*101325.*1000.);
				else OpenSMOKE::FatalErrorMessage("Unknown pressure coefficient units. Available: Pa/s || Pa/ms || bar/s || bar/ms || atm/s || atm/ms");


				volume_profile = true;
			}

			// Read volume profile
			if (dictionaries(main_dictionary_name_).CheckOption("@VolumeProfile") == true)
			{
				batchreactor_volumeprofile = new OpenSMOKE::BatchReactor_VolumeProfile();

				std::string name_of_profile_subdictionary;

				if (type_ != OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
					OpenSMOKE::FatalErrorMessage("The @VolumeProfile can be used only for NonIsothermal-UserDefinedVolume reactors");

				dictionaries(main_dictionary_name_).ReadDictionary("@VolumeProfile", name_of_profile_subdictionary);

				OpenSMOKE::OpenSMOKEVectorDouble x, y;
				std::string x_variable, y_variable;
				GetXYProfileFromDictionary(dictionaries(name_of_profile_subdictionary), x, y, x_variable, y_variable);

				if (x_variable != "time")
					OpenSMOKE::FatalErrorMessage("The @VolumeProfile must be defined versus the time");
				if (y_variable != "volume")
					OpenSMOKE::FatalErrorMessage("The @VolumeProfile must be defined versus the volume");

				if (y[1] != volume)
					OpenSMOKE::FatalErrorMessage("The @VolumeProfile and the @Volume options must be consistent");

				batchreactor_volumeprofile->SetProfile(x, y);

				volume_profile = true;
			}

			// Read pressure profile
			if (dictionaries(main_dictionary_name_).CheckOption("@PressureProfile") == true)
			{
				batchreactor_volumeprofile = new OpenSMOKE::BatchReactor_VolumeProfile();

				std::string name_of_profile_subdictionary;

				if (type_ != OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
					OpenSMOKE::FatalErrorMessage("The @PressureProfile can be used only for NonIsothermal-UserDefinedVolume reactors");

				dictionaries(main_dictionary_name_).ReadDictionary("@PressureProfile", name_of_profile_subdictionary);

				OpenSMOKE::OpenSMOKEVectorDouble x, y;
				std::string x_variable, y_variable;
				GetXYProfileFromDictionary(dictionaries(name_of_profile_subdictionary), x, y, x_variable, y_variable);

				if (x_variable != "time")
					OpenSMOKE::FatalErrorMessage("The @PressureProfile must be defined versus the time");
				if (y_variable != "pressure")
					OpenSMOKE::FatalErrorMessage("The @PressureProfile must be defined versus the volume");

				if (y[1] != P_Pa)
					OpenSMOKE::FatalErrorMessage("The @PressureProfile and the initial pressure of mixture must be consistent");

				batchreactor_volumeprofile->SetPressureProfile(x, y);

				volume_profile = true;
			}
		}

		// Polimi soot
		{
			polimi_soot_ = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML_);
			std::string name_of_polimisoot_analyzer_subdictionary;
			if (dictionaries(main_dictionary_name_).CheckOption("@PolimiSoot") == true)
			{
				dictionaries(main_dictionary_name_).ReadDictionary("@PolimiSoot", name_of_polimisoot_analyzer_subdictionary);
				polimi_soot_->SetupFromDictionary(dictionaries(name_of_polimisoot_analyzer_subdictionary));
			}
		}

/*		// Solve the ODE system: NonIsothermal, Constant Volume
		if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV)
		{
			batch_nonisothermal_constantv_ = new OpenSMOKE::BatchReactor_NonIsothermal_ConstantVolume(*thermodynamicsMapXML, *kineticsMapXML,
					*ode_parameters_, *batch_options_, *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_,  *polimi_soot_, volume, T, P_Pa, omega,
					global_thermal_exchange_coefficient, exchange_area, T_environment);
		//	batch_nonisothermal_constantv_->Solve(tStart_, tEnd_);
		}

		// Solve the ODE system: NonIsothermal, Volume assigned according to a specified law
		if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
		{
			batch_nonisothermal_userdefinedv_ = new OpenSMOKE::BatchReactor_NonIsothermal_UserDefinedVolume(*thermodynamicsMapXML, *kineticsMapXML,
					*ode_parameters_, *batch_options_, *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_,  *polimi_soot_, volume, T, P_Pa, omega, tStart_,
					global_thermal_exchange_coefficient, exchange_area, T_environment);

			batch_nonisothermal_userdefinedv_->SetVolumeProfile(*batchreactor_volumeprofile);

		//	batch_nonisothermal_userdefinedv_->Solve(tStart_, tEnd_);
		}

		// Solve the ODE system: Isothermal, Constant Volume
		if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV)
		{
			batch_isothermal_constantv_ = new OpenSMOKE::BatchReactor_Isothermal_ConstantVolume(*thermodynamicsMapXML, *kineticsMapXML,
					*ode_parameters_, *batch_options_, *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_,  *polimi_soot_, volume, T, P_Pa, omega);
		//	batch_isothermal_constantv_->Solve(tStart_, tEnd_);
		}

		// Solve the ODE system: NonIsothermal, Constant Pressure
		if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP)
		{
			batch_nonisothermal_constantp_ = new OpenSMOKE::BatchReactor_NonIsothermal_ConstantPressure(*thermodynamicsMapXML, *kineticsMapXML,
				*ode_parameters_, *batch_options_, *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_, *polimi_soot_, volume, T, P_Pa, omega,
				global_thermal_exchange_coefficient, exchange_area, T_environment);
		}

		// Solve the ODE system: Isothermal, Constant Pressure
		if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP)
		{
			batch_isothermal_constantp_ = new OpenSMOKE::BatchReactor_Isothermal_ConstantPressure(*thermodynamicsMapXML, *kineticsMapXML,
				*ode_parameters_, *batch_options_, *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_, *polimi_soot_, volume, T, P_Pa, omega);
		}
*/
	}

void BatchReactor_Plugin::Update_and_Solve_Batch(OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML, OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML)
{
	// Solve the ODE system: NonIsothermal, Constant Volume
	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV)
	{
		batch_nonisothermal_constantv_ = new OpenSMOKE::BatchReactor_NonIsothermal_ConstantVolume(*thermodynamicsMapXML, *kineticsMapXML,
				*ode_parameters_, *batch_options_, *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_, *idt, *polimi_soot_, volume, T, P_Pa, omega,
				global_thermal_exchange_coefficient, exchange_area, T_environment);
		batch_nonisothermal_constantv_->Solve(tStart_, tEnd_);
		
	}
		// Solve the ODE system: NonIsothermal, Volume assigned according to a specified law
	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
	{
		batch_nonisothermal_userdefinedv_ = new OpenSMOKE::BatchReactor_NonIsothermal_UserDefinedVolume(*thermodynamicsMapXML, *kineticsMapXML,
				*ode_parameters_, *batch_options_, *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_, *idt, *polimi_soot_, volume, T, P_Pa, omega, tStart_,
				global_thermal_exchange_coefficient, exchange_area, T_environment);

		batch_nonisothermal_userdefinedv_->SetVolumeProfile(*batchreactor_volumeprofile);
		batch_nonisothermal_userdefinedv_->Solve(tStart_, tEnd_);
		
	}

	// Solve the ODE system: Isothermal, Constant Volume
	if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV)
	{
		batch_isothermal_constantv_ = new OpenSMOKE::BatchReactor_Isothermal_ConstantVolume(*thermodynamicsMapXML, *kineticsMapXML,
				*ode_parameters_, *batch_options_, *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_, *idt,  *polimi_soot_, volume, T, P_Pa, omega);
		batch_isothermal_constantv_->Solve(tStart_, tEnd_);
		
	}

	// Solve the ODE system: NonIsothermal, Constant Pressure
	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP)
	{
		batch_nonisothermal_constantp_ = new OpenSMOKE::BatchReactor_NonIsothermal_ConstantPressure(*thermodynamicsMapXML, *kineticsMapXML,
			*ode_parameters_, *batch_options_, *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_, *idt, *polimi_soot_, volume, T, P_Pa, omega,
			global_thermal_exchange_coefficient, exchange_area, T_environment);
		batch_nonisothermal_constantp_->Solve(tStart_, tEnd_);
		
	}

	// Solve the ODE system: Isothermal, Constant Pressure
	if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP)
	{
		batch_isothermal_constantp_ = new OpenSMOKE::BatchReactor_Isothermal_ConstantPressure(*thermodynamicsMapXML, *kineticsMapXML,
			*ode_parameters_, *batch_options_, *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_, *idt, *polimi_soot_, volume, T, P_Pa, omega);
		batch_isothermal_constantp_->Solve(tStart_, tEnd_);
		
	}

}

double BatchReactor_Plugin::Solve_tau(std::string tau_calc_type_temp)
{
	double tau_ign_temp;

	if(tau_calc_type_temp == "Temperature")
	{
		tau_ign_temp = idt->temperature_increase_tau();
	} 
	else if(tau_calc_type_temp == "MaxTemp")
	{
		tau_ign_temp = idt->temperature_max_tau();
	}
	else if(tau_calc_type_temp == "MaxDTemp")
	{
		tau_ign_temp = idt->temperature_slope_tau();
	}
	else if(tau_calc_type_temp == "MaxP")
	{
		tau_ign_temp = idt->pressure_max_tau();
	}
	else if(tau_calc_type_temp == "MaxDp")
	{
		tau_ign_temp =idt->pressure_slope_tau();
	}
	else if(tau_calc_type_temp == "PressureIncrease")
	{
		tau_ign_temp = idt->pressure_increase_tau();
	}
	else if(tau_calc_type_temp.substr(0,12) == "MaxIntercept")
	{
		std::vector<unsigned int> species_index_temp;
		species_index_temp = idt->species_intercept_max_index();
		if (std::find(species_index_temp.begin(),species_index_temp.end(),thermodynamicsMapXML_->IndexOfSpecies(tau_calc_type_temp.substr(12,tau_calc_type_temp.length()))-1) == species_index_temp.end())
		{
			OpenSMOKE::FatalErrorMessage("Species " + tau_calc_type_temp.substr(12,tau_calc_type_temp.length()) + " not specified in the OpenSMOKE++ input file properly!");
		}
		int index_Species = std::distance(species_index_temp.begin(),std::find(species_index_temp.begin(),species_index_temp.end(),thermodynamicsMapXML_->IndexOfSpecies(tau_calc_type_temp.substr(12,tau_calc_type_temp.length()))-1));
		//std::cout<<"index_Species: "<<index_Species<<std::endl;
		tau_ign_temp = idt->species_intercept_max_tau(index_Species);
	} 
	else if(tau_calc_type_temp == "MaxOH")
	{
		std::vector<double> tau_ign_temp_vec;
		std::vector<unsigned int> species_index_temp;
		species_index_temp = idt->species_index();
		int index_OH = std::distance(species_index_temp.begin(),std::find(species_index_temp.begin(),species_index_temp.end(),thermodynamicsMapXML_->IndexOfSpecies("OH")))-1;
		tau_ign_temp_vec = idt->species_max_tau();
		tau_ign_temp = tau_ign_temp_vec[index_OH];
	} 
	else if(tau_calc_type_temp == "MaxDOH")
	{
		std::vector<double> tau_ign_temp_vec;
		std::vector<unsigned int> species_index_temp;
		species_index_temp = idt->species_index();
		int index_OH = std::distance(species_index_temp.begin(),std::find(species_index_temp.begin(),species_index_temp.end(),thermodynamicsMapXML_->IndexOfSpecies("OH")))-1;
		tau_ign_temp_vec = idt->species_slope_tau();
		tau_ign_temp = tau_ign_temp_vec[index_OH];
	}
    else if(tau_calc_type_temp == "MaxOH*")
    {
        std::vector<double> tau_ign_temp_vec;
        std::vector<unsigned int> species_index_temp;
        species_index_temp = idt->species_index();
        int index_OH_St = std::distance(species_index_temp.begin(),std::find(species_index_temp.begin(),species_index_temp.end(),thermodynamicsMapXML_->IndexOfSpecies("OH*")))-1;
        tau_ign_temp_vec = idt->species_max_tau();
        tau_ign_temp = tau_ign_temp_vec[index_OH_St];
    }
    else if(tau_calc_type_temp == "MaxDOH*")
    {
        std::vector<double> tau_ign_temp_vec;
        std::vector<unsigned int> species_index_temp;
        species_index_temp = idt->species_index();
        int index_OH_St = std::distance(species_index_temp.begin(),std::find(species_index_temp.begin(),species_index_temp.end(),thermodynamicsMapXML_->IndexOfSpecies("OH*")))-1;
        tau_ign_temp_vec = idt->species_slope_tau();
        tau_ign_temp = tau_ign_temp_vec[index_OH_St];
    }
	else 
	{
		OpenSMOKE::FatalErrorMessage("Specified method for calculating ignition delay time has not been implimented yet! (" + tau_calc_type_temp + ")");
	}

	//ciao
	//idt->Reset();

	clean_up();


	return tau_ign_temp;
}


std::vector<double> BatchReactor_Plugin::Solve_Species(std::vector<double> Abscissa, std::string Species)
{
	
	OpenSMOKE::OpenSMOKEVectorDouble omega_Final(thermodynamicsMapXML_->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble x_Final(thermodynamicsMapXML_->NumberOfSpecies());

	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV)
	{
		std::vector<double> collector;

		for (int i=0; i<Abscissa.size(); i ++) //Abscissa.size()
		{
			double interpolated_point;
			for (int j=0; batch_nonisothermal_constantv_->time_vector.size(); j++)
			{
				
				if(Abscissa[i] <= batch_nonisothermal_constantv_->time_vector[j])
				{
					interpolated_point = batch_nonisothermal_constantv_->species_matrix[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species)) + (batch_nonisothermal_constantv_->species_matrix[j](thermodynamicsMapXML_->IndexOfSpecies(Species))-batch_nonisothermal_constantv_->species_matrix[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species)))/(batch_nonisothermal_constantv_->time_vector[j]-batch_nonisothermal_constantv_->time_vector[j-1])*(Abscissa[i]-batch_nonisothermal_constantv_->time_vector[j-1]);
					//std::cout<<"Index of the time vector where the condition is satisfied is: "<< j << std::endl;
					collector.push_back(interpolated_point);

					//std::cout<<"The molar fraction of the species "<< Species <<" is " << interpolated_point << std::endl;
					break;
				}				
			}
		}
		
		clean_up();
		return collector;
	}

	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
	{
		std::vector<double> collector;

		for (int i=0; i<Abscissa.size(); i ++) //Abscissa.size()
		{
			double interpolated_point;

			for (int j=0; batch_nonisothermal_constantv_->time_vector.size(); j++)
			{
				
				if(Abscissa[i] <= batch_nonisothermal_userdefinedv_->time_vector[j])
				{
					interpolated_point = batch_nonisothermal_userdefinedv_->species_matrix[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species)) + (batch_nonisothermal_userdefinedv_->species_matrix[j](thermodynamicsMapXML_->IndexOfSpecies(Species))-batch_nonisothermal_userdefinedv_->species_matrix[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species)))/(batch_nonisothermal_userdefinedv_->time_vector[j]-batch_nonisothermal_userdefinedv_->time_vector[j-1])*(Abscissa[i]-batch_nonisothermal_userdefinedv_->time_vector[j-1]);
					//std::cout<<"Index of the time vector where the condition is satisfied is: "<< j << std::endl;
					collector.push_back(interpolated_point);
					//std::cout<<"The molar fraction of the species "<< Species <<" is " << batch_nonisothermal_constantv_->species_matrix[j](thermodynamicsMapXML_->IndexOfSpecies(Species)) << std::endl;
					break;
				}				
			}
		}

		clean_up();
		return collector;
	}

	if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV)
	{
		std::vector<double> collector;

		for (int i=0; i<Abscissa.size(); i ++) //Abscissa.size()
		{
			double interpolated_point;

			for (int j=0; batch_isothermal_constantv_->time_vector.size(); j++)
			{
				
				if(Abscissa[i] <= batch_isothermal_constantv_->time_vector[j])
				{
					interpolated_point = batch_isothermal_constantv_->species_matrix[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species)) + (batch_isothermal_constantv_->species_matrix[j](thermodynamicsMapXML_->IndexOfSpecies(Species))-batch_isothermal_constantv_->species_matrix[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species)))/(batch_isothermal_constantv_->time_vector[j]-batch_isothermal_constantv_->time_vector[j-1])*(Abscissa[i]-batch_isothermal_constantv_->time_vector[j-1]);
					//std::cout<<"Index of the time vector where the condition is satisfied is: "<< j << std::endl;
					collector.push_back(interpolated_point);
					//std::cout<<"The molar fraction of the species "<< Species <<" is " << batch_nonisothermal_constantv_->species_matrix[j](thermodynamicsMapXML_->IndexOfSpecies(Species)) << std::endl;
					break;
				}				
			}
		}
		clean_up();
		return collector;
	}

	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP)
	{
		std::vector<double> collector;

		for (int i=0; i<Abscissa.size(); i ++) //Abscissa.size()
		{
			double interpolated_point;

			for (int j=0; batch_nonisothermal_constantp_->time_vector.size(); j++)
			{
				
				if(Abscissa[i] <= batch_nonisothermal_constantp_->time_vector[j])
				{
					interpolated_point = batch_nonisothermal_constantp_->species_matrix[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species)) + (batch_nonisothermal_constantp_->species_matrix[j](thermodynamicsMapXML_->IndexOfSpecies(Species))-batch_nonisothermal_constantp_->species_matrix[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species)))/(batch_nonisothermal_constantp_->time_vector[j]-batch_nonisothermal_constantp_->time_vector[j-1])*(Abscissa[i]-batch_nonisothermal_constantp_->time_vector[j-1]);
					//std::cout<<"Index of the time vector where the condition is satisfied is: "<< j << std::endl;
					collector.push_back(interpolated_point);
					//std::cout<<"The molar fraction of the species "<< Species <<" is " << batch_nonisothermal_constantv_->species_matrix[j](thermodynamicsMapXML_->IndexOfSpecies(Species)) << std::endl;
					break;
				}				
			}
		}
		clean_up();
		return collector;
	}

	if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP)
	{
		std::vector<double> collector;

		for (int i=0; i<Abscissa.size(); i ++) //Abscissa.size()
		{
			double interpolated_point;

			for (int j=0; batch_isothermal_constantp_->time_vector.size(); j++)
			{
				
				if(Abscissa[i] <= batch_isothermal_constantp_->time_vector[j])
				{
					interpolated_point = batch_isothermal_constantp_->species_matrix[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species)) + (batch_isothermal_constantp_->species_matrix[j](thermodynamicsMapXML_->IndexOfSpecies(Species))-batch_isothermal_constantp_->species_matrix[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species)))/(batch_isothermal_constantp_->time_vector[j]-batch_isothermal_constantp_->time_vector[j-1])*(Abscissa[i]-batch_isothermal_constantp_->time_vector[j-1]);
					//std::cout<<"Index of the time vector where the condition is satisfied is: "<< j << std::endl;
					collector.push_back(interpolated_point);
					//std::cout<<"The molar fraction of the species "<< Species <<" is " << batch_nonisothermal_constantv_->species_matrix[j](thermodynamicsMapXML_->IndexOfSpecies(Species)) << std::endl;
					break;
				}				
			}
		}
		clean_up();
		return collector;
	}

}

std::vector<std::vector<double>> BatchReactor_Plugin::Solve_Multiple_Species(std::vector<std::vector<std::vector<double>>> &Exp_data_temp, std::vector<std::string> Species_vec)
{

	OpenSMOKE::OpenSMOKEVectorDouble omega_Final(thermodynamicsMapXML_->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble x_Final(thermodynamicsMapXML_->NumberOfSpecies());
	
	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV)
	{
		std::vector<std::vector<double>> collector;
		collector.resize(Species_vec.size());
		for(int z=0; z<Species_vec.size(); z++)
		{
			for (int i=0; i<Exp_data_temp[z][0].size(); i ++) //Abscissa.size()
			{
				double temporal_data = Interpolate(batch_nonisothermal_constantv_->time_vector, batch_nonisothermal_constantv_->species_matrix, Exp_data_temp[z][0][i], Species_vec[z]);

				if (temporal_data == 10000000) {
					collector[z].push_back(collector[z][i-1]);
					std::cout<< "This happen in the position " << i <<std::endl;
				} else {
					collector[z].push_back(temporal_data);
				}
					
			}
		}
		clean_up();
		return collector;
	}
	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
	{
		std::vector<std::vector<double>> collector;
		collector.resize(Species_vec.size());
		for(int z=0; z<Species_vec.size(); z++)
		{
			for (int i=0; i<Exp_data_temp[z][0].size(); i ++) //Abscissa.size()
			{
				double temporal_data = Interpolate(batch_nonisothermal_userdefinedv_->time_vector, batch_nonisothermal_userdefinedv_->species_matrix, Exp_data_temp[z][0][i], Species_vec[z]);

				if (temporal_data == 10000000) {
					collector[z].push_back(collector[z][i-1]);
					std::cout<< "This happen in the position " << i <<std::endl;
				} else {
					collector[z].push_back(temporal_data);
				}

			}
		}
		clean_up();
		return collector;
	}
	if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV)
	{
		std::vector<std::vector<double>> collector;
		collector.resize(Species_vec.size());
		for(int z=0; z<Species_vec.size(); z++)
		{
			for (int i=0; i<Exp_data_temp[z][0].size(); i ++) //Abscissa.size()
			{
				double temporal_data = Interpolate(batch_isothermal_constantv_->time_vector, batch_isothermal_constantv_->species_matrix, Exp_data_temp[z][0][i], Species_vec[z]);

				if (temporal_data == 10000000) {
					collector[z].push_back(collector[z][i-1]);
					std::cout<< "This happen in the position " << i <<std::endl;
				} else {
					collector[z].push_back(temporal_data);
				}

			}
		}
		clean_up();
		return collector;
	}
	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP)
	{
		std::vector<std::vector<double>> collector;
		collector.resize(Species_vec.size());
		for(int z=0; z<Species_vec.size(); z++)
		{
			for (int i=0; i<Exp_data_temp[z][0].size(); i ++) //Abscissa.size()
			{
				double temporal_data = Interpolate(batch_nonisothermal_constantp_->time_vector, batch_nonisothermal_constantp_->species_matrix, Exp_data_temp[z][0][i], Species_vec[z]);

				if (temporal_data == 10000000) {
					collector[z].push_back(collector[z][i-1]);
					std::cout<< "This happen in the position " << i <<std::endl;
				} else {
					collector[z].push_back(temporal_data);
				}

			}
		}
		clean_up();
		return collector;
	}
	if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP)
	{
		std::vector<std::vector<double>> collector;
		collector.resize(Species_vec.size());
		for(int z=0; z<Species_vec.size(); z++)
		{
			for (int i=0; i<Exp_data_temp[z][0].size(); i ++) //Abscissa.size()
			{
				double temporal_data = Interpolate(batch_isothermal_constantp_->time_vector, batch_isothermal_constantp_->species_matrix, Exp_data_temp[z][0][i], Species_vec[z]);

				if (temporal_data == 10000000) {
					collector[z].push_back(collector[z][i-1]);
					std::cout<< "This happen in the position " << i <<std::endl;
				} else {
					collector[z].push_back(temporal_data);
				}
			}
		}

		clean_up();
		return collector;

	}
}


double BatchReactor_Plugin::Interpolate(std::vector<double> Time_vec_temp, std::vector<OpenSMOKE::OpenSMOKEVectorDouble> Species_matrix_temp, double Abscissa_temp, std::string Species_name)
{
	double interpolated_point;

	if (Abscissa_temp > Time_vec_temp[Time_vec_temp.size()-1])
	{
		std::cout<<"WARNING_OptiSMOKE: The time (" + std::to_string(Abscissa_temp) + ") that you specified into the EXP_DATA is outside the runtime of this simulation (" + std::to_string(Time_vec_temp[Time_vec_temp.size()-1]) + ")! Extrapolation was not possible! Please be sure that this is what you want."<<std::endl;
		interpolated_point=10000000;

	} else {

		for (int j=0; j<Time_vec_temp.size(); j++)
		{
			if(Abscissa_temp <= Time_vec_temp[j])
			{
				interpolated_point = Species_matrix_temp[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species_name)) + (Species_matrix_temp[j](thermodynamicsMapXML_->IndexOfSpecies(Species_name))-Species_matrix_temp[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species_name)))/(Time_vec_temp[j]-Time_vec_temp[j-1])*(Abscissa_temp-Time_vec_temp[j-1]);
				//std::cout<<"Index of the time vector where the condition is satisfied is: "<< j << std::endl;
				//std::cout<<"The molar fraction of the species "<< Species_name <<" is " << interpolated_point <<" at abscissa "<<Abscissa_temp<< std::endl;
				break;
			}
		}
		
	}

	return interpolated_point;
}

void BatchReactor_Plugin::clean_up()
{
	delete batch_options_ ;
	delete ode_parameters_ ;	
	delete sensitivity_options_ ;
	delete onTheFlyROPA_ ;
	delete onTheFlyCEMA_ ;
	delete on_the_fly_post_processing_ ;
	delete polimi_soot_;

	if(volume_profile)
	{
		delete batchreactor_volumeprofile;
	}


	delete idt;

	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV)
	{
		delete batch_nonisothermal_constantv_ ;
	}
	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME)
	{
		delete batch_nonisothermal_userdefinedv_ ;
	}
	if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV)
	{
		delete batch_isothermal_constantv_ ;
	}
	if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP)
	{
		delete batch_nonisothermal_constantp_ ;
	}
	if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP)
	{
		delete batch_isothermal_constantp_ ;
	}
}

}
