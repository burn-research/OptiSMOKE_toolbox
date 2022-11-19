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
void PlugFlowReactor_Plugin::Setup(const std::string input_file_name,
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML,
		OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML)
{
		// Pointers
		thermodynamicsMapXML_ = thermodynamicsMapXML;
		kineticsMapXML_ = kineticsMapXML;

		// Defines the grammar rules
		OpenSMOKE::Grammar_PlugFlowReactor grammar_plugflowreactor;
		
		// Define the dictionaries
		std::string main_dictionary_name_ = "PlugFlowReactor";
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(main_dictionary_name_).SetGrammar(grammar_plugflowreactor);

		// Read plug-flow reactor type
		{
			std::string value;
			if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true)
			{
				dictionaries(main_dictionary_name_).ReadString("@Type", value);
				if (value == "Isothermal")	type_ = OpenSMOKE::PLUGFLOW_REACTOR_ISOTHERMAL;
				else if (value == "NonIsothermal")	type_ = OpenSMOKE::PLUGFLOW_REACTOR_NONISOTHERMAL;
				else OpenSMOKE::FatalErrorMessage("Unknown plug flow reactor type: " + value);

			}
		}

		constant_pressure = true;
		{
			if (dictionaries(main_dictionary_name_).CheckOption("@ConstantPressure") == true)
				dictionaries(main_dictionary_name_).ReadBool("@ConstantPressure", constant_pressure);
		}

		// Read initial conditions
		{
			std::string name_of_gas_status_subdictionary;
			if (dictionaries(main_dictionary_name_).CheckOption("@InletStatus") == true)
				dictionaries(main_dictionary_name_).ReadDictionary("@InletStatus", name_of_gas_status_subdictionary);

			GetGasStatusFromDictionary(dictionaries(name_of_gas_status_subdictionary), *thermodynamicsMapXML_, T, P_Pa, omega);
		}

		// Read integration limits
		end_value_ = 0.;
		time_independent_variable = false;

		// Read residence time
		{

			if (dictionaries(main_dictionary_name_).CheckOption("@ResidenceTime") == true)
			{
				double value;
				std::string units;
				dictionaries(main_dictionary_name_).ReadMeasure("@ResidenceTime", value, units);
				if (units == "s")		  end_value_ = value;
				else if (units == "ms")   end_value_ = value / 1000.;
				else if (units == "min")  end_value_ = value * 60.;
				else if (units == "h")    end_value_ = value * 3600.;
				else OpenSMOKE::FatalErrorMessage("Unknown time units");
				time_independent_variable = true;
			}
		}

		// Read length
		{
			if (dictionaries(main_dictionary_name_).CheckOption("@Length") == true)
			{
				double value;
				std::string units;
				dictionaries(main_dictionary_name_).ReadMeasure("@Length", value, units);
				if (units == "m")		  end_value_ = value;
				else if (units == "cm")   end_value_ = value / 100.;
				else if (units == "mm")   end_value_ = value / 1000.;
				else OpenSMOKE::FatalErrorMessage("Unknown length units");
				time_independent_variable = false;
			}
		}

		// Read diameter
		double diameter = 0.;
		{
			if (dictionaries(main_dictionary_name_).CheckOption("@Diameter") == true)
			{
				double value;
				std::string units;
				dictionaries(main_dictionary_name_).ReadMeasure("@Diameter", value, units);
				if (units == "m")		  diameter = value;
				else if (units == "cm")   diameter = value / 100.;
				else if (units == "mm")   diameter = value / 1000.;
				else OpenSMOKE::FatalErrorMessage("Unknown length units");
			}
		}

		// Read velocity
		velocity = 0.;
		{
			if (dictionaries(main_dictionary_name_).CheckOption("@Velocity") == true)
			{
				double value;
				std::string units;
				dictionaries(main_dictionary_name_).ReadMeasure("@Velocity", value, units);
				if (units == "m/s")		   velocity = value;
				else if (units == "cm/s")  velocity = value / 1.e2;
				else if (units == "mm/s")  velocity = value / 1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown velocity units");
			}

			if (dictionaries(main_dictionary_name_).CheckOption("@VolumetricFlowRate") == true)
			{
				double value;
				std::string units;
				double volumetric_flow_rate;
				dictionaries(main_dictionary_name_).ReadMeasure("@VolumetricFlowRate", value, units);
				if (units == "m3/s")	    volumetric_flow_rate = value;
				else if (units == "cm3/s")  volumetric_flow_rate = value / 1.e6;
				else if (units == "mm3/s")  volumetric_flow_rate = value / 1.e9;
				else OpenSMOKE::FatalErrorMessage("Unknown volumetric flow rate units");

				const double A = PhysicalConstants::pi_over_4*diameter*diameter;
				velocity = volumetric_flow_rate / A;
			}

			if (dictionaries(main_dictionary_name_).CheckOption("@MassFlowRate") == true)
			{
				double value;
				std::string units;
				double mass_flow_rate;
				dictionaries(main_dictionary_name_).ReadMeasure("@MassFlowRate", value, units);
				if (units == "kg/s")	  mass_flow_rate = value;
				else if (units == "g/s")  mass_flow_rate = value / 1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown mass flow rate units");

				const double A = PhysicalConstants::pi_over_4*diameter*diameter;
				const double MW = thermodynamicsMapXML_->MolecularWeight_From_MassFractions(omega.GetHandle());
				const double rho = P_Pa * MW / PhysicalConstants::R_J_kmol / T;
				velocity = mass_flow_rate / A / rho;
			}

			if (dictionaries(main_dictionary_name_).CheckOption("@MoleFlowRate") == true)
			{
				double value;
				std::string units;
				double mole_flow_rate;
				dictionaries(main_dictionary_name_).ReadMeasure("@MoleFlowRate", value, units);
				if (units == "kmol/s")		mole_flow_rate = value;
				else if (units == "mol/s")  mole_flow_rate = value / 1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown mole flow rate units");

				const double A = PhysicalConstants::pi_over_4*diameter*diameter;
				const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
				velocity = mole_flow_rate / A / cTot;
			}
		}

		// Read ratio between cross section and perimeter
		cross_section_over_perimeter = 0.;
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@CrossSectionOverPerimeter") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@CrossSectionOverPerimeter", value, units);
				if (units == "m")        cross_section_over_perimeter = value;
				else if (units == "dm")  cross_section_over_perimeter = value / 1.e1;
				else if (units == "cm")  cross_section_over_perimeter = value / 1.e2;
				else if (units == "mm")  cross_section_over_perimeter = value / 1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown length units");
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

		// Profile
		temperature_profile = false;
		{
			std::string name_of_gas_status_subdictionary;
			if (dictionaries(main_dictionary_name_).CheckOption("@TemperatureProfile") == true)
			{
				if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_NONISOTHERMAL)
					OpenSMOKE::FatalErrorMessage("The @TemperatureProfile cannot be used for adiabatic reactors");

				dictionaries(main_dictionary_name_).ReadDictionary("@TemperatureProfile", name_of_gas_status_subdictionary);

				OpenSMOKE::OpenSMOKEVectorDouble x, y;
				std::string x_variable, y_variable;
				GetXYProfileFromDictionary(dictionaries(name_of_gas_status_subdictionary), x, y, x_variable, y_variable);

				if (x_variable == "time" && time_independent_variable == false)
					OpenSMOKE::FatalErrorMessage("The @TemperatureProfile must be defined versus the time");
				if (x_variable == "length" && time_independent_variable == true)
					OpenSMOKE::FatalErrorMessage("The @TemperatureProfile must be defined versus the length");
				if (std::fabs(y[1] - T) > 1e-6)
					OpenSMOKE::FatalErrorMessage("Please check the @TemperatureProfile: the inlet values do not match");
				if (std::fabs(x[x.Size()] - end_value_) / end_value_ > 1.e-6)
					OpenSMOKE::FatalErrorMessage("Please check the @TemperatureProfile: the x domains do not match");

				temperature_profile = true;
				profile = new OpenSMOKE::PlugFlowReactor_Profile(x, y, x_variable);
			}
		}

		// Options
		{
			plugflow_options_ = new OpenSMOKE::PlugFlowReactor_Options();
			if (dictionaries(main_dictionary_name_).CheckOption("@Options") == true)
			{
				std::string name_of_options_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@Options", name_of_options_subdictionary);
				plugflow_options_->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
			}
			plugflow_options_->SetVerboseVideo(false);
			plugflow_options_->SetVerboseOutput(false);
			plugflow_options_->SetSensitivityAnalysis(false);
			//plugflow_options_->output_species();
		}

		// ODE Parameters
		{
			ode_parameters_= new OpenSMOKE::ODE_Parameters();
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

				//plugflow_options_->SetSensitivityAnalysis(true);
				//sensitivity_options_->SetupFromDictionary(dictionaries(name_of_sensitivity_options_subdictionary));
			}
		}

		// On the fly ROPA
		{
			onTheFlyROPA_ = new OpenSMOKE::OnTheFlyROPA(*thermodynamicsMapXML_, *kineticsMapXML_);
			if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyROPA") == true)
			{
				std::string name_of_options_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyROPA", name_of_options_subdictionary);
				//	onTheFlyROPA_->SetupFromDictionary(dictionaries(name_of_options_subdictionary), path_kinetics_output);
			}
		}

		// On the fly PostProcessing
		{
			on_the_fly_post_processing_ = new OpenSMOKE::OnTheFlyPostProcessing(*thermodynamicsMapXML_, *kineticsMapXML_, plugflow_options_->output_path());

			if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyPostProcessing") == true)
			{
				std::string name_of_options_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyPostProcessing", name_of_options_subdictionary);
				on_the_fly_post_processing_->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
			}
		}

		// Ignition Delay Times
		idt = new OpenSMOKE::IgnitionDelayTimes_Analyzer();
		if (dictionaries(main_dictionary_name_).CheckOption("@IgnitionDelayTimes") == true)
		{
			if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_ISOTHERMAL)
				OpenSMOKE::FatalErrorMessage("The @IgnitionDelayTimes can be used only for NonIsothermal reactors");
	
			std::string name_of_idt_subdictionary;
			dictionaries(main_dictionary_name_).ReadDictionary("@IgnitionDelayTimes", name_of_idt_subdictionary);
			idt->SetupFromDictionary(dictionaries(name_of_idt_subdictionary), *thermodynamicsMapXML_);
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
}

void PlugFlowReactor_Plugin::Update_and_Solve_PFR(OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML, OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML)
{
	if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_ISOTHERMAL)
	{
		plugflow_isothermal_ = new OpenSMOKE::PlugFlowReactor_Isothermal(*thermodynamicsMapXML, *kineticsMapXML, *ode_parameters_, *plugflow_options_, *onTheFlyROPA_, 							*on_the_fly_post_processing_, *idt, *polimi_soot_, time_independent_variable, constant_pressure, velocity, T, P_Pa, omega);
			if (temperature_profile == true)
			plugflow_isothermal_->SetTemperatureProfile(*profile);
		plugflow_isothermal_->Solve(end_value_);
	}

	else if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_NONISOTHERMAL)
	{
		plugflow_non_isothermal_ = new OpenSMOKE::PlugFlowReactor_NonIsothermal(*thermodynamicsMapXML, *kineticsMapXML,
			*ode_parameters_, *plugflow_options_, *onTheFlyROPA_, *on_the_fly_post_processing_, *idt, *polimi_soot_, time_independent_variable, constant_pressure, velocity, T, P_Pa, omega,global_thermal_exchange_coefficient, cross_section_over_perimeter, T_environment);
		plugflow_non_isothermal_->Solve(end_value_);
	}
	
}

double PlugFlowReactor_Plugin::Solve_tau(std::string tau_calc_type_temp)
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
	else 
	{
		OpenSMOKE::FatalErrorMessage("Specified method for calculating ignition delay time has not been implimented yet! (" + tau_calc_type_temp + ")");
	}
	//idt->Reset();
	//std::cout<<"Before clean_up()"<<std::endl;
	clean_up();
	//std::cout<<"After clean_up()"<<std::endl;
	return tau_ign_temp;
}

double PlugFlowReactor_Plugin::Solve_Species(std::string Species)
{
	OpenSMOKE::OpenSMOKEVectorDouble omega_Final(thermodynamicsMapXML_->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble x_Final(thermodynamicsMapXML_->NumberOfSpecies());
	
	if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_ISOTHERMAL)
	{
		plugflow_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	else if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_NONISOTHERMAL)
	{
		plugflow_non_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	MW_Final = thermodynamicsMapXML_->MolecularWeight_From_MassFractions(omega_Final.GetHandle());
	thermodynamicsMapXML_->MoleFractions_From_MassFractions(x_Final.GetHandle(), MW_Final, omega_Final.GetHandle());
	Mole_frac_temp = x_Final(thermodynamicsMapXML_->IndexOfSpecies(Species));

	clean_up();
	return Mole_frac_temp;

}

std::vector<std::vector<double>> PlugFlowReactor_Plugin::Solve_Multipl_Species_time_profile(std::vector<std::vector<std::vector<double>>> &Exp_data_temp, std::vector<std::string> Species_vec)
{
	OpenSMOKE::OpenSMOKEVectorDouble omega_Final(thermodynamicsMapXML_->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble x_Final(thermodynamicsMapXML_->NumberOfSpecies());

	//get conversion vector from  exp
	std::vector<double> conversion_exp;
	for (int i=0; i<Exp_data_temp[0][0].size(); i ++){

		double temp = 1 - Exp_data_temp[0][1][i]/Exp_data_temp[0][1][0];
		conversion_exp.push_back(temp);
	}
	
	//compute the time with 50% conv in the experiments
	for (int i=0; i<conversion_exp.size(); i ++){
		
		if (conversion_exp[i] >= 0.5){
			// linear interpolation to get the time with 50% conv in experiments
			time_fifty_percent_conv_exp = Exp_data_temp[0][0][i-1]+(0.5-conversion_exp[i-1])*(Exp_data_temp[0][0][i]-Exp_data_temp[0][0][i-1])/(conversion_exp[i]-conversion_exp[i-1]);
			break;
		}
	}

	// WARNING: In the following we assume that the first species listed among the targets is the FUEL!
	if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_ISOTHERMAL)
	{
		//get conversion vector from sim
		std::vector<double> conversion_sim;
		for (int i=0; i<plugflow_isothermal_->time_vector.size(); i ++){

			double temp = 1 - plugflow_isothermal_->species_matrix[i](thermodynamicsMapXML_->IndexOfSpecies(Species_vec[0]))/plugflow_isothermal_->species_matrix[0](thermodynamicsMapXML_->IndexOfSpecies(Species_vec[0]));

			conversion_sim.push_back(temp);
		}

		//compute the time with 50% conv in the simulation
		for (int i=0; i<conversion_sim.size(); i ++){
			
			if (conversion_sim[i] >= 0.5){
				// linear interpolation to get the time with 50% conv in experiments
				time_fifty_percent_conv_sim = plugflow_isothermal_->time_vector[i-1]+(0.5-conversion_sim[i-1])*(plugflow_isothermal_->time_vector[i]-plugflow_isothermal_->time_vector[i-1])/(conversion_sim[i]-conversion_sim[i-1]);
				break;
			}
		}

		time_shift = time_fifty_percent_conv_sim - time_fifty_percent_conv_exp;

		std::cout<<"The time_fifty_percent_conv_sim is "<< time_fifty_percent_conv_sim <<" [s]"<<std::endl;
		std::cout<<"The time_fifty_percent_conv_exp is "<< time_fifty_percent_conv_exp <<" [s]"<<std::endl;
		std::cout<<"The time-shift is "<< time_shift <<" [s]"<<std::endl;

		std::vector<double> time_vector_new;
		time_vector_new.push_back(0.0);

		std::vector<OpenSMOKE::OpenSMOKEVectorDouble> species_matrix_new;
		species_matrix_new.push_back(plugflow_isothermal_->species_matrix[0]);

		if (time_shift >0){
			

			for (int i=0; i<plugflow_isothermal_->time_vector.size(); i ++){
				
				double temp_time = plugflow_isothermal_->time_vector[i]-time_shift;

				if (temp_time>0){
					time_vector_new.push_back(temp_time);
					species_matrix_new.push_back(plugflow_isothermal_->species_matrix[i]);
				}
			}


		} else {

			for (int i=0; i<plugflow_isothermal_->time_vector.size(); i ++){
				
				double temp_time = plugflow_isothermal_->time_vector[i]-time_shift;

				time_vector_new.push_back(temp_time);
				species_matrix_new.push_back(plugflow_isothermal_->species_matrix[i]);
			}

		}

		std::vector<std::vector<double>> collector;
		collector.resize(Species_vec.size());
		for(int z=0; z<Species_vec.size(); z++)
		{
			for (int i=0; i<Exp_data_temp[z][0].size(); i ++) //Abscissa.size()
			{
				if (Exp_data_temp[z][0][i] > time_vector_new[time_vector_new.size()-1])
				{
					time_vector_new.push_back(Exp_data_temp[z][0][i]);
					species_matrix_new.push_back(species_matrix_new[species_matrix_new.size()-1]);
				}

				collector[z].push_back(Interpolate(time_vector_new, species_matrix_new, Exp_data_temp[z][0][i], Species_vec[z]));
			}
		}

		clean_up();
		return collector;
	}
	else if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_NONISOTHERMAL)
	{
		//get conversion vector from sim
		std::vector<double> conversion_sim;
		for (int i=0; i<plugflow_non_isothermal_->time_vector.size(); i ++){
			double temp = 1 - plugflow_non_isothermal_->species_matrix[i](thermodynamicsMapXML_->IndexOfSpecies(Species_vec[0]))/plugflow_non_isothermal_->species_matrix[0](thermodynamicsMapXML_->IndexOfSpecies(Species_vec[0]));
			conversion_sim.push_back(temp);
		}

		//compute the time with 50% conv in the simulation
		for (int i=0; i<conversion_sim.size(); i ++){
			
			if (conversion_sim[i] >= 0.5){
				// linear interpolation to get the time with 50% conv in experiments
				time_fifty_percent_conv_sim = plugflow_non_isothermal_->time_vector[i-1]+(0.5-conversion_sim[i-1])*(plugflow_non_isothermal_->time_vector[i]-plugflow_non_isothermal_->time_vector[i-1])/(conversion_sim[i]-conversion_sim[i-1]);
				break;
			}
		}


		time_shift = time_fifty_percent_conv_sim - time_fifty_percent_conv_exp;

		std::cout<<"The time_fifty_percent_conv_sim is "<< time_fifty_percent_conv_sim <<" [s]"<<std::endl;
		std::cout<<"The time_fifty_percent_conv_exp is "<< time_fifty_percent_conv_exp <<" [s]"<<std::endl;
		std::cout<<"The time-shift is "<< time_shift <<" [s]"<<std::endl;

		std::vector<double> time_vector_new;
		time_vector_new.push_back(0.0);

		std::vector<OpenSMOKE::OpenSMOKEVectorDouble> species_matrix_new;
		species_matrix_new.push_back(plugflow_non_isothermal_->species_matrix[0]);

		if (time_shift >0){

			for (int i=0; i<plugflow_non_isothermal_->time_vector.size(); i ++){
				
				double temp_time = plugflow_non_isothermal_->time_vector[i]-time_shift;

				if (temp_time>0){
					time_vector_new.push_back(temp_time);
					species_matrix_new.push_back(plugflow_non_isothermal_->species_matrix[i]);
				}
			}


		} else {

			for (int i=0; i<plugflow_non_isothermal_->time_vector.size(); i ++){
				
				double temp_time = plugflow_non_isothermal_->time_vector[i]-time_shift;

				time_vector_new.push_back(temp_time);
				species_matrix_new.push_back(plugflow_non_isothermal_->species_matrix[i]);
			}

		}
		
		std::vector<std::vector<double>> collector;
		collector.resize(Species_vec.size());

		for(int z=0; z<Species_vec.size(); z++)
		{
			for (int i=0; i<Exp_data_temp[z][0].size(); i ++) //Abscissa.size()
			{
				if (Exp_data_temp[z][0][i] > time_vector_new[time_vector_new.size()-1])
				{
					time_vector_new.push_back(Exp_data_temp[z][0][i]);
					species_matrix_new.push_back(species_matrix_new[species_matrix_new.size()-1]);
				}

				collector[z].push_back(Interpolate(time_vector_new, species_matrix_new, Exp_data_temp[z][0][i], Species_vec[z]));
			}
		}

		clean_up();
		return collector;
	}

}

std::vector<double> PlugFlowReactor_Plugin::Solve_Multipl_Species_outlet(std::vector<std::string> Species_vec)
{

	OpenSMOKE::OpenSMOKEVectorDouble omega_Final(thermodynamicsMapXML_->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble x_Final(thermodynamicsMapXML_->NumberOfSpecies());
	
	if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_ISOTHERMAL)
	{
		plugflow_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	else if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_NONISOTHERMAL)
	{
		plugflow_non_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}


	MW_Final = thermodynamicsMapXML_->MolecularWeight_From_MassFractions(omega_Final.GetHandle());
	thermodynamicsMapXML_->MoleFractions_From_MassFractions(x_Final.GetHandle(), MW_Final, omega_Final.GetHandle());


	if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_ISOTHERMAL)
	{
		std::vector<double> Mole_frac_temp;
		for(int z=0; z<Species_vec.size(); z++)
		{
			Mole_frac_temp.push_back(x_Final(thermodynamicsMapXML_->IndexOfSpecies(Species_vec[z])));
		}

		clean_up();
		return Mole_frac_temp;
	}
	else if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_NONISOTHERMAL)
	{
		std::vector<double> Mole_frac_temp;
		for(int z=0; z<Species_vec.size(); z++)
		{
			Mole_frac_temp.push_back(x_Final(thermodynamicsMapXML_->IndexOfSpecies(Species_vec[z])));
		}

		clean_up();
		return Mole_frac_temp;
	}

}

double PlugFlowReactor_Plugin::Solve_Outlet_Conversion(std::string specie)
{
	OpenSMOKE::OpenSMOKEVectorDouble omega_Final(thermodynamicsMapXML_->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble x_Final(thermodynamicsMapXML_->NumberOfSpecies());
	
	if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_ISOTHERMAL)
	{
		plugflow_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	else if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_NONISOTHERMAL)
	{
		plugflow_non_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	MW_Final = thermodynamicsMapXML_->MolecularWeight_From_MassFractions(omega_Final.GetHandle());
	thermodynamicsMapXML_->MoleFractions_From_MassFractions(x_Final.GetHandle(), MW_Final, omega_Final.GetHandle());

	if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_ISOTHERMAL)
	{
		double outlet_conversion = 1 - plugflow_isothermal_->species_matrix[plugflow_isothermal_->time_vector.size()-1](thermodynamicsMapXML_->IndexOfSpecies(specie))/plugflow_isothermal_->species_matrix[0](thermodynamicsMapXML_->IndexOfSpecies(specie));
		clean_up();
		return outlet_conversion;
	}
	else if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_NONISOTHERMAL)
	{
		// Alberto la calcola cosi(omega0_[i] - omega_[i]) / omega0_[i] * 100.
		double outlet_conversion = 1 - plugflow_non_isothermal_->species_matrix[plugflow_non_isothermal_->time_vector.size()-1](thermodynamicsMapXML_->IndexOfSpecies(specie))/plugflow_non_isothermal_->species_matrix[0](thermodynamicsMapXML_->IndexOfSpecies(specie));
		clean_up();
		return outlet_conversion;
	}

}

double PlugFlowReactor_Plugin::Interpolate(std::vector<double> Time_vec_temp, std::vector<OpenSMOKE::OpenSMOKEVectorDouble> Species_matrix_temp, double Abscissa_temp, std::string Species_name)
{
	if (Abscissa_temp > Time_vec_temp[Time_vec_temp.size()-1])
	{
		OpenSMOKE::FatalErrorMessage("The time (" + std::to_string(Abscissa_temp) + ") is outside the runtime of the simulation (" + std::to_string(Time_vec_temp[Time_vec_temp.size()-1]) + ")! Extrapolation not possible! Please check OpenSMOKE++ input file to correct this.");
	}
	double interpolated_point;
	for (int j=1; j<Time_vec_temp.size(); j++)
	{

		if(Abscissa_temp <= Time_vec_temp[j])
		{
			interpolated_point = Species_matrix_temp[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species_name)) + (Species_matrix_temp[j](thermodynamicsMapXML_->IndexOfSpecies(Species_name))-Species_matrix_temp[j-1](thermodynamicsMapXML_->IndexOfSpecies(Species_name)))/(Time_vec_temp[j]-Time_vec_temp[j-1])*(Abscissa_temp-Time_vec_temp[j-1]);
			//std::cout<<"Index of the time vector where the condition is satisfied is: "<< j << std::endl;
			//std::cout<<"The molar fraction of the species "<< Species_name <<" is " << interpolated_point <<" at abscissa "<<Abscissa_temp<< std::endl;
			break;
		}				
	}
	return interpolated_point;
}

/*std::vector<double> PlugFlowReactor_Plugin::Solve_Species(std::vector<std::string> Species)
{
	Mole_frac_temp_vec.resize(Species.size());

	OpenSMOKE::OpenSMOKEVectorDouble omega_Final(thermodynamicsMapXML_->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble x_Final(thermodynamicsMapXML_->NumberOfSpecies());
	
	if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_ISOTHERMAL)
	{
		plugflow_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	else if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_NONISOTHERMAL)
	{
		plugflow_non_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	MW_Final = thermodynamicsMapXML_->MolecularWeight_From_MassFractions(omega_Final.GetHandle());
	thermodynamicsMapXML_->MoleFractions_From_MassFractions(x_Final.GetHandle(), MW_Final, omega_Final.GetHandle());
			
	for (int i=0; i< Species.size(); i++)
	{
		Mole_frac_temp_vec[i] = x_Final(thermodynamicsMapXML_->IndexOfSpecies(Species[i]));
    	} 
return Mole_frac_temp_vec;

}*/
void PlugFlowReactor_Plugin::clean_up()
{
	delete plugflow_options_ ;
	delete ode_parameters_ ;	
	delete sensitivity_options_ ;
	delete onTheFlyROPA_ ;
	delete on_the_fly_post_processing_ ;
	delete polimi_soot_;
	delete idt;
	if (temperature_profile)
	{
		delete profile;
	}

	if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_ISOTHERMAL)
	{
		delete plugflow_isothermal_ ;	
	}
	else if (type_ == OpenSMOKE::PLUGFLOW_REACTOR_NONISOTHERMAL)
	{
		delete plugflow_non_isothermal_ ;	
	}
}
}
