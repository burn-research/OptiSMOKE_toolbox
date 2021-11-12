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
|            Author: Magnus Fürst <magnus.furst@ulb.ac.be>                |
|-------------------------------------------------------------------------|
|   License                                                               |
|                                                                         |
|   This file is part of OptiSMOKE.                                       |
|   Copyright (C) 2019 by Magnus Fürst                                    |
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
void PerfectlyStirredReactor_Plugin::Setup(const std::string input_file_name,
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML,
		OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML)
{
		// Pointers
		thermodynamicsMapXML_ = thermodynamicsMapXML;
		kineticsMapXML_ = kineticsMapXML;

		// Defines the grammar rules
		OpenSMOKE::Grammar_PerfectlyStirredReactor grammar_psr;
		
		// Define the dictionaries
		std::string main_dictionary_name_ = "PerfectlyStirredReactor";
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(main_dictionary_name_).SetGrammar(grammar_psr);

		// Read PSR reactor type
		{
			std::string value;
			if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true)
			{
				dictionaries(main_dictionary_name_).ReadString("@Type", value);
				if (value == "Isothermal-ConstantPressure")	type_ = OpenSMOKE::PERFECTLYSTIRRED_REACTOR_ISOTHERMAL_CONSTANTP;
				else if (value == "NonIsothermal-ConstantPressure")	type_ = OpenSMOKE::PERFECTLYSTIRRED_REACTOR_NONISOTHERMAL_CONSTANTP;
				else OpenSMOKE::FatalErrorMessage("Unknown perfectly stirred reactor type: " + value);

			}
		}


		// Read inlet conditions
		{
			std::string name_of_inlet_gas_status_subdictionary;
			if (dictionaries(main_dictionary_name_).CheckOption("@InletStatus") == true)
				dictionaries(main_dictionary_name_).ReadDictionary("@InletStatus", name_of_inlet_gas_status_subdictionary);

			GetGasStatusFromDictionary(dictionaries(name_of_inlet_gas_status_subdictionary), *thermodynamicsMapXML_, T_Inlet, P_Pa_Inlet, omega_Inlet);
		}
		// Read initial conditions
		{
			std::string name_of_initial_gas_status_subdictionary;
			if (dictionaries(main_dictionary_name_).CheckOption("@InitialStatus") == true)
			{
				dictionaries(main_dictionary_name_).ReadDictionary("@InitialStatus", name_of_initial_gas_status_subdictionary);
				GetGasStatusFromDictionary(dictionaries(name_of_initial_gas_status_subdictionary), *thermodynamicsMapXML_, T_Initial, P_Pa_Initial, omega_Initial);
			}
			else
			{
				T_Initial = T_Inlet;
				P_Pa_Initial = P_Pa_Inlet;
				omega_Initial = omega_Inlet;
			}
		}


		// Read end time (for transient simulations)
		tEnd = 1.e8;
		{
			double value;
		 	std::string units;
		 	if (dictionaries(main_dictionary_name_).CheckOption("@EndTime") == true)
		 	{
		 		dictionaries(main_dictionary_name_).ReadMeasure("@EndTime", value, units);
		 		if (units == "s")		  tEnd = value;
		 		else if (units == "ms")   tEnd = value/1000.;
		 		else if (units == "min")  tEnd = value*60.;
		 		else if (units == "h")    tEnd = value*3600.;
		 		else OpenSMOKE::FatalErrorMessage("Unknown time units");
		 	}
		 }

		unsigned int count_constraints = 0;

		// Read residence time
		residence_time = -1.;
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@ResidenceTime") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@ResidenceTime", value, units);
				if (units == "s")		  residence_time = value;
				else if (units == "ms")   residence_time = value/1000.;
				else if (units == "min")  residence_time = value*60.;
				else if (units == "h")    residence_time = value*3600.;
				else OpenSMOKE::FatalErrorMessage("Unknown time units");
				count_constraints++;
			}
		}		

		// Read mass flow rate
		mass_flow_rate = -1.;
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@MassFlowRate") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@MassFlowRate", value, units);
				if (units == "kg/s")	  mass_flow_rate = value;
				else if (units == "g/s")  mass_flow_rate = value/1000.;
				else OpenSMOKE::FatalErrorMessage("Unknown time units");
				count_constraints++;
			}
		}

		// Read volume
		volume = -1.;
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@Volume") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@Volume", value, units);
				if (units == "m3")		  volume = value;
				else if (units == "dm3")  volume = value/1.e3;
				else if (units == "cm3")  volume = value/1.e6;
				else if (units == "mm3")  volume = value/1.e9;
				else if (units == "l")    volume = value/1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown volume units");
				count_constraints++;
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
				if (units == "m2")		  exchange_area = value;
				else if (units == "dm2")  exchange_area = value/1.e2;
				else if (units == "cm2")  exchange_area = value/1.e4;
				else if (units == "mm2")  exchange_area = value/1.e6;
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
				else if (units == "kcal/m2/K")		global_thermal_exchange_coefficient = value*4186.8;
				else if (units == "kcal/m2/C")		global_thermal_exchange_coefficient = value*4186.8;
				else OpenSMOKE::FatalErrorMessage("Unknown global thermal exchange coefficient units");
			}
		}

		// Environment temperature
		T_environment = T_Inlet;
		{
			double value;
			std::string units;
			if (dictionaries(main_dictionary_name_).CheckOption("@EnvironmentTemperature") == true)
			{
				dictionaries(main_dictionary_name_).ReadMeasure("@EnvironmentTemperature", value, units);
				if (units == "K")		T_environment = value;
				else if (units == "C")	T_environment = value+273.15;
				else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
			}
		}
	
		if (count_constraints !=2)
			OpenSMOKE::FatalErrorMessage("Any two of @ResidenceTime @Volume @MassFlowRate must be specified");
	
		// Options
		{
			psr_options_ = new OpenSMOKE::PerfectlyStirredReactor_Options();
			if (dictionaries(main_dictionary_name_).CheckOption("@Options") == true)
			{
				std::string name_of_options_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@Options", name_of_options_subdictionary);
				psr_options_->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
			}
			psr_options_->SetVerboseVideo(false);
			psr_options_->SetVerboseOutput(false);
			psr_options_->SetSensitivityAnalysis(false);
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

				//psr_options_->SetSensitivityAnalysis(true);
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
			on_the_fly_post_processing_ = new OpenSMOKE::OnTheFlyPostProcessing(*thermodynamicsMapXML_, *kineticsMapXML_, psr_options_->output_path());

			if (dictionaries(main_dictionary_name_).CheckOption("@OnTheFlyPostProcessing") == true)
			{
				std::string name_of_options_subdictionary;
				dictionaries(main_dictionary_name_).ReadDictionary("@OnTheFlyPostProcessing", name_of_options_subdictionary);
				on_the_fly_post_processing_->SetupFromDictionary(dictionaries(name_of_options_subdictionary));
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


		// Virtual chemistry
		//virtual_chemistry_ = new OpenSMOKE::VirtualChemistry(*thermodynamicsMapXML_, false);

		// ------------------------------------------------------------------------------------------- //
		//                              Non parametric analysis                                        //
		// ------------------------------------------------------------------------------------------- //
		// Solve the ODE system: NonIsothermal, Constant Volume
		/*if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_NONISOTHERMAL_CONSTANTP)
		{
			psr_non_isothermal_ = new OpenSMOKE::PerfectlyStirredReactor_NonIsothermal_ConstantPressure(*thermodynamicsMapXML_, *kineticsMapXML_,
				*ode_parameters_, *psr_options_, *onTheFlyROPA_, *on_the_fly_post_processing_, *polimi_soot_, T_Initial, P_Pa_Initial, omega_Initial, T_Inlet, P_Pa_Inlet, omega_Inlet,					  residence_time, volume, mass_flow_rate, global_thermal_exchange_coefficient, exchange_area, T_environment);
		}


		else if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_ISOTHERMAL_CONSTANTP)
		{
			psr_isothermal_ = new OpenSMOKE::PerfectlyStirredReactor_Isothermal_ConstantPressure(*thermodynamicsMapXML_, *kineticsMapXML_,
				*ode_parameters_, *psr_options_, *onTheFlyROPA_, *on_the_fly_post_processing_, *polimi_soot_,T_Initial, P_Pa_Initial, omega_Initial,T_Inlet, P_Pa_Inlet, omega_Inlet, 					  residence_time, volume, mass_flow_rate);
		}*/
}

void PerfectlyStirredReactor_Plugin::Update_and_Solve_PSR(OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML, OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_updated)
{
	if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_ISOTHERMAL_CONSTANTP)
	{
		psr_isothermal_ = new OpenSMOKE::PerfectlyStirredReactor_Isothermal_ConstantPressure(*thermodynamicsMapXML, *kineticsMapXML_updated,
			*ode_parameters_, *psr_options_, *onTheFlyROPA_, *on_the_fly_post_processing_, *polimi_soot_,T_Initial, P_Pa_Initial, omega_Initial,T_Inlet, P_Pa_Inlet, omega_Inlet, 					  residence_time, volume, mass_flow_rate);
		psr_isothermal_->Solve(tEnd);
		
	}
	else if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_NONISOTHERMAL_CONSTANTP)
	{
		psr_non_isothermal_ = new OpenSMOKE::PerfectlyStirredReactor_NonIsothermal_ConstantPressure(*thermodynamicsMapXML, *kineticsMapXML_updated,
			*ode_parameters_, *psr_options_, *onTheFlyROPA_, *on_the_fly_post_processing_, *polimi_soot_, T_Initial, P_Pa_Initial, omega_Initial, T_Inlet, P_Pa_Inlet, omega_Inlet,					  residence_time, volume, mass_flow_rate, global_thermal_exchange_coefficient, exchange_area, T_environment);
		psr_non_isothermal_->Solve(tEnd);
		
	}
	
}
double PerfectlyStirredReactor_Plugin::Solve_Species(std::string Species)
{
	OpenSMOKE::OpenSMOKEVectorDouble omega_Final(thermodynamicsMapXML_->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble x_Final(thermodynamicsMapXML_->NumberOfSpecies());
	
	if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_ISOTHERMAL_CONSTANTP)
	{
		psr_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	else if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_NONISOTHERMAL_CONSTANTP)
	{
		psr_non_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	MW_Final = thermodynamicsMapXML_->MolecularWeight_From_MassFractions(omega_Final.GetHandle());
	thermodynamicsMapXML_->MoleFractions_From_MassFractions(x_Final.GetHandle(), MW_Final, omega_Final.GetHandle());
		
	Mole_frac_temp = x_Final(thermodynamicsMapXML_->IndexOfSpecies(Species));

	clean_up();


return Mole_frac_temp;

}

std::vector<double> PerfectlyStirredReactor_Plugin::Solve_Species(std::vector<std::string> Species)
{
	Mole_frac_temp_vec.resize(Species.size());

	OpenSMOKE::OpenSMOKEVectorDouble omega_Final(thermodynamicsMapXML_->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble x_Final(thermodynamicsMapXML_->NumberOfSpecies());
	
	if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_ISOTHERMAL_CONSTANTP)
	{
		psr_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	else if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_NONISOTHERMAL_CONSTANTP)
	{
		psr_non_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	//cioa
	MW_Final = thermodynamicsMapXML_->MolecularWeight_From_MassFractions(omega_Final.GetHandle());
	thermodynamicsMapXML_->MoleFractions_From_MassFractions(x_Final.GetHandle(), MW_Final, omega_Final.GetHandle());
			
	for (int i=0; i< Species.size(); i++)
	{
		if (Species[i]=="Temp" || Species[i]=="DeltaTemp" ){
			Mole_frac_temp_vec[i] = T_Final;
		} else {
			Mole_frac_temp_vec[i] = x_Final(thermodynamicsMapXML_->IndexOfSpecies(Species[i]));
		} 
    }

	clean_up();
	return Mole_frac_temp_vec;
}

std::vector<double> PerfectlyStirredReactor_Plugin::Solve_Multipl_Species(std::vector<std::string> Species_vec)
{
	Mole_frac_temp_vec.resize(Species_vec.size());

	OpenSMOKE::OpenSMOKEVectorDouble omega_Final(thermodynamicsMapXML_->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble x_Final(thermodynamicsMapXML_->NumberOfSpecies());
	
	if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_ISOTHERMAL_CONSTANTP)
	{
		psr_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	else if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_NONISOTHERMAL_CONSTANTP)
	{
		psr_non_isothermal_->GetFinalStatus(T_Final, P_Pa_Final, omega_Final);	
	}
	MW_Final = thermodynamicsMapXML_->MolecularWeight_From_MassFractions(omega_Final.GetHandle());
	thermodynamicsMapXML_->MoleFractions_From_MassFractions(x_Final.GetHandle(), MW_Final, omega_Final.GetHandle());
			
	for (int i=0; i< Species_vec.size(); i++)
	{
		if (Species_vec[i]=="Temp" || Species_vec[i]=="DeltaTemp"){
			Mole_frac_temp_vec[i] = T_Final;
		} else {
			Mole_frac_temp_vec[i] = x_Final(thermodynamicsMapXML_->IndexOfSpecies(Species_vec[i]));
		} 
    }

	clean_up();
	return Mole_frac_temp_vec;
}

void PerfectlyStirredReactor_Plugin::clean_up()
{
	delete psr_options_ ;
	delete ode_parameters_ ;	
	delete sensitivity_options_ ;
	delete onTheFlyROPA_ ;
	delete on_the_fly_post_processing_ ;
	delete polimi_soot_;

	if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_ISOTHERMAL_CONSTANTP)
	{
		delete psr_isothermal_;
	}
	else if (type_ == OpenSMOKE::PERFECTLYSTIRRED_REACTOR_NONISOTHERMAL_CONSTANTP)
	{
		delete psr_non_isothermal_;	
	}
}


}
