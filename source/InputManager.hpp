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
|					  Timoteo Dinelli <timoteo.dinelli@polimi.it>	      |
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

namespace OptiSMOKE{

    InputManager::InputManager(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary) : dictionary_(dictionary){
        input_file_name_ = "input.dic";
        main_dictionary_ = "OptiSMOKEpp";
        output_folder_ = "Output";
        kinetics_folder_ = "kinetics";
        optimized_kinetics_folder_ = "Optimized_kinetics";

        iDebug_ = false;
        iDebugSimulations_ = false;
        iXml_ = false;
        iNominalXml_ = false;
        iTransport_ = false;
    }
    
    InputManager::~InputManager(){}

    void InputManager::SetInputOptions(int argc, char* argv[]){
        //Input Options
        {
            po::options_description desc("Allowed options");
            desc.add_options()
            ("help", "Help Message")
            ("input", po::value<std::string>(), "Input File Path (default: \"input.dic\")");

            po::variables_map vm;
            po::store(po::parse_command_line(argc, argv, desc), vm);
            po::notify(vm);

            if (vm.count("help")){
                std::cout << desc << std::endl;
            }

            if (vm.count("input")){
                input_file_name_ = vm["input"].as<std::string>();
            }
        }
    }

    void InputManager::ReadDictionary(){
        // It can be trivial however this is for future 
        // parallelization with MPI see:
        // https://github.com/astagni/DoctorSMOKEpp/blob/main/src/DataManager.hpp
        // Remember that this all goes under rank=0
        ReadMainDictionary();
        
        // This to process or not kinetics folder
        if(!iXml_ || !iNominalXml_){

            if(!iTransport_){
                OpenSMOKE::RapidKineticMechanismWithoutTransport(
                    kinetics_data_.chemkin_output(),
                    kinetics_data_.chemkin_thermodynamics(),
                    kinetics_data_.chemkin_kinetics());
            }
            else{
                OpenSMOKE::RapidKineticMechanismWithTransport(
                    kinetics_data_.chemkin_output(),
                    kinetics_data_.chemkin_transport(),
                    kinetics_data_.chemkin_thermodynamics(),
                    kinetics_data_.chemkin_kinetics());
            }
        }

        CreateMaps();
    }

    void InputManager::ReadMainDictionary(){
        dictionary_.ReadDictionariesFromFile(input_file_name_);
        dictionary_(main_dictionary_).SetGrammar(main_grammar_);

        // kinetics folder
        if(dictionary_(main_dictionary_).CheckOption("@KineticsFolder")){
            iXml_ = true;
            dictionary_(main_dictionary_).ReadPath("@KineticsFolder", kinetics_folder_);
            if(!fs::exists(kinetics_folder_)){
                OptiSMOKE::FatalErrorMessage("The @KineticsFolder path does not exists!");
            }
            OpenSMOKE::CheckKineticsFolder(kinetics_folder_);
        
            // nominal kinetic pre-processor
            dictionary_(main_dictionary_).ReadDictionary("@NominalKineticsPreProcessor", preprocessor_dictionary_);
            kinetics_data_.SetupFromDictionary(dictionary_, preprocessor_dictionary_, iTransport_);
        }
        else if(dictionary_(main_dictionary_).CheckOption("@NominalKineticsFolder")){
            iNominalXml_ = true;
            dictionary_(main_dictionary_).ReadPath("@NominalKineticsFolder", kinetics_folder_);
            if(!fs::exists(kinetics_folder_)){
                OptiSMOKE::FatalErrorMessage("The @NominalKineticsFolder path does not exists!");
            }
            OpenSMOKE::CheckKineticsFolder(kinetics_folder_);
            
            // kinetic pre-processor
            dictionary_(main_dictionary_).ReadDictionary("@KineticsPreProcessor", preprocessor_dictionary_);
            kinetics_data_.SetupFromDictionary(dictionary_, preprocessor_dictionary_, iTransport_);
        }
        else{
            OptiSMOKE::FatalErrorMessage("Please provide a folder the kinetic mechanism available are: @NominalKineticsFolder | @KineticsFolder");
        }

        // name of optimized kinetic folder
        dictionary_(main_dictionary_).ReadPath("@NameOfOptimizedKineticsFolder", optimized_kinetics_folder_);
       
        // debug
        if(dictionary_(main_dictionary_).CheckOption("@Debug")){
            dictionary_(main_dictionary_).ReadBool("@Debug", iDebug_);
        }

        // debug simulations
        if(dictionary_(main_dictionary_).CheckOption("@DebugSim")){
            dictionary_(main_dictionary_).ReadBool("@DebugSim", iDebug_);
        }

        // path OpenSMOKE input files
        dictionary_(main_dictionary_).ReadPath("@PathDatasetInputFiles", path_opensmoke_input_files_);

        // path data set input files
        dictionary_(main_dictionary_).ReadPath("@PathExperimentalDataFiles", path_experimental_data_files_);

        // Dictionaries

        // Dakota options
        dictionary_(main_dictionary_).ReadDictionary("@DakotaOptions", dakota_dictionary_);
        dakota_options_.SetupFromDictionary(dictionary_, dakota_dictionary_);

        // CM options
        dictionary_(main_dictionary_).ReadDictionary("@CurveMatchingOptions", curvematching_dictionary_);
        curvematching_options_.SetupFromDictionary(dictionary_, curvematching_dictionary_);

        // Optimization setup
        dictionary_(main_dictionary_).ReadDictionary("@OptimizationSetup", optimization_setup_dictionary_);
        optimization_setup_.SetupFromDictionary(dictionary_, optimization_setup_dictionary_);

        // Optimization target
        dictionary_(main_dictionary_).ReadDictionary("@OptimizationTarget", optimization_target_dictionary_);
        optimization_target_.SetupFromDictionary(dictionary_, optimization_target_dictionary_);
        
    }

    void InputManager::DakotaInputString(){
		dakota_input_string_ = "\tenvironment,"
		                       "\n\t\ttabular_data";
		dakota_input_string_.append("\n\t\t\ttabular_data_file '" + dakota_options_.tabular_data_file() + "'");
		
        dakota_input_string_.append("\n	method,"); 
		dakota_input_string_.append("\n 		" + dakota_options_.method());
		dakota_input_string_.append("\n		  max_iterations = " + dakota_options_.max_iterations());
		dakota_input_string_.append("\n 		  max_function_evaluations = " + dakota_options_.max_function_evaluations());
        dakota_input_string_.append("\n		  convergence_tolerance = " + dakota_options_.convergence_tolerance());
		dakota_input_string_.append("\n 		  solution_target = " + dakota_options_.solution_target());
		dakota_input_string_.append("\n 		  seed = " + dakota_options_.seed());

		if(dakota_options_.iDiverseInput()){
		    dakota_input_string_.append("\n");
			for (int i = 0; i < dakota_options_.diverse_dakota_input().size(); i++)
			{
				dakota_input_string_.append( " " + dakota_options_.diverse_dakota_input()[i]);
			}
		} 
        else if(dakota_options_.method() == "coliny_ea"){
		  	dakota_input_string_.append("\n 		  population_size = " + dakota_options_.population_size());
			dakota_input_string_.append("\n		  fitness_type " + dakota_options_.fitness_type());
			dakota_input_string_.append("\n		  mutation_type " + dakota_options_.mutation_type());
			dakota_input_string_.append("\n		  mutation_rate " + dakota_options_.mutation_rate());
			dakota_input_string_.append("\n		  crossover_type " + dakota_options_.crossover_type());
			dakota_input_string_.append("\n		  crossover_rate " + dakota_options_.crossover_rate());
			dakota_input_string_.append("\n		  replacement_type " + dakota_options_.replacement_type());
		} 
        else if(dakota_options_.method() == "coliny_direct"){
			dakota_input_string_.append("\n                 division " + dakota_options_.division());
			dakota_input_string_.append("\n                 max_boxsize_limit " + dakota_options_.max_boxsize_limit());
			dakota_input_string_.append("\n                 min_boxsize_limit " + dakota_options_.min_boxsize_limit());
		}
        else{
            OptiSMOKE::FatalErrorMessage("Available implemented methods are coliny_ea | coliny_direct");
        }
			
        dakota_input_string_.append("\n	variables,");

		if(optimization_setup_.parameter_distribution() == "uniform"){
			dakota_input_string_.append("\n		continuous_design = " + std::to_string(optimization_target_.number_of_parameters()));
			dakota_input_string_.append("\n		  descriptors " + param_name_string);
			dakota_input_string_.append("\n 		  initial_point " + initial_values_string);
			dakota_input_string_.append("\n 	 	  lower_bounds " + lower_bounds_string);
			dakota_input_string_.append("\n 		  upper_bounds " + upper_bounds_string);
		}
        else if (optimization_setup_.parameter_distribution() == "normal"){
			dakota_input_string_.append("\n               active uncertain " );
			dakota_input_string_.append("\n		normal_uncertain = " + std::to_string(optimization_target_.number_of_parameters()));
			dakota_input_string_.append("\n		  descriptors " + param_name_string);
			dakota_input_string_.append("\n 		  means " + initial_values_string);
			dakota_input_string_.append("\n 		  std_deviations " + std_deviations_string);
		}// Da fare check su consistenza nell' input sul tipo di parameter boundary

		dakota_input_string_.append("\n	interface,");
		dakota_input_string_.append("\n		direct");
		dakota_input_string_.append("\n		  analysis_driver = 'plugin_opensmoke'");
		dakota_input_string_.append("\n	responses,");
		dakota_input_string_.append("\n		num_objective_functions = 1");
			
		// Options to use other optimization method (e.g. gradient-based)
        // Qua forse va messa la possibilità di fare altri tipi di gradienti 
        // accordingly to dakota sicuro lo faccio ora non c'ho voglia
		if (dakota_options_.iGradient() == true){
			dakota_input_string_.append("\n               numerical_gradients");
			dakota_input_string_.append("\n               	method_source dakota");
			dakota_input_string_.append("\n               	interval_type forward");
			dakota_input_string_.append("\n               	fd_step_size = 1.e-5");
		}
		else{
			dakota_input_string_.append("\n		no_gradients");
		}

		dakota_input_string_.append("\n		no_hessians");
    }

    void InputManager::CreateMaps(){
        // This goes under kinetics map
        fs::path path_kinetics_output;
	
        if (!iXml_) // To be interpreted on-the-fly
            path_kinetics_output = kinetics_data_.chemkin_output();
        else if (iXml_) // Already in XML format
            path_kinetics_output = kinetics_folder_;

        std::cout.setstate(std::ios_base::failbit); // Disable video output
        boost::property_tree::ptree ptree;
    	boost::property_tree::read_xml( (path_kinetics_output / "kinetics.xml").string(), ptree );

        thermodynamicsMapXML_ = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
        kineticsMapXML_ = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML_, ptree);
        if(iTransport_)
            transportMapXML_ = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree);
        std::cout.clear(); // Re-enable video output

        // This goes under nominal kinetics map
        fs::path path_nominal_kinetics_output;
	
        if (!iNominalXml_) // To be interpreted on-the-fly
            path_nominal_kinetics_output = kinetics_data_.chemkin_output();
        else if (iNominalXml_) // Already in XML format
            path_nominal_kinetics_output = kinetics_folder_;

        std::cout.setstate(std::ios_base::failbit); // Disable video output
        boost::property_tree::ptree nominal_ptree;
    	boost::property_tree::read_xml( (path_nominal_kinetics_output / "kinetics.xml").string(), nominal_ptree);

        nominalthermodynamicsMapXML_ = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(nominal_ptree);
        nominalkineticsMapXML_ = new OpenSMOKE::KineticsMap_CHEMKIN(*nominalthermodynamicsMapXML_, nominal_ptree);
        if(iTransport_)
            nominaltransportMapXML_ = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(nominal_ptree);
        std::cout.clear(); // Re-enable video output

    }

    void InputManager::InitialParameters(){

    }     

    void InputManager::ComputeBoundaries(){

		double T_low = 300;
		double T_high = 2500;

		if(optimization_setup_.parameter_boundaries() == "Furst"){
		    // Initialize needed values at the specific size
		    std::vector<double> list_of_nominal_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> list_of_nominal_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> list_of_nominal_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());

		    std::vector<double> list_of_min_abs_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> list_of_max_abs_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> list_of_min_abs_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> list_of_max_abs_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> list_of_min_abs_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> list_of_max_abs_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());

		    std::vector<double> kappa_lower_T_low(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> kappa_upper_T_low(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> kappa_lower_T_high(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> kappa_upper_T_high(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> Beta_1(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> Beta_2(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> lnA_1(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> lnA_2(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> E_over_R_1(optimization_target_.list_of_target_uncertainty_factors().size());
		    std::vector<double> E_over_R_2(optimization_target_.list_of_target_uncertainty_factors().size());

		    for (unsigned int i = 0; i < optimization_target_.list_of_target_uncertainty_factors().size(); i++){
                // Nominal values of parameters
                list_of_nominal_lnA_double[i] = std::log(nominalkineticsMapXML_->A(optimization_target_.list_of_target_uncertainty_factors()[i]-1));
                list_of_nominal_Beta_double[i] = nominalkineticsMapXML_->Beta(optimization_target_.list_of_target_uncertainty_factors()[i]-1);
                list_of_nominal_E_over_R_double[i] = nominalkineticsMapXML_->E_over_R(optimization_target_.list_of_target_uncertainty_factors()[i]-1);
                
                // Min and Max of lnA
                list_of_min_abs_lnA_double[i] = list_of_nominal_lnA_double[i] + std::log(std::pow(10, -optimization_target_.list_of_target_uncertainty_factors()[i]));
                list_of_max_abs_lnA_double[i] = list_of_nominal_lnA_double[i] + std::log(std::pow(10, optimization_target_.list_of_target_uncertainty_factors()[i]));
                if (std::find(optimization_target_.list_of_target_lnA().begin(),
                            optimization_target_.list_of_target_lnA().end(),
                            optimization_target_.list_of_target_uncertainty_factors()[i]) != optimization_target_.list_of_target_lnA().end()){
                    list_of_min_abs_lnA_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_double[i]));
                    list_of_max_abs_lnA_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_double[i]));
                }
                
			    // Limiting values for the rate coefficient
			    kappa_lower_T_low[i] = list_of_min_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_low) - list_of_nominal_E_over_R_double[i]*std::pow(T_low,-1);
			    kappa_upper_T_low[i] = list_of_max_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_low) - list_of_nominal_E_over_R_double[i]*std::pow(T_low,-1);
			    kappa_lower_T_high[i] = list_of_min_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_high) - list_of_nominal_E_over_R_double[i]*std::pow(T_high,-1);
			    kappa_upper_T_high[i] = list_of_max_abs_lnA_double[i] + list_of_nominal_Beta_double[i]*std::log(T_high) - list_of_nominal_E_over_R_double[i]*std::pow(T_high,-1);
				
                // Calculating extreme values for Beta
				Beta_1[i] = (kappa_upper_T_low[i] - kappa_lower_T_high[i] - list_of_nominal_E_over_R_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));
				Beta_2[i] = (kappa_lower_T_low[i] - kappa_upper_T_high[i] - list_of_nominal_E_over_R_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));

				list_of_min_abs_Beta_double[i] = std::min(Beta_1[i],Beta_2[i]);
				list_of_max_abs_Beta_double[i] = std::max(Beta_1[i],Beta_2[i]);
				
                if (std::find(optimization_target_.list_of_target_Beta().begin(),
                              optimization_target_.list_of_target_Beta().end(),
                              optimization_target_.list_of_target_uncertainty_factors()[i]) != optimization_target_.list_of_target_Beta().end()){
					list_of_min_abs_Beta_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_double[i]));
					list_of_max_abs_Beta_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_double[i]));
				}

				// Calculting extreame values of E_over_R
				lnA_1[i] = ( kappa_lower_T_high[i] - (T_low/T_high) * kappa_upper_T_low[i] - list_of_nominal_Beta_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
				E_over_R_1[i] = lnA_1[i] * T_low + T_low * list_of_nominal_Beta_double[i] * std::log(T_low) - kappa_upper_T_low[i] * T_low;
                lnA_2[i] = ( kappa_upper_T_high[i] - (T_low/T_high) * kappa_lower_T_low[i] - list_of_nominal_Beta_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
				E_over_R_2[i] = lnA_2[i] * T_low + T_low * list_of_nominal_Beta_double[i] * std::log(T_low) - kappa_lower_T_low[i] * T_low;    
                list_of_min_abs_E_over_R_double[i] = std::min(E_over_R_1[i],E_over_R_2[i]);
				list_of_max_abs_E_over_R_double[i] = std::max(E_over_R_1[i],E_over_R_2[i]);
				
                if (std::find(optimization_target_.list_of_target_E_over_R().begin(),
                              optimization_target_.list_of_target_E_over_R().end(),
                              optimization_target_.list_of_target_uncertainty_factors()[i]) != optimization_target_.list_of_target_E_over_R().end()){
					list_of_min_abs_E_over_R_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_double[i]));
					list_of_max_abs_E_over_R_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_double[i]));
				}
		    }

		    // Initialize needed values
		    std::vector<double> list_of_nominal_lnA_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> list_of_nominal_Beta_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> list_of_nominal_E_over_R_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());

		    std::vector<double> list_of_min_abs_lnA_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> list_of_max_abs_lnA_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> list_of_min_abs_Beta_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> list_of_max_abs_Beta_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> list_of_min_abs_E_over_R_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> list_of_max_abs_E_over_R_inf_double(optimization_target_.list_of_target_uncertainty_factors_inf().size());

		    std::vector<double> kappa_lower_T_low_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> kappa_upper_T_low_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> kappa_lower_T_high_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> kappa_upper_T_high_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> Beta_1_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> Beta_2_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> lnA_1_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> lnA_2_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> E_over_R_1_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());
		    std::vector<double> E_over_R_2_inf(optimization_target_.list_of_target_uncertainty_factors_inf().size());

            std::vector<unsigned int> indices_of_falloff_reactions = nominalkineticsMapXML_->IndicesOfFalloffReactions();

            for (unsigned int i=0; i < optimization_target_.list_of_target_uncertainty_factors_inf().size(); i++){
                
                int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(), indices_of_falloff_reactions.end(), optimization_target_.list_of_target_uncertainty_factors_inf()[i])-indices_of_falloff_reactions.begin();
                
                // Nominal values of inf parameters
                list_of_nominal_lnA_inf_double[i] = std::log(nominalkineticsMapXML_->A_falloff_inf(pos_FallOff_Reaction));
                list_of_nominal_Beta_inf_double[i] = nominalkineticsMapXML_->Beta_falloff_inf(pos_FallOff_Reaction);
                list_of_nominal_E_over_R_inf_double[i] = nominalkineticsMapXML_->E_over_R_falloff_inf(pos_FallOff_Reaction);

                // Min and Max of lnA_inf
                list_of_min_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,-optimization_target_.list_of_target_uncertainty_factors_inf()[i]));
                list_of_max_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,optimization_target_.list_of_target_uncertainty_factors_inf()[i]));
                if (std::find(optimization_target_.list_of_target_lnA_inf().begin(), 
                            optimization_target_.list_of_target_lnA_inf().end(),
                            optimization_target_.list_of_target_uncertainty_factors_inf()[i]) != optimization_target_.list_of_target_lnA_inf().end()){
                    list_of_min_abs_lnA_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_inf_double[i]));
                    list_of_max_abs_lnA_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_inf_double[i]));
                }
                
                // Limiting values for the rate coefficient
                kappa_lower_T_low_inf[i]  = list_of_min_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_low) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_low,-1);
                kappa_upper_T_low_inf[i]  = list_of_max_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_low) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_low,-1);
                kappa_lower_T_high_inf[i] = list_of_min_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_high) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_high,-1);
                kappa_upper_T_high_inf[i] = list_of_max_abs_lnA_inf_double[i] + list_of_nominal_Beta_inf_double[i]*std::log(T_high) - list_of_nominal_E_over_R_inf_double[i]*std::pow(T_high,-1);

                Beta_1_inf[i] = (kappa_upper_T_low_inf[i] - kappa_lower_T_high_inf[i] - list_of_nominal_E_over_R_inf_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));
                Beta_2_inf[i] = (kappa_lower_T_low_inf[i] - kappa_upper_T_high_inf[i] - list_of_nominal_E_over_R_inf_double[i] * (1/T_high - 1/T_low))/(std::log(T_low) - std::log(T_high));
                list_of_min_abs_Beta_inf_double[i] = std::min(Beta_1_inf[i],Beta_2_inf[i]);
                list_of_max_abs_Beta_inf_double[i] = std::max(Beta_1_inf[i],Beta_2_inf[i]);
                if (std::find(optimization_target_.list_of_target_Beta_inf().begin(),
                            optimization_target_.list_of_target_Beta_inf().end(),
                            optimization_target_.list_of_target_uncertainty_factors_inf()[i]) != optimization_target_.list_of_target_Beta_inf().end()){
                    list_of_min_abs_Beta_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_inf_double[i]));
                    list_of_max_abs_Beta_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_inf_double[i]));
                }

                // Calculting extreame values of E_over_R
                lnA_1_inf[i] = ( kappa_lower_T_high_inf[i] - (T_low/T_high) * kappa_upper_T_low_inf[i] - list_of_nominal_Beta_inf_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
                E_over_R_1_inf[i] = lnA_1_inf[i] * T_low + T_low * list_of_nominal_Beta_inf_double[i] * std::log(T_low) - kappa_upper_T_low_inf[i] * T_low;
                lnA_2_inf[i] = ( kappa_upper_T_high_inf[i] - (T_low/T_high) * kappa_lower_T_low_inf[i] - list_of_nominal_Beta_inf_double[i] * (std::log(T_high) - (T_low/T_high) * std::log(T_low)) ) / (1 - (T_low/T_high));
                E_over_R_2_inf[i] = lnA_2_inf[i] * T_low + T_low * list_of_nominal_Beta_inf_double[i] * std::log(T_low) - kappa_lower_T_low_inf[i] * T_low;
                list_of_min_abs_E_over_R_inf_double[i] = std::min(E_over_R_1_inf[i],E_over_R_2_inf[i]);
                list_of_max_abs_E_over_R_inf_double[i] = std::max(E_over_R_1_inf[i],E_over_R_2_inf[i]);
                if (std::find(optimization_target_.list_of_target_E_over_R_inf().begin(),
                            optimization_target_.list_of_target_E_over_R_inf().end(),
                            optimization_target_.list_of_target_uncertainty_factors_inf()[i]) != optimization_target_.list_of_target_E_over_R_inf().end()){
                    list_of_min_abs_E_over_R_inf_.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_inf_double[i]));
                    list_of_max_abs_E_over_R_inf_.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_inf_double[i]));
                }
            }
		}
		
        if(optimization_setup_.parameter_boundaries() == "Narrow")
		{
			std::vector<double> list_of_nominal_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> list_of_nominal_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> list_of_nominal_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());

			std::vector<double> list_of_min_abs_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> list_of_max_abs_lnA_double(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> list_of_min_abs_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> list_of_max_abs_Beta_double(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> list_of_min_abs_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> list_of_max_abs_E_over_R_double(optimization_target_.list_of_target_uncertainty_factors().size());
			
			std::vector<double> Beta_1(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> Beta_2(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> lnA_1(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> lnA_2(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> E_over_R_1(optimization_target_.list_of_target_uncertainty_factors().size());
			std::vector<double> E_over_R_2(optimization_target_.list_of_target_uncertainty_factors().size());
			
			for (unsigned int i=0; i < optimization_target_.list_of_target_uncertainty_factors().size(); i++){
				list_of_nominal_lnA_double[i] = std::log(nominalkineticsMapXML_->A(optimization_target_.list_of_target_uncertainty_factors()[i]-1));
				list_of_nominal_Beta_double[i] = nominalkineticsMapXML_->Beta(optimization_target_.list_of_target_uncertainty_factors()[i]-1);
				list_of_nominal_E_over_R_double[i] = nominalkineticsMapXML_->E_over_R(optimization_target_.list_of_target_uncertainty_factors()[i]-1);
				

				list_of_min_abs_lnA_double[i] = list_of_nominal_lnA_double[i]+std::log(std::pow(10, -optimization_target_.list_of_target_uncertainty_factors()[i]));
				list_of_max_abs_lnA_double[i] = list_of_nominal_lnA_double[i]+std::log(std::pow(10, optimization_target_.list_of_target_uncertainty_factors()[i]));
				
				if (std::find(list_of_target_lnA.begin(),list_of_target_lnA.end(),list_of_target_uncertainty_factors[i]) != list_of_target_lnA.end()){	
					list_of_min_abs_lnA.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_double[i]));
					list_of_max_abs_lnA.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_double[i]));
				}

				Beta_1[i] = list_of_nominal_Beta_double[i]+std::log(std::pow(10,list_of_uncertainty_factors[i])) / std::log(T_high);
				Beta_2[i] = list_of_nominal_Beta_double[i]-std::log(std::pow(10,list_of_uncertainty_factors[i])) / std::log(T_high);
				
				list_of_min_abs_Beta_double[i] = std::min(Beta_1[i],Beta_2[i]);
				list_of_max_abs_Beta_double[i] = std::max(Beta_1[i],Beta_2[i]);

				if (std::find(list_of_target_Beta.begin(),list_of_target_Beta.end(),list_of_target_uncertainty_factors[i]) != list_of_target_Beta.end()){
					list_of_min_abs_Beta.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_double[i]));
					list_of_max_abs_Beta.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_double[i]));
				}

				E_over_R_1[i] = list_of_nominal_E_over_R_double[i]-std::log(std::pow(10,list_of_uncertainty_factors[i])) * T_low;
				E_over_R_2[i] = list_of_nominal_E_over_R_double[i]+std::log(std::pow(10,list_of_uncertainty_factors[i])) * T_low;
				
				list_of_min_abs_E_over_R_double[i] = std::min(E_over_R_1[i],E_over_R_2[i]);
				list_of_max_abs_E_over_R_double[i] = std::max(E_over_R_1[i],E_over_R_2[i]);
				
				if (std::find(list_of_target_E_over_R.begin(),list_of_target_E_over_R.end(),list_of_target_uncertainty_factors[i]) != list_of_target_E_over_R.end()){
					list_of_min_abs_E_over_R.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_double[i]));
					list_of_max_abs_E_over_R.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_double[i]));
				}
			}

            std::vector<double> list_of_nominal_lnA_inf_double;
            std::vector<double> list_of_nominal_Beta_inf_double;
            std::vector<double> list_of_nominal_E_over_R_inf_double;

            std::vector<double> list_of_min_abs_lnA_inf_double;
            std::vector<double> list_of_max_abs_lnA_inf_double;
            std::vector<double> list_of_min_abs_Beta_inf_double;
            std::vector<double> list_of_max_abs_Beta_inf_double;
            std::vector<double> list_of_min_abs_E_over_R_inf_double;
            std::vector<double> list_of_max_abs_E_over_R_inf_double;

            std::vector<double> Beta_1_inf;
            std::vector<double> Beta_2_inf;
            std::vector<double> lnA_1_inf;
            std::vector<double> lnA_2_inf;
            std::vector<double> E_over_R_1_inf;
            std::vector<double> E_over_R_2_inf;

            // Resize them
            list_of_nominal_lnA_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
            list_of_nominal_Beta_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
            list_of_nominal_E_over_R_inf_double.resize(list_of_target_uncertainty_factors_inf.size());

            list_of_min_abs_lnA_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
            list_of_max_abs_lnA_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
            list_of_min_abs_Beta_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
            list_of_max_abs_Beta_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
            list_of_min_abs_E_over_R_inf_double.resize(list_of_target_uncertainty_factors_inf.size());
            list_of_max_abs_E_over_R_inf_double.resize(list_of_target_uncertainty_factors_inf.size());

            Beta_1_inf.resize(list_of_target_uncertainty_factors_inf.size());
            Beta_2_inf.resize(list_of_target_uncertainty_factors_inf.size());
            lnA_1_inf.resize(list_of_target_uncertainty_factors_inf.size());
            lnA_2_inf.resize(list_of_target_uncertainty_factors_inf.size());
            E_over_R_1_inf.resize(list_of_target_uncertainty_factors_inf.size());
            E_over_R_2_inf.resize(list_of_target_uncertainty_factors_inf.size());

		for (int i=0; i < list_of_target_uncertainty_factors_inf.size(); i++)
		{
			int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(),indices_of_falloff_reactions.end(),list_of_target_uncertainty_factors_inf[i])-indices_of_falloff_reactions.begin();
			// Nominal values of inf parameters
			list_of_nominal_lnA_inf_double[i] = std::log(nominalkineticsMapXML->A_falloff_inf(pos_FallOff_Reaction));
			list_of_nominal_Beta_inf_double[i] = nominalkineticsMapXML->Beta_falloff_inf(pos_FallOff_Reaction);
			list_of_nominal_E_over_R_inf_double[i] = nominalkineticsMapXML->E_over_R_falloff_inf(pos_FallOff_Reaction);

			// Min and Max of lnA_inf
			list_of_min_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,-list_of_uncertainty_factors_inf[i]));
			list_of_max_abs_lnA_inf_double[i] = list_of_nominal_lnA_inf_double[i]+std::log(std::pow(10,list_of_uncertainty_factors_inf[i]));
			if (std::find(list_of_target_lnA_inf.begin(),list_of_target_lnA_inf.end(),list_of_target_uncertainty_factors_inf[i]) != list_of_target_lnA_inf.end())
			{
				list_of_min_abs_lnA_inf.push_back(boost::lexical_cast<std::string>(list_of_min_abs_lnA_inf_double[i]));
				list_of_max_abs_lnA_inf.push_back(boost::lexical_cast<std::string>(list_of_max_abs_lnA_inf_double[i]));
			}
			
			// Calculating extreme values for Beta
			Beta_1_inf[i] = list_of_nominal_Beta_inf_double[i]+std::log(std::pow(10,list_of_target_uncertainty_factors_inf[i])) / std::log(T_high);
			Beta_2_inf[i] = list_of_nominal_Beta_inf_double[i]-std::log(std::pow(10,list_of_target_uncertainty_factors_inf[i])) / std::log(T_high);;
			list_of_min_abs_Beta_inf_double[i] = std::min(Beta_1_inf[i],Beta_2_inf[i]);
			list_of_max_abs_Beta_inf_double[i] = std::max(Beta_1_inf[i],Beta_2_inf[i]);
			if (std::find(list_of_target_Beta_inf.begin(),list_of_target_Beta_inf.end(),list_of_target_uncertainty_factors_inf[i]) != list_of_target_Beta_inf.end())
			{
					list_of_min_abs_Beta_inf.push_back(boost::lexical_cast<std::string>(list_of_min_abs_Beta_inf_double[i]));
					list_of_max_abs_Beta_inf.push_back(boost::lexical_cast<std::string>(list_of_max_abs_Beta_inf_double[i]));
			}

			// Calculting extreame values of E_over_R
			E_over_R_1_inf[i] = list_of_nominal_E_over_R_inf_double[i]-std::log(std::pow(10,list_of_target_uncertainty_factors_inf[i])) * T_low;
			E_over_R_2_inf[i] = list_of_nominal_E_over_R_inf_double[i]+std::log(std::pow(10,list_of_target_uncertainty_factors_inf[i])) * T_low;
			list_of_min_abs_E_over_R_inf_double[i] = std::min(E_over_R_1_inf[i],E_over_R_2_inf[i]);
			list_of_max_abs_E_over_R_inf_double[i] = std::max(E_over_R_1_inf[i],E_over_R_2_inf[i]);
			if (std::find(list_of_target_E_over_R_inf.begin(),list_of_target_E_over_R_inf.end(),list_of_target_uncertainty_factors_inf[i]) != list_of_target_E_over_R_inf.end())
			{
				list_of_min_abs_E_over_R_inf.push_back(boost::lexical_cast<std::string>(list_of_min_abs_E_over_R_inf_double[i]));
				list_of_max_abs_E_over_R_inf.push_back(boost::lexical_cast<std::string>(list_of_max_abs_E_over_R_inf_double[i]));
			}
	
		}

		}

		if(optimization_setup_.parameter_boundaries() == "Re-parametrization"){
            OptiSMOKE::ErrorMessage("Compute Boundaries", "Re-Implementation not done yet");
        }

		// EPLR - Alpha, Beta, Eps 
		for (int i=0; i < list_of_target_EPLR.size(); i++ )
		{						
			list_of_nominal_lnA_EPLR.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_lnA_EPLR.push_back(boost::lexical_cast<std::string>(-list_of_uncertainty_factors_EPLR[i]));
			list_of_max_lnA_EPLR.push_back(boost::lexical_cast<std::string>( list_of_uncertainty_factors_EPLR[i]));
		}

		for (int i=0; i < list_of_target_EPLR.size(); i++ )
		{					
			list_of_nominal_ER_EPLR.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_ER_EPLR.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_EPLR[i]))*T_low));
			list_of_max_ER_EPLR.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_EPLR[i]))*T_low));
		}

		for (int i=0; i < list_of_target_EPLR.size(); i++ )
		{
			list_of_nominal_Beta_EPLR.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_Beta_EPLR.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_EPLR[i])) / std::log(T_high)));
			list_of_max_Beta_EPLR.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_EPLR[i])) / std::log(T_high)));
		}
		
		// EPLR - Alpha, Beta, Eps 

		// Extended PLOG - Alpha, Beta, Eps 
		for (int i=0; i < list_of_target_extplog.size(); i++ )
		{						
			list_of_nominal_lnA_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_lnA_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(-list_of_uncertainty_factors_extplog[i]));
			list_of_max_lnA_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>( list_of_uncertainty_factors_extplog[i]));
		}

		for (int i=0; i < list_of_target_extplog.size(); i++ )
		{					
			list_of_nominal_ER_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_ER_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_extplog[i]))*T_low));
			list_of_max_ER_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_extplog[i]))*T_low));
		}

		for (int i=0; i < list_of_target_extplog.size(); i++ )
		{
			list_of_nominal_Beta_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_Beta_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_extplog[i])) / std::log(T_high)));
			list_of_max_Beta_ext_plog_coefficients.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_extplog[i])) / std::log(T_high)));
		}

		// Extended PLOGS - THIRD BODIES!
		int pos_extended_plog_reaction;
		int pos_extended_plog_species;
		for (int i=0; i < list_of_target_extended_plog_reactions.size(); i++ )
		{
			// find the position of the reaction having index list_of_target_extended_plog_reactions[i] within indices_of_extendedplogs
			pos_extended_plog_reaction = std::find(indices_of_extendedplogs_opt.begin(),indices_of_extendedplogs_opt.end(),list_of_target_extended_plog_reactions[i])-indices_of_extendedplogs_opt.begin();
			for (int k=0; k < nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).species().size(); k++ )
			{		
				std::string temp_string = nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).species()[k];
				if(!nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).species()[k].compare(0,temp_string.length(),list_of_target_extended_plog_species[i]))
			    {
               		pos_extended_plog_species = k;
               	        break;
			        
               	}
			}
			// get the value of third body			
			list_of_nominal_TB_ExtPLOG.push_back(boost::lexical_cast<std::string>(nominalkineticsMapXML->extendedplogopt_reactions(pos_extended_plog_reaction).ThirdBody(pos_extended_plog_species)[0]));
			list_of_min_TB_ExtPLOG.push_back(boost::lexical_cast<std::string>(list_of_min_tb_extplog[i]));
			list_of_max_TB_ExtPLOG.push_back(boost::lexical_cast<std::string>(list_of_max_tb_extplog[i]));
		}


 		// CLASSIC PLOG - Alpha, Beta, Eps
		for (int i=0; i < list_of_target_classic_plog_reactions.size(); i++ )
		{						
			list_of_nominal_lnA_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_lnA_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(-list_of_uncertainty_factors_classic_plog[i]));
			list_of_max_lnA_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>( list_of_uncertainty_factors_classic_plog[i]));
		}

		for (int i=0; i < list_of_target_classic_plog_reactions.size(); i++ )
		{			
			// the nominal random variable is 0, so that Eps_0 = Esp_0 + D is verified:		
			list_of_nominal_ER_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			// The minimum and the maximum values of the random variable are then computed as follows:
			list_of_min_ER_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_classic_plog[i]))*T_low));
			list_of_max_ER_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_classic_plog[i]))*T_low));
		}

		for (int i=0; i < list_of_target_classic_plog_reactions.size(); i++ )
		{
			list_of_nominal_Beta_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(0));
			list_of_min_Beta_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(-std::log(std::pow(10,list_of_uncertainty_factors_classic_plog[i])) / std::log(T_high)));
			list_of_max_Beta_classic_plog_coefficients.push_back(boost::lexical_cast<std::string>(+std::log(std::pow(10,list_of_uncertainty_factors_classic_plog[i])) / std::log(T_high)));
		}
	}

    }
} // namespace OptiSMOKE


