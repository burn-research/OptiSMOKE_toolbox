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

    InputManager::InputManager(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary) : dictionary_(dictionary)
    {
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

    void InputManager::SetInputOptions(int argc, char* argv[])
    {
        //Input Options
        {
            po::options_description desc("Allowed options");
            desc.add_options()
            ("help", "Help Message")
            ("input", po::value<std::string>(), "Input File Path (default: \"input.dic\")");

            po::variables_map vm;
            po::store(po::parse_command_line(argc, argv, desc), vm);
            po::notify(vm);

            if (vm.count("help"))
            {
                std::cout << desc << std::endl;
            }

            if (vm.count("input"))
            {
                input_file_name_ = vm["input"].as<std::string>();
            }
        }
    }

    void InputManager::ReadDictionary()
    {
        // It can be trivial however this is for future 
        // parallelization with MPI see:
        // https://github.com/astagni/DoctorSMOKEpp/blob/main/src/DataManager.hpp
        // Remember that this all goes under rank=0
        ReadMainDictionary();
        
        // This to process or not kinetics folder
        if(!iXml_ || !iNominalXml_)
        {

            if(!iTransport_){
                OpenSMOKE::RapidKineticMechanismWithoutTransport(
                    kinetics_data_.chemkin_output(),
                    kinetics_data_.chemkin_thermodynamics(),
                    kinetics_data_.chemkin_kinetics());
            }
            else
            {
                OpenSMOKE::RapidKineticMechanismWithTransport(
                    kinetics_data_.chemkin_output(),
                    kinetics_data_.chemkin_transport(),
                    kinetics_data_.chemkin_thermodynamics(),
                    kinetics_data_.chemkin_kinetics());
            }
        }

        CreateMaps();
    }

    void InputManager::ReadMainDictionary()
    {
        dictionary_.ReadDictionariesFromFile(input_file_name_);
        dictionary_(main_dictionary_).SetGrammar(main_grammar_);

        // kinetics folder
        if(dictionary_(main_dictionary_).CheckOption("@KineticsFolder"))
        {
            iXml_ = true;
            dictionary_(main_dictionary_).ReadPath("@KineticsFolder", kinetics_folder_);
            if(!fs::exists(kinetics_folder_))
            {
                OptiSMOKE::FatalErrorMessage("The @KineticsFolder path does not exists!");
            }
            OpenSMOKE::CheckKineticsFolder(kinetics_folder_);
        
            // nominal kinetic pre-processor
            dictionary_(main_dictionary_).ReadDictionary("@NominalKineticsPreProcessor", preprocessor_dictionary_);
            kinetics_data_.SetupFromDictionary(dictionary_, preprocessor_dictionary_, iTransport_);
        }
        else if(dictionary_(main_dictionary_).CheckOption("@NominalKineticsFolder"))
        {
            iNominalXml_ = true;
            dictionary_(main_dictionary_).ReadPath("@NominalKineticsFolder", kinetics_folder_);
            if(!fs::exists(kinetics_folder_))
            {
                OptiSMOKE::FatalErrorMessage("The @NominalKineticsFolder path does not exists!");
            }
            OpenSMOKE::CheckKineticsFolder(kinetics_folder_);
            
            // kinetic pre-processor
            dictionary_(main_dictionary_).ReadDictionary("@KineticsPreProcessor", preprocessor_dictionary_);
            kinetics_data_.SetupFromDictionary(dictionary_, preprocessor_dictionary_, iTransport_);
        }
        else
        {
            OptiSMOKE::FatalErrorMessage("Please provide a folder the kinetic mechanism available are: @NominalKineticsFolder | @KineticsFolder");
        }

        // name of optimized kinetic folder
        dictionary_(main_dictionary_).ReadPath("@NameOfOptimizedKineticsFolder", optimized_kinetics_folder_);
       
        // debug
        if(dictionary_(main_dictionary_).CheckOption("@Debug"))
        {
            dictionary_(main_dictionary_).ReadBool("@Debug", iDebug_);
        }

        // debug simulations
        if(dictionary_(main_dictionary_).CheckOption("@DebugSim"))
        {
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

    void InputManager::DakotaInputString()
    {
		dakota_input_string_ = "     environment,"
		                       "\n      tabular_data";
		dakota_input_string_.append("\n 		tabular_data_file '" + dakota_options_.tabular_data_file() + "'");
		
        dakota_input_string_.append("\n	method,"); 
		dakota_input_string_.append("\n 		" + dakota_options_.method());
		dakota_input_string_.append("\n		  max_iterations = " + dakota_options_.max_iterations());
		dakota_input_string_.append("\n 		  max_function_evaluations = " + dakota_options_.max_function_evaluations());
        dakota_input_string_.append("\n		  convergence_tolerance = " + dakota_options_.convergence_tolerance());
		dakota_input_string_.append("\n 		  solution_target = " + dakota_options_.solution_target());
		dakota_input_string_.append("\n 		  seed = " + dakota_options_.seed());

		if(dakota_options_.iDiverseInput())
		{
		    dakota_input_string_.append("\n");
			for (int i = 0; i < dakota_options_.diverse_dakota_input().size(); i++)
			{
				dakota_input_string_.append( " " + dakota_options_.diverse_dakota_input()[i]);
			}
		} 
        else if(dakota_options_.method() == "coliny_ea")
	   	{
		  	dakota_input_string_.append("\n 		  population_size = " + dakota_options_.population_size());
			dakota_input_string_.append("\n		  fitness_type " + dakota_options_.fitness_type());
			dakota_input_string_.append("\n		  mutation_type " + dakota_options_.mutation_type());
			dakota_input_string_.append("\n		  mutation_rate " + dakota_options_.mutation_rate());
			dakota_input_string_.append("\n		  crossover_type " + dakota_options_.crossover_type());
			dakota_input_string_.append("\n		  crossover_rate " + dakota_options_.crossover_rate());
			dakota_input_string_.append("\n		  replacement_type " + dakota_options_.replacement_type());
		} 
        else if(dakota_options_.method() == "coliny_direct")
		{
			dakota_input_string_.append("\n                 division " + dakota_options_.division());
			dakota_input_string_.append("\n                 max_boxsize_limit " + dakota_options_.max_boxsize_limit());
			dakota_input_string_.append("\n                 min_boxsize_limit " + dakota_options_.min_boxsize_limit());
		}
        else
        {
            OptiSMOKE::FatalErrorMessage("Available implemented methods are coliny_ea | coliny_direct");
        }
			
        dakota_input_string_.append("\n	variables,");
        /*
		if(optimization_setup_.parameter_distribution() == "uniform")
        {
			dakota_input_string_.append("\n		continuous_design = " + std::to_string(optimization_target_.number_of_parameters()));
			dakota_input_string_.append("\n		  descriptors " + param_name_string);
			dakota_input_string_.append("\n 		  initial_point " + initial_values_string);
			dakota_input_string_.append("\n 	 	  lower_bounds " + lower_bounds_string);
			dakota_input_string_.append("\n 		  upper_bounds " + upper_bounds_string);
		}
        else if (optimization_setup_.parameter_distribution() == "normal")
        {
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
        // Qua forse va messa la possibilità di fare altri tipi di gradienti accordingly to dakota
		if (dakota_options_.iGradient() == true)
        {
			dakota_input_string_.append("\n               numerical_gradients");
			dakota_input_string_.append("\n               	method_source dakota");
			dakota_input_string_.append("\n               	interval_type forward");
			dakota_input_string_.append("\n               	fd_step_size = 1.e-5");
		}
		else
        {
			dakota_input_string_.append("\n		no_gradients");
		}

		dakota_input_string_.append("\n		no_hessians");*/
    }

    void InputManager::CreateMaps()
    {
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

    void InputManager::InitialParameters()
    {

    }     

    void InputManager::ComputeBoundaries()
    {

    }
} // namespace OptiSMOKE


