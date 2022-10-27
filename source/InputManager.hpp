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
        optimized_kinetics_folder_ = "Optimal_kinetics";

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
        ReadMainDictionary();
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
        
        // penalty function
        /* L'ho spostata in optimization setup
        if(dictionary_(main_dictionary_).CheckOption("@PenaltyFunction"))
        {
            dictionary_(main_dictionary_).ReadBool("@PenaltyFunction", iPenaltyFunction_);
        }
        */
       
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
        dictionary_(main_dictionary_).ReadDictionary("@OptimizationTargets", optimization_target_dictionary_);
        optimization_target_.SetupFromDictionary(dictionary_, optimization_target_dictionary_);
        
    }
} // namespace OptiSMOKE


