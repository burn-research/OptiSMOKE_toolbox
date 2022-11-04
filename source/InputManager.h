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

#ifndef INPUTMANAGER_H
#define INPUTMANAGER_H

namespace OptiSMOKE{
    
    class InputManager
    {
    public:

        /// @brief Default Constructor
        InputManager(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary);

        /// @brief Default Destructor
        ~InputManager();

        /// @brief Function to handle the input file
        /// @param argc 
        /// @param argv 
        void SetInputOptions(int argc, char* argv[]);

        void ReadDictionary();
        // Da privatizzare
        void DakotaInputString();

        // Public access variables
        inline const OptiSMOKE::options_kinetics& kinetics_data() const {return kinetics_data_;};
        inline const bool& iXml() const {return iXml_;};
        inline const bool& iNominalXml() const {return iNominalXml_;};
        inline const std::string& input_file_name() const {return input_file_name_;};
        inline const std::string& main_dictionary() const {return main_dictionary_;};
        inline const fs::path& output_folder() const {return output_folder_;};
        inline const fs::path& kinetics_folder() const {return kinetics_folder_;};
        inline const fs::path& optimized_kinetics_folder() const {return optimized_kinetics_folder_;};
        inline const fs::path& path_opensmoke_input_files() const {return path_opensmoke_input_files_;};
        inline const fs::path& path_experimental_data_files() const {return path_experimental_data_files_;};
        inline const bool& iDebug() const {return iDebug_;};
        inline const bool& iDebugSimulations() const {return iDebugSimulations_;};
        inline const std::string& dakota_input_string() const {return dakota_input_string_;};

    private:

        OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_;

        // Grammar Allocation
        OptiSMOKE::grammar_optismoke main_grammar_;

        // Options
        OptiSMOKE::options_kinetics kinetics_data_;
        OptiSMOKE::options_optimization_target optimization_target_;
        OptiSMOKE::options_optimization_setup optimization_setup_;
        OptiSMOKE::options_curvematching curvematching_options_;
        OptiSMOKE::options_dakota dakota_options_;

        // Dictionaries string
        std::string preprocessor_dictionary_;
        std::string optimization_target_dictionary_;
        std::string optimization_setup_dictionary_;
        std::string curvematching_dictionary_;
        std::string dakota_dictionary_;
        
        // Variables of main dictionaries
        std::string input_file_name_;
        std::string main_dictionary_;
        fs::path output_folder_;
        fs::path kinetics_folder_;
        fs::path optimized_kinetics_folder_;
        fs::path path_opensmoke_input_files_;
        fs::path path_experimental_data_files_;
        bool iDebug_;
        bool iDebugSimulations_;
        bool iXml_; // Forse ne devo tenere solo uno tra iXml_ e iNominalXml_ sta roba da capire
        bool iNominalXml_;
        bool iTransport_;
        std::string dakota_input_string_;

        // Functions
        void ReadMainDictionary();

    };
} // namespace OptiSMOKE

#include "InputManager.hpp"
#endif // INPUTMANAGER_H