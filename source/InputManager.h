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

        // Da privatizzare
        void DakotaInputPreliminaryOptions();

        //
        void ReadExperimentalDataFiles();

        // Public access variables
        inline const OptiSMOKE::options_kinetics& kinetics_data() const {return kinetics_data_;};
        inline const OptiSMOKE::options_optimization_target& optimization_target() const {return optimization_target_;};
        inline const OptiSMOKE::options_optimization_setup& optimization_setup() const {return optimization_setup_;};
        inline const OptiSMOKE::options_curvematching& curvematching_options() const {return curvematching_options_;};
        inline const OptiSMOKE::options_dakota& dakota_options() const {return dakota_options_;};
        
        inline const bool& iXml() const {return iXml_;};
        inline const bool& iNominalXml() const {return iNominalXml_;};
        inline const std::string& input_file_name() const {return input_file_name_;};
        inline const std::string& main_dictionary() const {return main_dictionary_;};
        inline const std::string& optimization_library() const {return optimization_library_;};
        inline const fs::path& output_folder() const {return output_folder_;};
        inline const fs::path& kinetics_folder() const {return kinetics_folder_;};
        inline const fs::path& optimized_kinetics_folder() const {return optimized_kinetics_folder_;};
        inline const std::vector<std::string>& path_opensmoke_input_files() const {return path_opensmoke_input_files_;};
        inline const std::vector<std::string>& path_experimental_data_files() const {return path_experimental_data_files_;};
        inline const bool& iDebug() const {return iDebug_;};
        inline const bool& iDebugSimulations() const {return iDebugSimulations_;};
        inline const std::string& dakota_input_string() const {return dakota_input_string_;};
        inline const bool& iTransport() const {return iTransport_;};

        // Da checcare
        inline OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML() {return thermodynamicsMapXML_;};
        inline OpenSMOKE::ThermodynamicsMap_CHEMKIN* nominalthermodynamicsMapXML() {return nominalthermodynamicsMapXML_;};
        inline OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML() {return kineticsMapXML_;};
        inline OpenSMOKE::KineticsMap_CHEMKIN* nominalkineticsMapXML() {return nominalkineticsMapXML_;};
        inline OpenSMOKE::TransportPropertiesMap_CHEMKIN* transportMapXML() {return transportMapXML_;};
        inline OpenSMOKE::TransportPropertiesMap_CHEMKIN* nominaltransportMapXML() {return nominaltransportMapXML_;};

        inline const std::vector<std::vector<std::string>>& input_paths() const {return input_paths_;};
        inline const std::vector<std::string>& QoI() const {return QoI_;};
		inline const std::vector<std::string>& QoI_target() const {return QoI_target_;};        
        inline const std::vector<std::vector<std::vector<double>>>& expdata_x() const {return expdata_x_;};
        inline const std::vector<std::vector<std::vector<double>>>& expdata_y() const {return expdata_y_;};
        inline const std::vector<std::vector<std::vector<double>>>& uncertainty() const {return uncertainty_;};
        inline const std::vector<std::string>& dataset_names() const {return dataset_names_;};
        inline const std::vector<std::string>& solver_name() const {return solver_name_;};
    private:

        OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_;

        // Standard Map
        OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;
        OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;
        OpenSMOKE::TransportPropertiesMap_CHEMKIN* transportMapXML_;

        // Nominal Map
        OpenSMOKE::ThermodynamicsMap_CHEMKIN* nominalthermodynamicsMapXML_;
        OpenSMOKE::KineticsMap_CHEMKIN* nominalkineticsMapXML_;
        OpenSMOKE::TransportPropertiesMap_CHEMKIN* nominaltransportMapXML_;

        // Grammar Allocation
        OptiSMOKE::grammar_optismoke main_grammar_;

        // Options
        OptiSMOKE::options_kinetics kinetics_data_;
        OptiSMOKE::options_optimization_target optimization_target_;
        OptiSMOKE::options_optimization_setup optimization_setup_;
        OptiSMOKE::options_curvematching curvematching_options_;
        OptiSMOKE::options_dakota dakota_options_;

        //Reading data from json files
        OptiSMOKE::DataManager data_manager_;
        
        // Dictionaries string
        std::string preprocessor_dictionary_;
        std::string optimization_target_dictionary_;
        std::string optimization_setup_dictionary_;
        std::string curvematching_dictionary_;
        std::string dakota_dictionary_;
        
        // Variables of main dictionaries
        std::string input_file_name_;
        std::string main_dictionary_;
        std::string optimization_library_;
        fs::path output_folder_;
        fs::path kinetics_folder_;
        fs::path optimized_kinetics_folder_;
        std::vector<std::string> path_opensmoke_input_files_;
        std::vector<std::string> path_experimental_data_files_;
        bool iDebug_;
        bool iDebugSimulations_;
        bool iXml_;
        bool iNominalXml_;
        bool iTransport_;
        std::string dakota_input_string_;

        // Functions
        void ReadMainDictionary();

        void CreateMaps();
        
        void FromTargetToInitialParameter();
        
        void ComputeBoundaries();

        // Boundaries variables
        std::vector<std::string> list_of_min_abs_lnA_;
	    std::vector<std::string> list_of_max_abs_lnA_;
        
        std::vector<std::string> list_of_min_abs_Beta_;
	    std::vector<std::string> list_of_max_abs_Beta_;

        std::vector<std::string> list_of_min_abs_E_over_R_;
	    std::vector<std::string> list_of_max_abs_E_over_R_;

        std::vector<std::string> list_of_min_abs_lnA_inf_;
	    std::vector<std::string> list_of_max_abs_lnA_inf_;

        std::vector<std::string> list_of_min_abs_Beta_inf_;
	    std::vector<std::string> list_of_max_abs_Beta_inf_;

        std::vector<std::string> list_of_min_abs_E_over_R_inf_;
	    std::vector<std::string> list_of_max_abs_E_over_R_inf_;

        // EPLR
    	std::vector<std::string> list_of_nominal_lnA_EPLR_;
	    std::vector<std::string> list_of_min_lnA_EPLR_;
	    std::vector<std::string> list_of_max_lnA_EPLR_;
	
	    std::vector<std::string> list_of_nominal_ER_EPLR_;
	    std::vector<std::string> list_of_min_ER_EPLR_;
	    std::vector<std::string> list_of_max_ER_EPLR_;
	 
	    std::vector<std::string> list_of_nominal_Beta_EPLR_;
	    std::vector<std::string> list_of_min_Beta_EPLR_;
	    std::vector<std::string> list_of_max_Beta_EPLR_;

        std::vector<std::string> list_of_nominal_lnA_ext_plog_coefficients_;
	    std::vector<std::string> list_of_min_lnA_ext_plog_coefficients_;
	    std::vector<std::string> list_of_max_lnA_ext_plog_coefficients_;
	
	    std::vector<std::string> list_of_nominal_ER_ext_plog_coefficients_;
	    std::vector<std::string> list_of_min_ER_ext_plog_coefficients_;
	    std::vector<std::string> list_of_max_ER_ext_plog_coefficients_;
	
	    std::vector<std::string> list_of_nominal_Beta_ext_plog_coefficients_;
	    std::vector<std::string> list_of_min_Beta_ext_plog_coefficients_;
	    std::vector<std::string> list_of_max_Beta_ext_plog_coefficients_;

        std::vector<std::string> list_of_nominal_lnA_classic_plog_coefficients_;
	    std::vector<std::string> list_of_min_lnA_classic_plog_coefficients_;
	    std::vector<std::string> list_of_max_lnA_classic_plog_coefficients_;

	    std::vector<std::string> list_of_nominal_ER_classic_plog_coefficients_;
	    std::vector<std::string> list_of_min_ER_classic_plog_coefficients_;
	    std::vector<std::string> list_of_max_ER_classic_plog_coefficients_;

    	std::vector<std::string> list_of_nominal_Beta_classic_plog_coefficients_;
	    std::vector<std::string> list_of_min_Beta_classic_plog_coefficients_;
	    std::vector<std::string> list_of_max_Beta_classic_plog_coefficients_;

        // String needed by dakota in order to generate the input string
        std::string param_name_string_;
		std::string initial_values_string_;
		std::string lower_bounds_string_;
		std::string upper_bounds_string_;
		std::string std_deviations_string_;

        // Initial Parameters string
        std::vector<std::string> list_of_initial_lnA_;
        std::vector<std::string> list_of_initial_Beta_;
        std::vector<std::string> list_of_initial_E_over_R;
        std::vector<std::string> list_of_initial_lnA_inf_;
        std::vector<std::string> list_of_initial_Beta_inf_;
        std::vector<std::string> list_of_initial_E_over_R_inf_;
        std::vector<std::string> list_of_initial_thirdbody_eff_;

        // Experimental data files variables
        std::vector<std::string> dataset_names_;
		std::vector<std::vector<std::string>> input_paths_;
		std::vector<std::string> solver_name_;
		std::vector<std::string> QoI_;
		std::vector<std::string> QoI_target_;
		std::vector<bool> multiple_input_;
        std::vector<bool> save_simulations_;

		std::vector<std::vector<std::vector<std::string>>> ordinates_label_;
		std::vector<std::vector<std::vector<std::string>>> abscissae_label_;
		std::vector<std::vector<std::vector<std::string>>> uncertainty_kind_;
		std::vector<std::vector<std::vector<double>>> expdata_x_;
		std::vector<std::vector<std::vector<double>>> expdata_y_;
		std::vector<std::vector<std::vector<double>>> uncertainty_;
    };
} // namespace OptiSMOKE

#include "InputManager.hpp"
#endif // INPUTMANAGER_H