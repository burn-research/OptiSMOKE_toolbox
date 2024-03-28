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
|					  Timoteo Dinelli <timoteo.dinelli@polimi.it> |
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

namespace OptiSMOKE {
class InputManager {
 public:
  InputManager(OpenSMOKE::OpenSMOKE_DictionaryManager &dictionary);

  ~InputManager();

  void SetInputOptions(int argc, char *argv[]);

  void ReadDictionary();

  void DakotaInputString();

  void ReadExperimentalDataFiles();

  void SetUpNLOPT();

  // Public access variables
  const OptiSMOKE::options_kinetics &kinetics_data() const { return kinetics_data_; };

  const OptiSMOKE::options_optimization_target &optimization_target() const { return optimization_target_; };

  const OptiSMOKE::options_optimization_setup &optimization_setup() const { return optimization_setup_; };

  const OptiSMOKE::options_curvematching &curvematching_options() const { return curvematching_options_; };

  const OptiSMOKE::options_dakota &dakota_options() const { return dakota_options_; };

  const OptiSMOKE::options_nlopt &nlopt_options() const { return nlopt_options_; };

  const bool &iXml() const { return iXml_; };

  const std::string &input_file_name() const { return input_file_name_; };

  const std::string &main_dictionary() const { return main_dictionary_; };

  const std::string &optimization_library() const { return optimization_library_; };

  const fs::path &output_folder() const { return output_folder_; };

  const fs::path &kinetics_folder() const { return kinetics_folder_; };

  const fs::path &optimized_kinetics_folder() const { return optimized_kinetics_folder_; };

  const std::vector<std::string> &path_experimental_data_files() const { return path_experimental_data_files_; };

  const std::string &dakota_input_string() const { return dakota_input_string_; };

  const bool &iTransport() const { return iTransport_; };

  // Da checcare
  const OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML() const { return thermodynamicsMapXML_; };

  const OpenSMOKE::ThermodynamicsMap_CHEMKIN *nominalthermodynamicsMapXML() const {
    return nominalthermodynamicsMapXML_;
  };

  const OpenSMOKE::KineticsMap_CHEMKIN *kineticsMapXML() const { return kineticsMapXML_; };

  const OpenSMOKE::KineticsMap_CHEMKIN *nominalkineticsMapXML() const { return nominalkineticsMapXML_; };

  const OpenSMOKE::TransportPropertiesMap_CHEMKIN *transportMapXML() const { return transportMapXML_; };

  const OpenSMOKE::TransportPropertiesMap_CHEMKIN *nominaltransportMapXML() const { return nominaltransportMapXML_; };

  const std::vector<std::vector<std::string>> &input_paths() const { return input_paths_; };

  const std::vector<std::string> &QoI() const { return QoI_; };

  const std::vector<std::string> &QoI_target() const { return QoI_target_; };

  const std::vector<std::string> &reactor_mode() const { return reactor_mode_; };

  const std::vector<std::vector<std::vector<double>>> &expdata_x() const { return expdata_x_; };

  const std::vector<std::vector<std::vector<double>>> &expdata_y() const { return expdata_y_; };

  const std::vector<std::vector<std::vector<double>>> &uncertainty() const { return uncertainty_; };

  const std::vector<std::string> &dataset_names() const { return dataset_names_; };

  const std::vector<std::string> &solver_name() const { return solver_name_; };

  const std::vector<std::vector<std::string>> &ordinates_label() const { return ordinates_label_; };

  const std::vector<std::string> &param_str() const { return param_str_; };

  const std::vector<double> &initial_values() const { return initial_values_; };

  const std::vector<double> &lb() const { return lb_; };

  const std::vector<double> &ub() const { return ub_; };

  const fs::path &parametric_file_name() const { return parametric_file_name_; };

  // Standard Map
  OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML_;
  OpenSMOKE::KineticsMap_CHEMKIN *kineticsMapXML_;
  OpenSMOKE::TransportPropertiesMap_CHEMKIN *transportMapXML_;

  // Nominal Map
  OpenSMOKE::ThermodynamicsMap_CHEMKIN *nominalthermodynamicsMapXML_;
  OpenSMOKE::KineticsMap_CHEMKIN *nominalkineticsMapXML_;
  OpenSMOKE::TransportPropertiesMap_CHEMKIN *nominaltransportMapXML_;

 private:
  OpenSMOKE::OpenSMOKE_DictionaryManager &dictionary_;

  // Grammar Allocation
  OptiSMOKE::grammar_optismoke main_grammar_;

  // Options
  OptiSMOKE::options_kinetics kinetics_data_;
  OptiSMOKE::options_optimization_target optimization_target_;
  OptiSMOKE::options_optimization_setup optimization_setup_;
  OptiSMOKE::options_curvematching curvematching_options_;
  OptiSMOKE::options_dakota dakota_options_;
  OptiSMOKE::options_nlopt nlopt_options_;

  // Reading data from json files
  OptiSMOKE::DataManager data_manager_;

  // Dictionaries string
  std::string preprocessor_dictionary_;
  std::string optimization_target_dictionary_;
  std::string optimization_setup_dictionary_;
  std::string curvematching_dictionary_;
  std::string dakota_dictionary_;
  std::string nlopt_dictionary_;

  // Variables of main dictionaries
  std::string input_file_name_;
  std::string main_dictionary_;
  std::string optimization_library_;
  fs::path output_folder_;
  fs::path kinetics_folder_;
  fs::path optimized_kinetics_folder_;
  std::vector<std::string> path_experimental_data_files_;
  bool iXml_;
  bool iNominalXml_;
  bool iTransport_;
  std::string dakota_input_string_;

  // Functions
  void ReadMainDictionary();

  void CreateMaps();

  void FromTargetToInitialParameter();

  void ComputeBoundaries();

  void TargetsPreliminaryOptions();

  // Initial Parameters string assigned in
  // FromTargetToInitialParameters()
  std::vector<std::string> list_of_initial_lnA_;
  std::vector<std::string> list_of_initial_Beta_;
  std::vector<std::string> list_of_initial_E_over_R;
  std::vector<std::string> list_of_initial_lnA_inf_;
  std::vector<std::string> list_of_initial_Beta_inf_;
  std::vector<std::string> list_of_initial_E_over_R_inf_;
  std::vector<std::string> list_of_initial_thirdbody_eff_;

  // String needed by dakota in order to generate
  // the input string
  std::string param_name_string_;
  std::string initial_values_string_;
  std::string lower_bounds_string_;
  std::string upper_bounds_string_;
  std::string std_deviations_string_;

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

  std::vector<std::string> list_of_nominal_lnA_classic_plog_coefficients_;
  std::vector<std::string> list_of_min_lnA_classic_plog_coefficients_;
  std::vector<std::string> list_of_max_lnA_classic_plog_coefficients_;

  std::vector<std::string> list_of_nominal_ER_classic_plog_coefficients_;
  std::vector<std::string> list_of_min_ER_classic_plog_coefficients_;
  std::vector<std::string> list_of_max_ER_classic_plog_coefficients_;

  std::vector<std::string> list_of_nominal_Beta_classic_plog_coefficients_;
  std::vector<std::string> list_of_min_Beta_classic_plog_coefficients_;
  std::vector<std::string> list_of_max_Beta_classic_plog_coefficients_;

  std::vector<std::string> list_of_nominal_lnA_rpbmr_coefficients_;
  std::vector<std::string> list_of_min_lnA_rpbmr_coefficients_;
  std::vector<std::string> list_of_max_lnA_rpbmr_coefficients_;

  std::vector<std::string> list_of_nominal_Beta_rpbmr_coefficients_;
  std::vector<std::string> list_of_min_Beta_rpbmr_coefficients_;
  std::vector<std::string> list_of_max_Beta_rpbmr_coefficients_;

  std::vector<std::string> list_of_nominal_E_over_R_rpbmr_coefficients_;
  std::vector<std::string> list_of_min_E_over_R_rpbmr_coefficients_;
  std::vector<std::string> list_of_max_E_over_R_rpbmr_coefficients_;

  // Names vector to pass euther from dakota to nlopt
  // I know this is not the fanciest way
  std::vector<std::string> name_vec_lnA;
  std::vector<std::string> name_vec_lnA_inf;
  std::vector<std::string> name_vec_Beta;
  std::vector<std::string> name_vec_Beta_inf;
  std::vector<std::string> name_vec_E_over_R;
  std::vector<std::string> name_vec_E_over_R_inf;
  std::vector<std::string> name_vec_thirdbody;
  std::vector<std::string> name_vec_lnA_classic_plog;
  std::vector<std::string> name_vec_ER_classic_plog;
  std::vector<std::string> name_vec_Beta_classic_plog;
  std::vector<std::string> name_vec_lnA_rpbmr;
  std::vector<std::string> name_vec_ER_rpbmr;
  std::vector<std::string> name_vec_Beta_rpbmr;

  // Double vector needed by nlopt
  std::vector<double> initial_values_;
  std::vector<double> lb_;
  std::vector<double> ub_;
  std::vector<std::string> param_str_;
  fs::path parametric_file_name_;

  // Experimental data files variables
  std::vector<std::string> dataset_names_;
  std::vector<std::vector<std::string>> input_paths_;
  std::vector<std::string> solver_name_;
  std::vector<std::string> QoI_;
  std::vector<std::string> reactor_mode_;
  std::vector<std::string> QoI_target_;
  std::vector<bool> multiple_input_;

  std::vector<std::vector<std::string>> ordinates_label_;
  std::vector<std::vector<std::string>> abscissae_label_;
  std::vector<std::vector<std::string>> uncertainty_kind_;
  std::vector<std::vector<std::vector<double>>> expdata_x_;
  std::vector<std::vector<std::vector<double>>> expdata_y_;
  std::vector<std::vector<std::vector<double>>> uncertainty_;
};
}  // namespace OptiSMOKE

#include "InputManager.hpp"
#endif  // INPUTMANAGER_H
