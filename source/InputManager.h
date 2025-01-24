#pragma once

namespace OptiSMOKE {
class InputManager {
 public:
  InputManager(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary);

  ~InputManager();

  void SetInputOptions(int argc, char* argv[]);

  void ReadDictionary();

  void DakotaInputString();

  // ==================================================
  // Standard Map
  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;
  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;
  OpenSMOKE::TransportPropertiesMap_CHEMKIN* transportMapXML_;

  // ==================================================
  // Nominal Map
  OpenSMOKE::ThermodynamicsMap_CHEMKIN* nominalthermodynamicsMapXML_;
  OpenSMOKE::KineticsMap_CHEMKIN* nominalkineticsMapXML_;
  OpenSMOKE::TransportPropertiesMap_CHEMKIN* nominaltransportMapXML_;

 private:
  OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_;

  // ==================================================
  // Variables of main dictionaries
  std::string input_file_name_;
  std::string main_dictionary_;
  fs::path output_folder_;
  fs::path kinetics_folder_;
  fs::path optimized_kinetics_folder_;

  // ==================================================
  // Grammar Allocation
  OptiSMOKE::grammar_optismoke main_grammar_;

  // ==================================================
  // Options
  OptiSMOKE::options_kinetics kinetics_data_;
  OptiSMOKE::options_optimization_target optimization_target_;
  OptiSMOKE::options_optimization_setup optimization_setup_;
  OptiSMOKE::options_curvematching curvematching_options_;
  OptiSMOKE::options_dakota dakota_options_;

  // ==================================================
  // Reading data from json files
  // OptiSMOKE::DataManager data_manager_;

  std::string optimization_library_;
  std::vector<std::string> path_experimental_data_files_;
  std::string dakota_input_string_;
};
}  // namespace OptiSMOKE

#include "InputManager.hpp"
