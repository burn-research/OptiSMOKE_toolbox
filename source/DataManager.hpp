/* ----------------------------------------------------------------------------------- *\
|                                                                                       |
|                ____        __  _ _____ __  _______  __ __ ______                      |
|               / __ \____  / /_(_) ___//  |/  / __ \/ //_// ____/___  ____             |
|              / / / / __ \/ __/ /\__ \/ /|_/ / / / / ,<  / __/ / __ \/ __ \            |
|             / /_/ / /_/ / /_/ /___/ / /  / / /_/ / /| |/ /___/ /_/ / /_/ /            |
|             \____/ .___/\__/_//____/_/  /_/\____/_/ |_/_____/ .___/ .___/             |
|                 /_/                                        /_/   /_/                  |
|                                                                                       |
| ------------------------------------------------------------------------------------- |
|  See license and copyright at the end of this file.                                   |
| ------------------------------------------------------------------------------------- |
|                                                                                       |
|           Authors: Timoteo Dinelli  <timoteo.dinelli@polimi.it>                       |
|                    Andrea Bertolino <andrea.bertolino@ulb.be>                         |
|                    Magnus Fürst     <magnus.furst@ulb.ac.be>                          |
|                                                                                       |
|           [1] CRECK Modeling Lab <https://www.creckmodeling.polimi.it>                |
|               Department of Chemistry, Materials and Chemical Engineering             |
|               Politecnico di Milano, P.zza Leonardo da Vinci 32, 20133 Milano         |
|                                                                                       |
|           [2] BRITE Research Group <https://brite-research.be>                        |
|               Brussels Institute for Thermal-fluid systems and clean Energy           |
|               Avenue F.D. Rooseveltlaan 50, Bruxelles 1050 Brussel                    |
|                                                                                       |
\* ----------------------------------------------------------------------------------- */

namespace OptiSMOKE {

DataManager::DataManager(const std::string file_path) : filename_(file_path) {}

bool DataManager::LoadFile() {
  try {
    boost::property_tree::read_json(filename_, root_);
    return true;
  } catch (const boost::property_tree::json_parser_error& e) {
    std::cerr << "Error reading JSON file (" << filename_ << "): " << e.what() << std::endl;
    return false;
  }
}

bool DataManager::ParseSimulationInformations() {
  try {
    auto& sim_tree = root_.get_child("simulation_info");
    sim_info_.solver = sim_tree.get<std::string>("solver");

    if (auto reactor_mode = sim_tree.get_optional<std::string>("reactor_mode")) {
      sim_info_.reactor_mode = reactor_mode.get();
    } else {
      sim_info_.reactor_mode = "None";
    }

    sim_info_.QoI = sim_tree.get<std::string>("QoI");

    if (auto QoI_target = sim_tree.get_optional<std::string>("QoI_target")) {
      sim_info_.QoI_target = QoI_target.get();
    } else {
      sim_info_.QoI_target = "None";
    }

    if (auto multiple_input = sim_tree.get_optional<bool>("multiple_input")) {
      sim_info_.multiple_input = multiple_input.get();
    } else {
      sim_info_.multiple_input = false;
    }

    if (auto save_simulations_data = sim_tree.get_optional<bool>("save_simulations_data")) {
      sim_info_.save_simulations_data = save_simulations_data.get();
    } else {
      sim_info_.save_simulations_data = false;
    }

    for (const auto& item : sim_tree.get_child("OS_Input_File")) {
      sim_info_.OS_Input_File.push_back(item.second.get_value<std::string>());
    }
    return ValidateSimulationKeywords();
  } catch (const boost::property_tree::ptree_error& e) {
    std::cerr << "Error parsing simulation info (" << filename_ << "): " << e.what() << std::endl;
    return false;
  }
}

bool DataManager::ParseExperimentalData() {
  try {
    exp_data_.clear();  // Clear any existing data

    for (const auto& data_entry : root_.get_child("data")) {
      ExperimentalDataset data_set;
      auto& data = data_entry.second;  // Get the data object

      // Parse basic information
      data_set.abscissae_label = data.get<std::string>("abscissae_label");
      data_set.abscissae_unit = data.get<std::string>("abscissae_unit");
      data_set.ordinates_label = data.get<std::string>("ordinates_label");
      data_set.ordinates_unit = data.get<std::string>("ordinates_unit");

      // Parse arrays
      for (const auto& item : data.get_child("abscissae")) {
        data_set.abscissae.push_back(item.second.get_value<double>());
      }
      for (const auto& item : data.get_child("ordinates")) {
        data_set.ordinates.push_back(item.second.get_value<double>());
      }

      exp_data_.push_back(data_set);
    }
    return true;
  } catch (const boost::property_tree::ptree_error& e) {
    std::cerr << "Parsing Error parsing experimental data (" << filename_ << "): " << e.what() << std::endl;
    return false;
  }
}

bool DataManager::ParseBasicInformations() {
  try {
    dataset_name_ = root_.get<std::string>("name");
    return true;
  } catch (const boost::property_tree::ptree_error& e) {
    std::cerr << "Parsing Error parsing basic info (" << filename_ << "): " << e.what() << std::endl;
    return false;
  }
}

void DataManager::PrintSimulationInformations() const {
  std::cout << "Simulation Info:\n";
  std::cout << " Solver: " << sim_info_.solver << "\n";
  std::cout << " Reactor Mode: " << sim_info_.reactor_mode << "\n";
  std::cout << " QoI: " << sim_info_.QoI << "\n";
  std::cout << " QoI Target: " << sim_info_.QoI_target << "\n";
  std::cout << " Multiple Input: " << (sim_info_.multiple_input ? "true" : "false") << "\n";
  std::cout << " Save Simulations Data: " << (sim_info_.save_simulations_data ? "true" : "false") << "\n";
  std::cout << " Input Files:\n";
  for (const auto& file : sim_info_.OS_Input_File) {
    std::cout << "  - " << file << "\n";
  }
}

void DataManager::PrintExperimentalData() const {
  std::cout << "Experimental Data Sets (" << exp_data_.size() << " sets):\n";
  for (size_t dataset_idx = 0; dataset_idx < exp_data_.size(); ++dataset_idx) {
    const auto& data_set = exp_data_[dataset_idx];
    std::cout << "\n Data Set " << dataset_idx + 1 << ":\n";
    std::cout << "  Abscissae Label: " << data_set.abscissae_label << "\n";
    std::cout << "  Abscissae Unit: " << data_set.abscissae_unit << "\n";
    std::cout << "  Ordinates Label: " << data_set.ordinates_label << "\n";
    std::cout << "  Ordinates Unit: " << data_set.ordinates_unit << "\n";
    std::cout << "  Data Points:\n";
    for (size_t i = 0; i < data_set.abscissae.size(); ++i) {
      std::cout << "   " << data_set.abscissae[i] << " " << data_set.ordinates[i] << "\n";
    }
  }
}

void DataManager::PrintDataSet(const size_t index) const {
  if (index >= exp_data_.size()) {
    std::cerr << "Error: Invalid data set index\n";
    return;
  }

  const auto& data_set = exp_data_[index];
  std::cout << "Data Set " << index + 1 << ":\n";
  std::cout << "Abscissae Label: " << data_set.abscissae_label << "\n";
  std::cout << "Abscissae Unit: " << data_set.abscissae_unit << "\n";
  std::cout << "Ordinates Label: " << data_set.ordinates_label << "\n";
  std::cout << "Ordinates Unit: " << data_set.ordinates_unit << "\n";
  std::cout << "Data Points:\n";
  for (size_t i = 0; i < data_set.abscissae.size(); ++i) {
    std::cout << "  " << data_set.abscissae[i] << " " << data_set.abscissae_unit << " -> " << data_set.ordinates[i]
              << " " << data_set.ordinates_unit << "\n";
  }
}

bool DataManager::ValidateSimulationKeywords() {
  // Validate solver
  if (std::find(valid_solvers_.begin(), valid_solvers_.end(), sim_info_.solver) == valid_solvers_.end()) {
    std::cerr << "Parsing Error (" << filename_ << "):\nInvalid solver type '" << sim_info_.solver
              << "'.\nValid options are:\n ";
    for (const auto& solver : valid_solvers_) {
      std::cerr << solver << " ";
    }
    std::cerr << std::endl;
    return false;
  }

  // Validate reactor mode
  if (std::find(valid_reactor_modes_.begin(), valid_reactor_modes_.end(), sim_info_.reactor_mode)
      == valid_reactor_modes_.end()) {
    std::cerr << "Parsing Error (" << filename_ << "):\nInvalid reactor mode '" << sim_info_.reactor_mode
              << "'.\nValid options are:\n";
    for (const auto& mode : valid_reactor_modes_) {
      std::cerr << mode << " ";
    }
    std::cerr << std::endl;
    return false;
  }

  // Validate QoI type
  if (std::find(valid_QoI_types_.begin(), valid_QoI_types_.end(), sim_info_.QoI) == valid_QoI_types_.end()) {
    std::cerr << "Parsing Error(" << filename_ << "):\nInvalid QoI type '" << sim_info_.QoI
              << "'.\nValid options are:\n";
    for (const auto& qoi : valid_QoI_types_) {
      std::cerr << qoi << " ";
    }
    std::cerr << std::endl;
    return false;
  }

  // Validate QoI targets
  if (std::find(valid_QoI_targets_.begin(), valid_QoI_targets_.end(), sim_info_.QoI_target)
      == valid_QoI_targets_.end()) {
    std::cerr << "Parsing Error (" << filename_ << "):\nInvalid QoI targets '" << sim_info_.QoI_target
              << "'.\nValid options are:\n";
    for (const auto& qoi_target : valid_QoI_targets_) {
      std::cerr << qoi_target << " ";
    }
    std::cerr << std::endl;
    return false;
  }

  return true;
}

bool DataManager::ValidateExperimentalData() {
  for (size_t i = 0; i < exp_data_.size(); ++i) {
    const auto& dataset = exp_data_[i];

    // Validate array sizes
    if (dataset.abscissae.size() != dataset.ordinates.size()) {
      std::cerr << "Parsing Error in dataset " << filename_ << ": Mismatched sizes between abscissae and ordinates"
                << std::endl;
      return false;
    }

    // Validate that abscissae are monotonically increasing
    for (size_t j = 1; j < dataset.abscissae.size(); ++j) {
      if (dataset.abscissae[j] <= dataset.abscissae[j - 1]) {
        std::cerr << "Parsing Error in dataset " << filename_ << ": Abscissae values must be strictly increasing"
                  << std::endl;
        return false;
      }
    }
  }

  return true;
}
}  // namespace OptiSMOKE
/* ----------------------------------------------------------------------------------- *\
|                                                                                       |
|     MIT License                                                                       |
|                                                                                       |
|     Copyright (c) 2025 Timoteo Dinelli, Andrea Bertolino, Magnus Fürst                |
|                                                                                       |
|     Permission is hereby granted, free of charge, to any person obtaining a copy      |
|     of this software and associated documentation files (the "Software"), to deal     |
|     in the Software without restriction, including without limitation the rights      |
|     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell         |
|     copies of the Software, and to permit persons to whom the Software is             |
|     furnished to do so, subject to the following conditions:                          |
|                                                                                       |
|     The above copyright notice and this permission notice shall be included in all    |
|     copies or substantial portions of the Software.                                   |
|                                                                                       |
|     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR        |
|     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,          |
|     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       |
|     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER            |
|     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,     |
|     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE     |
|     SOFTWARE.                                                                         |
|                                                                                       |
\* ----------------------------------------------------------------------------------- */
