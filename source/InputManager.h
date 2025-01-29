/* ----------------------------------------------------------------------------------- *\
|                                                                                       |
|                 ____        __  _ _____ __  _______  __ __ ______                     |
|                / __ \____  / /_(_) ___//  |/  / __ \/ //_// ____/___  ____            |
|               / / / / __ \/ __/ /\__ \/ /|_/ / / / / ,<  / __/ / __ \/ __ \           |
|              / /_/ / /_/ / /_/ /___/ / /  / / /_/ / /| |/ /___/ /_/ / /_/ /           |
|              \____/ .___/\__/_//____/_/  /_/\____/_/ |_/_____/ .___/ .___/            |
|                  /_/                                        /_/   /_/                 |
|                                                                                       |
| ------------------------------------------------------------------------------------- |
|  See license and copyright at the end of this file.                                   |
| ------------------------------------------------------------------------------------- |
|                                                                                       |
|            Authors: Timoteo Dinelli  <timoteo.dinelli@polimi.it>                      |
|                     Andrea Bertolino <andrea.bertolino@ulb.be>                        |
|                     Magnus Fürst     <magnus.furst@ulb.ac.be>                         |
|                                                                                       |
|            [1] CRECK Modeling Lab <https://www.creckmodeling.polimi.it>               |
|                Department of Chemistry, Materials and Chemical Engineering            |
|                Politecnico di Milano, P.zza Leonardo da Vinci 32, 20133 Milano        |
|                                                                                       |
|            [2] BRITE Research Group <https://brite-research.be>                       |
|                Brussels Institute for Thermal-fluid systems and clean Energy          |
|                Avenue F.D. Rooseveltlaan 50, Bruxelles 1050 Brussel                   |
|                                                                                       |
\* ----------------------------------------------------------------------------------- */
#pragma once

namespace OptiSMOKE {
class InputManager {
 public:
  InputManager(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary);

  ~InputManager();

  void SetInputOptions(int argc, char* argv[]);

  void ReadDictionary();

  void ReadExperimentalData();

  void DakotaInputString();

  const std::string& optimization_library() const { return optimization_library_; };

 private:
  OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_;

  // ==================================================
  // Standard Map
  std::shared_ptr<OpenSMOKE::ThermodynamicsMap_CHEMKIN> tmd_map_;
  std::shared_ptr<OpenSMOKE::KineticsMap_CHEMKIN> kin_map_;
  std::shared_ptr<OpenSMOKE::TransportPropertiesMap_CHEMKIN> tran_map_;

  // ==================================================
  // Nominal Map
  std::shared_ptr<OpenSMOKE::ThermodynamicsMap_CHEMKIN> nominal_tmd_map_;
  std::shared_ptr<OpenSMOKE::KineticsMap_CHEMKIN> nominal_kin_map_;
  std::shared_ptr<OpenSMOKE::TransportPropertiesMap_CHEMKIN> nominal_tran_map_;

  // ==================================================
  // Variables of main dictionaries
  std::string input_file_name_;
  std::string main_dictionary_;
  fs::path output_folder_;
  fs::path kinetics_folder_;
  fs::path optimized_kinetics_folder_;

  // ==================================================
  // Grammar Allocation
  OptiSMOKE::GrammarOptismoke main_grammar_;

  // ==================================================
  // Reading data from json files
  // OptiSMOKE::DataManager data_manager_;

  std::string optimization_library_;
  std::vector<std::string> path_experimental_data_files_;
};
}  // namespace OptiSMOKE

#include "InputManager.hpp"
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
