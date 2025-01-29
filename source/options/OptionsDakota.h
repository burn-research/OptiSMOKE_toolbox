/* ------------------------------------------------------------------------------- *\
|                                                                                   |
|             ____        __  _ _____ __  _______  __ __ ______                     |
|            / __ \____  / /_(_) ___//  |/  / __ \/ //_// ____/___  ____            |
|           / / / / __ \/ __/ /\__ \/ /|_/ / / / / ,<  / __/ / __ \/ __ \           |
|          / /_/ / /_/ / /_/ /___/ / /  / / /_/ / /| |/ /___/ /_/ / /_/ /           |
|          \____/ .___/\__/_//____/_/  /_/\____/_/ |_/_____/ .___/ .___/            |
|              /_/                                        /_/   /_/                 |
|                                                                                   |
| --------------------------------------------------------------------------------- |
|  Please refer to the copyright statement and license                              |
|  information at the end of this file.                                             |
| --------------------------------------------------------------------------------- |
|                                                                                   |
|        Authors: Timoteo Dinelli  <timoteo.dinelli@polimi.it>                      |
|                 Andrea Bertolino <andrea.bertolino@ulb.be>                        |
|                 Magnus Fürst     <magnus.furst@ulb.ac.be>                         |
|                                                                                   |
|          [1] CRECK Modeling Lab <https://www.creckmodeling.polimi.it>             |
|              Department of Chemistry, Materials and Chemical Engineering          |
|              Politecnico di Milano, P.zza Leonardo da Vinci 32, 20133 Milano      |
|                                                                                   |
|          [2] BRITE Research Group <https://brite-research.be>                     |
|              Brussels Institute for Thermal-fluid systems and clean Energy        |
|              Avenue F.D. Rooseveltlaan 50, Bruxelles 1050 Brussel                 |
|                                                                                   |
\* ------------------------------------------------------------------------------- */
#pragma once

namespace OptiSMOKE {
class OptionsDakota {
 public:
  OptionsDakota();

  ~OptionsDakota() {};

  void SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager, std::string dictionary_name);

  const std::string& method() const { return method_; };

  const std::string& population_size() const { return population_size_; };

  const std::string& fitness_type() const { return fitness_type_; };

  const std::string& mutation_type() const { return mutation_type_; };

  const std::string& mutation_rate() const { return mutation_rate_; };

  const std::string& crossover_type() const { return crossover_type_; };

  const std::string& crossover_rate() const { return crossover_rate_; };

  const std::string& replacement_type() const { return replacement_type_; };

  const std::string& max_iterations() const { return max_iterations_; };

  const std::string& max_function_evaluations() const { return max_function_evaluations_; };

  const std::string& convergence_tolerance() const { return convergence_tolerance_; };

  const std::string& solution_target() const { return solution_target_; };

  const std::string& seed() const { return seed_; };

  const std::vector<std::string>& diverse_dakota_input() const { return diverse_dakota_input_; };

  const std::string& division() const { return division_; };

  const std::string& max_boxsize_limit() const { return max_boxsize_limit_; };

  const std::string& min_boxsize_limit() const { return min_boxsize_limit_; };

  const bool& dakota_gradient() const { return dakota_gradient_; };

  const bool& diverse_input() const { return diverse_input_; };

  const std::string& tabular_data_file() const { return tabular_data_file_; };

  const bool& echo_dakota_string() const { return echo_dakota_string_; };

 private:
  GrammarDakota grammar_;

  std::vector<std::string> diverse_dakota_input_;
  bool diverse_input_;
  bool echo_dakota_string_;
  std::string tabular_data_file_;
  std::string method_;
  bool dakota_gradient_;

  // ==================================================
  // Values for coliny_evolutionary
  std::string population_size_;
  std::string fitness_type_;
  std::string mutation_type_;
  std::string mutation_rate_;
  std::string crossover_type_;
  std::string crossover_rate_;
  std::string replacement_type_;
  std::string max_iterations_;
  std::string max_function_evaluations_;
  std::string convergence_tolerance_;
  std::string solution_target_;
  std::string seed_;

  // ==================================================
  // Values for coliny_direct
  std::string division_;
  std::string max_boxsize_limit_;
  std::string min_boxsize_limit_;
};
}  // namespace OptiSMOKE

#include "OptionsDakota.hpp"
/* ------------------------------------------------------------------------------- *\
|                                                                                   |
|   MIT License                                                                     |
|                                                                                   |
|   Copyright (c) 2025 Timoteo Dinelli, Andrea Bertolino, Magnus Fürst              |
|                                                                                   |
|   Permission is hereby granted, free of charge, to any person obtaining a copy    |
|   of this software and associated documentation files (the "Software"), to deal   |
|   in the Software without restriction, including without limitation the rights    |
|   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       |
|   copies of the Software, and to permit persons to whom the Software is           |
|   furnished to do so, subject to the following conditions:                        |
|                                                                                   |
|   The above copyright notice and this permission notice shall be included in all  |
|   copies or substantial portions of the Software.                                 |
|                                                                                   |
|   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      |
|   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        |
|   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     |
|   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          |
|   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   |
|   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   |
|   SOFTWARE.                                                                       |
|                                                                                   |
\* ------------------------------------------------------------------------------- */
