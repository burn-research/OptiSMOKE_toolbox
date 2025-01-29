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

namespace OptiSMOKE {
OptionsDakota::OptionsDakota() {
  method_ = "coliny_ea";
  population_size_ = "50";
  fitness_type_ = "merit_function";
  mutation_type_ = "offset_normal";
  mutation_rate_ = "1.0";
  crossover_type_ = "two_point";
  crossover_rate_ = "0.0";
  replacement_type_ = "chc = 10";

  division_ = "major_dimension";
  max_boxsize_limit_ = "0.0";
  min_boxsize_limit_ = "1.0e-4";

  dakota_gradient_ = false;
  diverse_input_ = false;

  tabular_data_file_ = "tabular_data.dat";

  echo_dakota_string_ = false;
}

void OptionsDakota::SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                                        std::string dictionary_name) {
  dictionary_manager(dictionary_name).SetGrammar(grammar_);

  if (dictionary_manager(dictionary_name).CheckOption("@TabularDataFile")) {
    dictionary_manager(dictionary_name).ReadString("@TabularDataFile", tabular_data_file_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@Method")) {
    dictionary_manager(dictionary_name).ReadString("@Method", method_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@MaxIterations")) {
    dictionary_manager(dictionary_name).ReadString("@MaxIterations", max_iterations_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@MaxFunctionEvaluations")) {
    dictionary_manager(dictionary_name).ReadString("@MaxFunctionEvaluations", max_function_evaluations_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@ConvergenceTolerance")) {
    dictionary_manager(dictionary_name).ReadString("@ConvergenceTolerance", convergence_tolerance_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@SolutionTarget")) {
    dictionary_manager(dictionary_name).ReadString("@SolutionTarget", solution_target_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@Seed")) {
    dictionary_manager(dictionary_name).ReadString("@Seed", seed_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@PopulationSize")) {
    dictionary_manager(dictionary_name).ReadString("@PopulationSize", population_size_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@FitnessType")) {
    dictionary_manager(dictionary_name).ReadString("@FitnessType", fitness_type_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@MutationType")) {
    dictionary_manager(dictionary_name).ReadString("@MutationType", mutation_type_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@MutationRate")) {
    dictionary_manager(dictionary_name).ReadString("@MutationRate", mutation_rate_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@CrossoverType")) {
    dictionary_manager(dictionary_name).ReadString("@CrossoverType", crossover_type_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@CrossoverRate")) {
    dictionary_manager(dictionary_name).ReadString("@CrossoverRate", crossover_rate_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@ReplacementType")) {
    dictionary_manager(dictionary_name).ReadString("@ReplacementType", replacement_type_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@Division")) {
    dictionary_manager(dictionary_name).ReadString("@Division", division_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@MaxBoxsizeLimit")) {
    dictionary_manager(dictionary_name).ReadString("@MaxBoxsizeLimit", max_boxsize_limit_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@MinBoxsizeLimit")) {
    dictionary_manager(dictionary_name).ReadString("@MinBoxsizeLimit", min_boxsize_limit_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@DiverseInput")) {
    diverse_input_ = true;
    dictionary_manager(dictionary_name).ReadOption("@DiverseInput", diverse_dakota_input_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@Gradient")) {
    dictionary_manager(dictionary_name).ReadBool("@Gradient", dakota_gradient_);
  }

  if (dictionary_manager(dictionary_name).CheckOption("@EchoDakotaInput")) {
    dictionary_manager(dictionary_name).ReadBool("@EchoDakotaInput", echo_dakota_string_);
  }
}
}  // namespace OptiSMOKE
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
