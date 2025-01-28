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

class SimulationsInterface {
 public:
  SimulationsInterface(const OptiSMOKE::InputManager& data);

  ~SimulationsInterface();

  void run();

  void Setup();

  double ComputeObjectiveFunction();

  bool CheckKineticConstasts();

  void SubstituteKineticParameters(const std::vector<double>& c_vars);

  void PrepareASCIIFile(std::ofstream& fOutput,
                        const fs::path output_file_ascii,
                        const std::vector<std::string>& names);

  void PrintASCIIFile(std::ofstream& fOutput, const int eval_nr, const std::vector<double>& b, const double fn_val);

 private:
  const OptiSMOKE::InputManager& data_;

  std::vector<OptiSMOKE::BatchReactor*> batch_reactors;
  std::vector<OptiSMOKE::PlugFlowReactor*> plugflow_reactors;
  std::vector<OptiSMOKE::PerfectlyStirredReactor*> perfectlystirred_reactors;
  std::vector<OptiSMOKE::PremixedLaminarFlame1D*> premixed1D;

  unsigned int n_batch;
  unsigned int n_pfr;
  unsigned int n_psr;
  unsigned int n_premixed;
  unsigned int n_counterflow;

  std::vector<std::vector<double>> k_upper;
  std::vector<std::vector<double>> k_lower;

  std::vector<std::vector<double>> k_upper_inf;
  std::vector<std::vector<double>> k_lower_inf;

  std::vector<std::vector<std::vector<double>>> k_upper_classic_plog;
  std::vector<std::vector<std::vector<double>>> k_lower_classic_plog;

  std::vector<std::vector<std::vector<double>>> k_upper_rpbrm;
  std::vector<std::vector<std::vector<double>>> k_lower_rpbrm;

  std::vector<std::vector<std::vector<double>>> simulations_results_;

  std::vector<double> T_span = {300,  400,  500,  600,  700,  800,  900,  1000, 1100, 1200, 1300, 1400,
                                1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500};

  void ChangeDirectParamaters(std::string type, int index, double parameter);

  void ChangeFallOffParamaters(std::string type, int index, double parameter);

  void ChangeThirdBodyEfficiencies(unsigned int i, std::string name, double parameter);

  void ChangePLOGReactions(std::string type, unsigned int index, double parameter);

  void ChangeRPBRMReactions(std::string type, unsigned int index, double parameter, unsigned int index_coll);

  void ChangeReactionOrder(const std::string& type,
                           const int reaction_index,
                           const unsigned int& species_idx,
                           const double parameter);
};
}  // namespace OptiSMOKE

#include "SimulationsInterface.hpp"
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
