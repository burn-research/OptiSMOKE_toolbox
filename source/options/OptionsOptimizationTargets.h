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
|           Authors: Timoteo Dinelli  <timoteo.dinelli@polimi.it>                   |
|                    Andrea Bertolino <andrea.bertolino@ulb.be>                     |
|                    Magnus Fürst     <magnus.furst@ulb.ac.be>                      |
|                                                                                   |
|             [1] CRECK Modeling Lab <https://www.creckmodeling.polimi.it>          |
|                 Department of Chemistry, Materials and Chemical Engineering       |
|                 Politecnico di Milano, P.zza Leonardo da Vinci 32, 20133 Milano   |
|                                                                                   |
|             [2] BRITE Research Group <https://brite-research.be>                  |
|                 Brussels Institute for Thermal-fluid systems and clean Energy     |
|                 Avenue F.D. Rooseveltlaan 50, Bruxelles 1050 Brussel              |
|                                                                                   |
\* ------------------------------------------------------------------------------- */
#pragma once

namespace OptiSMOKE {
class OptionsOptimizationTargets {
 public:
  OptionsOptimizationTargets();

  ~OptionsOptimizationTargets();

  void SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager, std::string dictionary_name);

  const size_t& NumberOfBatchReactors() const { return numberOfBatchReactor_; };
  const size_t& NumberOfPlugFlowReactors() const { return numberOfPlugFlowReactor_; };
  const size_t& NumberOfPerfectlyStirredReactors() const { return numberOfPerfectlyStirredReactor_; };
  const size_t& NumberOfPremixedFlames() const { return numberOfPremixedLaminarFlame_; };
  const size_t& NumberOfCounterFlowFlames() const { return numberOfCounterFlowFlame_; };
  const size_t& NumberOfKTExperiments() const { return number_of_KT_experiments_; };
  const size_t& NumberOfKTPExperiments() const { return number_of_KTP_experiments_; };

  const size_t& NumberOfParameters() const { return numberOfParameters_; };

  const std::vector<int>& list_of_target_lnA() const { return list_of_target_lnA_; };
  const std::vector<int>& list_of_target_Beta() const { return list_of_target_Beta_; };
  const std::vector<int>& list_of_target_E_over_R() const { return list_of_target_E_over_R_; };

  const std::vector<int>& list_of_target_lnA_inf() const { return list_of_target_lnA_inf_; };
  const std::vector<int>& list_of_target_Beta_inf() const { return list_of_target_Beta_inf_; };
  const std::vector<int>& list_of_target_E_over_R_inf() const { return list_of_target_E_over_R_inf_; };

  const std::vector<int>& list_of_target_thirdbody_reactions() const { return list_of_target_thirdbody_reactions_; };
  const std::vector<std::string>& list_of_target_thirdbody_species() const {
    return list_of_target_thirdbody_species_;
  };

  const std::vector<int>& list_of_target_classic_plog_reactions() const {
    return list_of_target_classic_plog_reactions_;
  };
  const std::vector<double>& list_of_uncertainty_factors_classic_plog() const {
    return list_of_uncertainty_factors_classic_plog_;
  };

  const std::vector<int>& list_of_target_uncertainty_factors() const { return list_of_target_uncertainty_factors_; };
  const std::vector<double>& list_of_uncertainty_factors() const { return list_of_uncertainty_factors_; };

  const std::vector<int>& list_of_target_uncertainty_factors_inf() const {
    return list_of_target_uncertainty_factors_inf_;
  };
  const std::vector<double>& list_of_uncertainty_factors_inf() const { return list_of_uncertainty_factors_inf_; };

  const std::vector<double>& list_of_min_rel_lnA() const { return list_of_min_rel_lnA_; };
  const std::vector<double>& list_of_max_rel_lnA() const { return list_of_max_rel_lnA_; };

  const std::vector<double>& list_of_min_rel_Beta() const { return list_of_min_rel_Beta_; };
  const std::vector<double>& list_of_max_rel_Beta() const { return list_of_max_rel_Beta_; };

  const std::vector<double>& list_of_min_rel_E_over_R() const { return list_of_min_rel_E_over_R_; };
  const std::vector<double>& list_of_max_rel_E_over_R() const { return list_of_max_rel_E_over_R_; };

  const std::vector<double>& list_of_min_rel_lnA_inf() const { return list_of_min_rel_lnA_inf_; };
  const std::vector<double>& list_of_max_rel_lnA_inf() const { return list_of_max_rel_lnA_inf_; };

  const std::vector<double>& list_of_min_rel_Beta_inf() const { return list_of_min_rel_Beta_inf_; };
  const std::vector<double>& list_of_max_rel_Beta_inf() const { return list_of_max_rel_Beta_inf_; };

  const std::vector<double>& list_of_min_rel_E_over_R_inf() const { return list_of_min_rel_E_over_R_inf_; };
  const std::vector<double>& list_of_max_rel_E_over_R_inf() const { return list_of_max_rel_E_over_R_inf_; };

  const std::vector<double>& list_of_min_rel_thirdbody_eff() const { return list_of_min_rel_thirdbody_eff_; };
  const std::vector<double>& list_of_max_rel_thirdbody_eff() const { return list_of_max_rel_thirdbody_eff_; };

  const std::vector<std::string>& list_of_min_abs_thirdbody_eff() const { return list_of_min_abs_thirdbody_eff_; };
  const std::vector<std::string>& list_of_max_abs_thirdbody_eff() const { return list_of_max_abs_thirdbody_eff_; };

  const std::vector<int>& list_of_target_rpbmr_reactions() const { return list_of_target_rpbmr_reactions_; };

  const std::vector<double>& list_of_uncertainty_factors_rpbmr() const { return list_of_uncertainty_factors_rpbmr_; };

  const std::vector<std::string>& list_of_target_rpbmr_bathgases() const { return list_of_target_rpbmr_bathgases_; }

  const std::vector<int>& list_of_ford() const { return list_of_ford_; };
  const std::vector<std::string>& list_of_species_ford() const { return list_of_species_ford_; };
  const std::vector<double>& list_of_max_abs_FORD() const { return list_of_max_abs_FORD_; };
  const std::vector<double>& list_of_min_abs_FORD() const { return list_of_min_abs_FORD_; };

  const std::vector<int>& list_of_rord() const { return list_of_rord_; };
  const std::vector<std::string>& list_of_species_rord() const { return list_of_species_rord_; };
  const std::vector<double>& list_of_max_abs_RORD() const { return list_of_max_abs_RORD_; };
  const std::vector<double>& list_of_min_abs_RORD() const { return list_of_min_abs_RORD_; };

 private:
  GrammarOptimizationTargets grammar_;

  void ReadReactionClassesDefinition(fs::path reaction_classes);

  int numberOfBatchReactor_;
  int numberOfPlugFlowReactor_;
  int numberOfPerfectlyStirredReactor_;
  int numberOfPremixedLaminarFlame_;
  int numberOfCounterFlowFlame_;
  int numberOfKTExperiments_;
  int numberOfKTPExperiments_;
  int numberOfParameters_;

  fs::path reactions_classes_definition_;

  std::vector<int> list_of_target_lnA_;
  std::vector<int> list_of_target_Beta_;
  std::vector<int> list_of_target_E_over_R_;

  std::vector<int> list_of_target_lnA_inf_;
  std::vector<int> list_of_target_Beta_inf_;
  std::vector<int> list_of_target_E_over_R_inf_;

  std::vector<int> list_of_target_thirdbody_reactions_;
  std::vector<std::string> list_of_target_thirdbody_species_;

  std::vector<int> list_of_target_classic_plog_reactions_;
  std::vector<double> list_of_uncertainty_factors_classic_plog_;

  std::vector<int> list_of_target_uncertainty_factors_;
  std::vector<double> list_of_uncertainty_factors_;

  std::vector<int> list_of_target_uncertainty_factors_inf_;
  std::vector<double> list_of_uncertainty_factors_inf_;

  std::vector<int> list_of_target_rpbmr_reactions_;
  std::vector<double> list_of_uncertainty_factors_rpbmr_;
  std::vector<std::string> list_of_target_rpbmr_bathgases_;

  std::vector<double> list_of_min_rel_lnA_;
  std::vector<double> list_of_max_rel_lnA_;

  std::vector<double> list_of_min_rel_Beta_;
  std::vector<double> list_of_max_rel_Beta_;

  std::vector<double> list_of_min_rel_E_over_R_;
  std::vector<double> list_of_max_rel_E_over_R_;

  std::vector<double> list_of_min_rel_lnA_inf_;
  std::vector<double> list_of_max_rel_lnA_inf_;

  std::vector<double> list_of_min_rel_Beta_inf_;
  std::vector<double> list_of_max_rel_Beta_inf_;

  std::vector<double> list_of_min_rel_E_over_R_inf_;
  std::vector<double> list_of_max_rel_E_over_R_inf_;

  std::vector<double> list_of_min_rel_thirdbody_eff_;
  std::vector<double> list_of_max_rel_thirdbody_eff_;

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

  std::vector<std::string> list_of_min_abs_thirdbody_eff_;
  std::vector<std::string> list_of_max_abs_thirdbody_eff_;

  // FORD
  std::vector<int> list_of_ford_;
  std::vector<std::string> list_of_species_ford_;
  std::vector<double> list_of_max_abs_FORD_;
  std::vector<double> list_of_min_abs_FORD_;

  // RORD
  std::vector<int> list_of_rord_;
  std::vector<std::string> list_of_species_rord_;
  std::vector<double> list_of_max_abs_RORD_;
  std::vector<double> list_of_min_abs_RORD_;
};
}  // namespace OptiSMOKE
#include "OptionsOptimizationTargets.hpp"
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
