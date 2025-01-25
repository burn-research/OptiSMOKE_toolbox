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

namespace OptiSMOKE {
OptionsOptimizationTargets::OptionsOptimizationTargets() {
  numberOfBatchReactor_ = 0;
  numberOfPlugFlowReactor_ = 0;
  numberOfPerfectlyStirredReactor_ = 0;
  numberOfPremixedLaminarFlame_ = 0;
  numberOfCounterFlowFlame_ = 0;
  numberOfKTExperiments_ = 0;
  numberOfKTPExperiments_ = 0;

  numberOfParameters_ = 0;
}

OptionsOptimizationTargets::~OptionsOptimizationTargets() {}

void OptionsOptimizationTargets::SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                                                     std::string dictionary_name) {
  dictionary_manager(dictionary_name).SetGrammar(grammar_);

  if (dictionary_manager(dictionary_name).CheckOption("@NumberOfBatchReactor"))
    dictionary_manager(dictionary_name).ReadInt("@NumberOfBatchReactor", numberOfBatchReactor_);

  if (dictionary_manager(dictionary_name).CheckOption("@NumberOfPlugFlowReactor"))
    dictionary_manager(dictionary_name).ReadInt("@NumberOfPlugFlowReactor", numberOfPlugFlowReactor_);

  if (dictionary_manager(dictionary_name).CheckOption("@NumberOfPerfectlyStirredReactor"))
    dictionary_manager(dictionary_name).ReadInt("@NumberOfPerfectlyStirredReactor", numberOfPerfectlyStirredReactor_);

  if (dictionary_manager(dictionary_name).CheckOption("@NumberOfPremixedLaminarFlame"))
    dictionary_manager(dictionary_name).ReadInt("@NumberOfPremixedLaminarFlame", numberOfPremixedLaminarFlame_);

  if (dictionary_manager(dictionary_name).CheckOption("@NumberOfCounterFlowFlame"))
    dictionary_manager(dictionary_name).ReadInt("@NumberOfCounterFlowFlame", numberOfCounterFlowFlame_);

  if (dictionary_manager(dictionary_name).CheckOption("@NumberOfKTExperiments"))
    dictionary_manager(dictionary_name).ReadInt("@NumberOfKTExperiments", numberOfKTExperiments_);

  if (dictionary_manager(dictionary_name).CheckOption("@NumberOfKTPExperiments"))
    dictionary_manager(dictionary_name).ReadInt("@NumberOfKTPExperiments", numberOfKTPExperiments_);

  // Direct reactions - specific parameters to be optimized
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_lnA"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_lnA", list_of_target_lnA_);

  numberOfParameters_ += list_of_target_lnA_.size();

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_Beta"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_Beta", list_of_target_Beta_);

  numberOfParameters_ += list_of_target_Beta_.size();

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_E_over_R"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_E_over_R", list_of_target_E_over_R_);

  numberOfParameters_ += list_of_target_E_over_R_.size();

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_lnA_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_lnA_inf", list_of_target_lnA_inf_);

  numberOfParameters_ += list_of_target_lnA_inf_.size();

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_Beta_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_Beta_inf", list_of_target_Beta_inf_);

  numberOfParameters_ += list_of_target_Beta_inf_.size();

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_E_over_R_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_E_over_R_inf", list_of_target_E_over_R_inf_);

  numberOfParameters_ += list_of_target_E_over_R_inf_.size();

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ThirdBody_Reactions"))
    dictionary_manager(dictionary_name)
        .ReadOption("@ListOfTarget_ThirdBody_Reactions", list_of_target_thirdbody_reactions_);

  numberOfParameters_ += list_of_target_thirdbody_reactions_.size();

  // Classic PLOG - which
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_classic_PLOG_Reactions"))
    dictionary_manager(dictionary_name)
        .ReadOption("@ListOfTarget_classic_PLOG_Reactions", list_of_target_classic_plog_reactions_);

  numberOfParameters_ += list_of_target_classic_plog_reactions_.size() * 3;

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_classic_PLOG"))
    dictionary_manager(dictionary_name)
        .ReadOption("@ListOfUncertaintyFactors_classic_PLOG", list_of_uncertainty_factors_classic_plog_);

  // Reduced Pressure Based Reactions
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_RPBMR_Reactions"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_RPBMR_Reactions", list_of_target_rpbmr_reactions_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_BathGases_RPBMR"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_BathGases_RPBMR", list_of_target_rpbmr_bathgases_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_RPBMR"))
    dictionary_manager(dictionary_name)
        .ReadOption("@ListOfUncertaintyFactors_RPBMR", list_of_uncertainty_factors_rpbmr_);

  numberOfParameters_ += list_of_target_rpbmr_bathgases_.size() * 3;

  // List of third body species
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ThirdBody_Species"))
    dictionary_manager(dictionary_name)
        .ReadOption("@ListOfTarget_ThirdBody_Species", list_of_target_thirdbody_species_);

  // List of target reactions for uncertainty factors
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTargetUncertaintyFactors"))
    dictionary_manager(dictionary_name)
        .ReadOption("@ListOfTargetUncertaintyFactors", list_of_target_uncertainty_factors_);

  // List of uncertainty factors
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors", list_of_uncertainty_factors_);

  // List of target inf reactions for uncertainty factors
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTargetUncertaintyFactors_inf"))
    dictionary_manager(dictionary_name)
        .ReadOption("@ListOfTargetUncertaintyFactors_inf", list_of_target_uncertainty_factors_inf_);

  // List of inf uncertainty factors
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors_inf", list_of_uncertainty_factors_inf_);

  // List of relative maximum parameters
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxRel_lnA"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxRel_lnA", list_of_max_rel_lnA_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxRel_Beta"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxRel_Beta", list_of_max_rel_Beta_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxRel_E_over_R"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxRel_E_over_R", list_of_max_rel_E_over_R_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxRel_lnA_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxRel_lnA_inf", list_of_max_rel_lnA_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxRel_Beta_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxRel_Beta_inf", list_of_max_rel_Beta_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxRel_E_over_R_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxRel_E_over_R_inf", list_of_max_rel_E_over_R_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxRel_ThirdBody_Eff"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxRel_ThirdBody_Eff", list_of_max_rel_thirdbody_eff_);

  // List of relative minimum parameters
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinRel_lnA"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinRel_lnA", list_of_min_rel_lnA_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinRel_Beta"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinRel_Beta", list_of_min_rel_Beta_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinRel_E_over_R"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinRel_E_over_R", list_of_min_rel_E_over_R_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinRel_lnA_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinRel_lnA_inf", list_of_min_rel_lnA_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinRel_Beta_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinRel_Beta_inf", list_of_min_rel_Beta_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinRel_E_over_R_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinRel_E_over_R_inf", list_of_min_rel_E_over_R_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinRel_ThirdBody_Eff"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinRel_ThirdBody_Eff", list_of_min_rel_thirdbody_eff_);

  // List of absolute maximum parameters
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxAbs_lnA"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxAbs_lnA", list_of_max_abs_lnA_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxAbs_Beta"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxAbs_Beta", list_of_max_abs_Beta_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxAbs_E_over_R"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxAbs_E_over_R", list_of_max_abs_E_over_R_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxAbs_lnA_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxAbs_lnA_inf", list_of_max_abs_lnA_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxAbs_Beta_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxAbs_Beta_inf", list_of_max_abs_Beta_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxAbs_E_over_R_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxAbs_E_over_R_inf", list_of_max_abs_E_over_R_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxAbs_ThirdBody_Eff"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxAbs_ThirdBody_Eff", list_of_max_abs_thirdbody_eff_);

  // List of absolute minimum parameters
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinAbs_Beta"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinAbs_Beta", list_of_min_abs_Beta_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinAbs_E_over_R"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinAbs_E_over_R", list_of_min_abs_E_over_R_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinAbs_lnA_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinAbs_lnA_inf", list_of_min_abs_lnA_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinAbs_Beta_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinAbs_Beta_inf", list_of_min_abs_Beta_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinAbs_E_over_R_inf"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinAbs_E_over_R_inf", list_of_min_abs_E_over_R_inf_);

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinAbs_ThirdBody_Eff"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinAbs_ThirdBody_Eff", list_of_min_abs_thirdbody_eff_);

  // FORD
  if (dictionary_manager(dictionary_name).CheckOption("@ListOf_reactions_FORD"))
    dictionary_manager(dictionary_name).ReadOption("@ListOf_reactions_FORD", list_of_ford_);
  if (dictionary_manager(dictionary_name).CheckOption("@ListOf_species_FORD"))
    dictionary_manager(dictionary_name).ReadOption("@ListOf_species_FORD", list_of_species_ford_);
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxAbs_FORD"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxAbs_FORD", list_of_max_abs_FORD_);
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinAbs_FORD"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinAbs_FORD", list_of_min_abs_FORD_);

  numberOfParameters_ += list_of_ford_.size();

  // RORD
  if (dictionary_manager(dictionary_name).CheckOption("@ListOf_reactions_RORD"))
    dictionary_manager(dictionary_name).ReadOption("@ListOf_reactions_RORD", list_of_rord_);
  if (dictionary_manager(dictionary_name).CheckOption("@ListOf_species_RORD"))
    dictionary_manager(dictionary_name).ReadOption("@ListOf_species_RORD", list_of_species_rord_);
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxAbs_RORD"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMaxAbs_RORD", list_of_max_abs_RORD_);
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinAbs_RORD"))
    dictionary_manager(dictionary_name).ReadOption("@ListOfMinAbs_RORD", list_of_min_abs_RORD_);

  numberOfParameters_ += list_of_rord_.size();

  if (dictionary_manager(dictionary_name).CheckOption("@ReactionsClassesDefinitions")) {
    dictionary_manager(dictionary_name).ReadPath("@ReactionsClassesDefinitions", reactions_classes_definition_);
    if (!fs::exists(reactions_classes_definition_))
      OptiSMOKE::FatalErrorMessage("The @ReactionsClassesDefinitions path does not exists!");

    ReadReactionClassesDefinition(reactions_classes_definition_);
  }
}

void OptionsOptimizationTargets::ReadReactionClassesDefinition(fs::path classes_definition) {}

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
