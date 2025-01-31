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
OptionsOptimizationTargets::OptionsOptimizationTargets() { n_parameters_ = 0; }

OptionsOptimizationTargets::~OptionsOptimizationTargets() {}

void OptionsOptimizationTargets::SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                                                     std::string dictionary_name) {
  dictionary_manager(dictionary_name).SetGrammar(grammar_);

  // ==================================================
  // Direct reactions - specific parameters to be optimized
  if (dictionary_manager(dictionary_name).CheckOption("@TargetA")) {
    dictionary_manager(dictionary_name).ReadOption("@TargetA", list_of_target_lnA_);
  }

  n_parameters_ += list_of_target_lnA_.size();

  if (dictionary_manager(dictionary_name).CheckOption("@TargetBeta")) {
    dictionary_manager(dictionary_name).ReadOption("@TargetBeta", list_of_target_Beta_);
  }

  n_parameters_ += list_of_target_Beta_.size();

  if (dictionary_manager(dictionary_name).CheckOption("@TargetEact")) {
    dictionary_manager(dictionary_name).ReadOption("@TargetEact", list_of_target_E_over_R_);
  }

  n_parameters_ += list_of_target_E_over_R_.size();

  // ==================================================
  // HPL reactions - which
  if (dictionary_manager(dictionary_name).CheckOption("@TargetAinf")) {
    dictionary_manager(dictionary_name).ReadOption("@TargetAinf", list_of_target_lnA_inf_);
  }

  n_parameters_ += list_of_target_lnA_inf_.size();

  if (dictionary_manager(dictionary_name).CheckOption("@TargetBetaInf")) {
    dictionary_manager(dictionary_name).ReadOption("@TargetBetaInf", list_of_target_Beta_inf_);
  }

  n_parameters_ += list_of_target_Beta_inf_.size();

  if (dictionary_manager(dictionary_name).CheckOption("@TargetEactInf")) {
    dictionary_manager(dictionary_name).ReadOption("@TargetEactInf", list_onumberOfParameters_ f_target_E_over_R_inf_);
  }

  n_parameters_ += list_of_target_E_over_R_inf_.size();

  // ==================================================
  // Third body reactions - which
  if (dictionary_manager(dictionary_name).CheckOption("@Target3BodyReactions")) {
    dictionary_manager(dictionary_name).ReadOption("@Target3BodyReactions", list_of_target_thirdbody_reactions_);
  }

  n_parameters_ += list_of_target_thirdbody_reactions_.size();

  // ==================================================
  // Classic PLOG - which
  if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_classic_PLOG_Reactions")) {
    dictionary_manager(dictionary_name)
        .ReadOption("@ListOfTarget_classic_PLOG_Reactions", list_of_target_classic_plog_reactions_);
  }

  n_parameters_ += list_of_target_classic_plog_reactions_.size() * 3;

  if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_classic_PLOG")) {
    dictionary_manager(dictionary_name)
        .ReadOption("@ListOfUncertaintyFactors_classic_PLOG", list_of_uncertainty_factors_classic_plog_);
  }

  // ==================================================
  // Reduced Pressure Based Reactions
  // if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_RPBMR_Reactions"))
  //   dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_RPBMR_Reactions", list_of_target_rpbmr_reactions_);
  //
  // if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_BathGases_RPBMR"))
  //   dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_BathGases_RPBMR", list_of_target_rpbmr_bathgases_);
  //
  // if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_RPBMR"))
  //   dictionary_manager(dictionary_name)
  //       .ReadOption("@ListOfUncertaintyFactors_RPBMR", list_of_uncertainty_factors_rpbmr_);
  //
  // n_parameters_ += list_of_target_rpbmr_bathgases_.size() * 3;

  // ==================================================
  // List of third body species
  if (dictionary_manager(dictionary_name).CheckOption("@Target3BodySpecies")) {
    dictionary_manager(dictionary_name).ReadOption("@Target3BodySpecies", list_of_target_thirdbody_species_);
  }

  // List of target reactions for uncertainty factors
  if (dictionary_manager(dictionary_name).CheckOption("@TargetUncertaintyFactors")) {
    dictionary_manager(dictionary_name).ReadOption("@TargetUncertaintyFactors", list_of_target_uncertainty_factors_);
  }

  // List of uncertainty factors
  if (dictionary_manager(dictionary_name).CheckOption("@UncertaintyFactors")) {
    dictionary_manager(dictionary_name).ReadOption("@UncertaintyFactors", list_of_uncertainty_factors_);
  }

  // List of target inf reactions for uncertainty factors
  if (dictionary_manager(dictionary_name).CheckOption("@UncertaintyFactorsInf")) {
    dictionary_manager(dictionary_name).ReadOption("@UncertaintyFactorsInf", list_of_target_uncertainty_factors_inf_);
  }

  // List of inf uncertainty factors
  if (dictionary_manager(dictionary_name).CheckOption("@UncertaintyFactorsInf")) {
    dictionary_manager(dictionary_name).ReadOption("@UncertaintyFactorsInf", list_of_uncertainty_factors_inf_);
  }

  // FORD
  // if (dictionary_manager(dictionary_name).CheckOption("@ListOf_reactions_FORD"))
  //   dictionary_manager(dictionary_name).ReadOption("@ListOf_reactions_FORD", list_of_ford_);
  // if (dictionary_manager(dictionary_name).CheckOption("@ListOf_species_FORD"))
  //   dictionary_manager(dictionary_name).ReadOption("@ListOf_species_FORD", list_of_species_ford_);
  // if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxAbs_FORD"))
  //   dictionary_manager(dictionary_name).ReadOption("@ListOfMaxAbs_FORD", list_of_max_abs_FORD_);
  // if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinAbs_FORD"))
  //   dictionary_manager(dictionary_name).ReadOption("@ListOfMinAbs_FORD", list_of_min_abs_FORD_);
  // n_parameters_ += list_of_ford_.size();
  // RORD
  // if (dictionary_manager(dictionary_name).CheckOption("@ListOf_reactions_RORD"))
  //   dictionary_manager(dictionary_name).ReadOption("@ListOf_reactions_RORD", list_of_rord_);
  // if (dictionary_manager(dictionary_name).CheckOption("@ListOf_species_RORD"))
  //   dictionary_manager(dictionary_name).ReadOption("@ListOf_species_RORD", list_of_species_rord_);
  // if (dictionary_manager(dictionary_name).CheckOption("@ListOfMaxAbs_RORD"))
  //   dictionary_manager(dictionary_name).ReadOption("@ListOfMaxAbs_RORD", list_of_max_abs_RORD_);
  // if (dictionary_manager(dictionary_name).CheckOption("@ListOfMinAbs_RORD"))
  //   dictionary_manager(dictionary_name).ReadOption("@ListOfMinAbs_RORD", list_of_min_abs_RORD_);
  // n_parameters_ += list_of_rord_.size();

  // if (dictionary_manager(dictionary_name).CheckOption("@ReactionsClassesDefinition")) {
  //   dictionary_manager(dictionary_name).ReadPath("@ReactionsClassesDefinition", reactions_classes_definition_);
  //   if (!fs::exists(reactions_classes_definition_)) {
  //     OptiSMOKE::FatalErrorMessage("The file containing the reactions classes definition does not exists!");
  //   }
  //
  //   ReadReactionClassesDefinition(reactions_classes_definition_);
  // }
}

void OptionsOptimizationTargets::ReadReactionClassesDefinition(fs::path classes_definition) {}

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
