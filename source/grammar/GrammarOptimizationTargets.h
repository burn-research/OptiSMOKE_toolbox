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
class GrammarOptimizationTargets : public OpenSMOKE::OpenSMOKE_DictionaryGrammar {
 protected:
  virtual void DefineRules() {
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfBatchReactor",
                                                      OpenSMOKE::SINGLE_INT,
                                                      "Number of batch reactor datasets",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfPlugFlowReactor",
                                                      OpenSMOKE::SINGLE_INT,
                                                      "Number of plug flow reactor datasets",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfPerfectlyStirredReactor",
                                                      OpenSMOKE::SINGLE_INT,
                                                      "Number of perfectly stirred reactor datasets",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfPremixedLaminarFlame",
                                                      OpenSMOKE::SINGLE_INT,
                                                      "Number of laminar flame datasets",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfCounterFlowFlame",
                                                      OpenSMOKE::SINGLE_INT,
                                                      "Number of counterflow flame datasets",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfConstraints",
                                                      OpenSMOKE::SINGLE_STRING,
                                                      "File including the path for list of input files",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@ListOfTargetUncertaintyFactors",
        OpenSMOKE::VECTOR_INT,
        "List of reaction indices (starting from 1) for which uncertainty factors are defined",
        false,
        "none",
        "@ListOfUncertaintyFactors",
        "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfUncertaintyFactors",
                                                      OpenSMOKE::VECTOR_DOUBLE,
                                                      "List of uncertainty factors",
                                                      false,
                                                      "none",
                                                      "@ListOfTargetUncertaintyFactors",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@ListOfTargetUncertaintyFactors_inf",
        OpenSMOKE::VECTOR_INT,
        "List of reaction indices (starting from 1) for which uncertainty factors are defined",
        false,
        "none",
        "@ListOfUncertaintyFactors_inf",
        "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfUncertaintyFactors_inf",
                                                      OpenSMOKE::VECTOR_DOUBLE,
                                                      "List of uncertainty factors",
                                                      false,
                                                      "none",
                                                      "@ListOfTargetUncertaintyFactors_inf",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@ListOfTarget_lnA",
        OpenSMOKE::VECTOR_INT,
        "List of target reactions for frequency factors (indices starting from 1)",
        false,
        "none",
        "none",
        "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@ListOfTarget_Beta",
        OpenSMOKE::VECTOR_INT,
        "List of target reactions for temperature exponents (indices starting from 1)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@ListOfTarget_E_over_R",
        OpenSMOKE::VECTOR_INT,
        "List of target reactions for activation temperatures (indices starting from 1)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@ListOfTarget_lnA_inf",
        OpenSMOKE::VECTOR_INT,
        "List of target reactions for frequency factors (inf) (indices starting from 1)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@ListOfTarget_Beta_inf",
        OpenSMOKE::VECTOR_INT,
        "List of target reactions for temperature exponents (inf) (indices starting from 1)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@ListOfTarget_E_over_R_inf",
        OpenSMOKE::VECTOR_INT,
        "List of target reactions for activation temperatures (inf) (indices starting from 1)",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_ThirdBody_Reactions",
                                                      OpenSMOKE::VECTOR_INT,
                                                      "List of target third body reactions",
                                                      false,
                                                      "none",
                                                      "@ListOfTarget_ThirdBody_Species",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_ThirdBody_Species",
                                                      OpenSMOKE::VECTOR_STRING,
                                                      "List of target third body species",
                                                      false,
                                                      "none",
                                                      "@ListOfTarget_ThirdBody_Reactions",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_ThirdBody_Eff",
                                                      OpenSMOKE::VECTOR_STRING,
                                                      "List of maximum values for the third body efficiencies",
                                                      false,
                                                      "none",
                                                      "@ListOfTarget_ThirdBody_Reactions",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_ThirdBody_Eff",
                                                      OpenSMOKE::VECTOR_STRING,
                                                      "List of minimum values for the third body efficiencies",
                                                      false,
                                                      "none",
                                                      "@ListOfTarget_ThirdBody_Reactions",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_classic_PLOG_Reactions",
                                                      OpenSMOKE::VECTOR_INT,
                                                      "List of target classic pressure logarithmic reactions",
                                                      false,
                                                      "none",
                                                      "@ListOfUncertaintyFactors_classic_PLOG",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfUncertaintyFactors_classic_PLOG",
                                                      OpenSMOKE::VECTOR_DOUBLE,
                                                      "List of uncertainty factors for classic PLOG",
                                                      false,
                                                      "none",
                                                      "@ListOfTarget_classic_PLOG_Reactions",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_RPBMR_Reactions",
                                                      OpenSMOKE::VECTOR_INT,
                                                      "TODO",
                                                      false,
                                                      "none",
                                                      "@ListOfUncertaintyFactors_RPBMR",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfUncertaintyFactors_RPBMR",
                                                      OpenSMOKE::VECTOR_DOUBLE,
                                                      "TODO",
                                                      false,
                                                      "none",
                                                      "@ListOfTarget_RPBMR_Reactions",
                                                      "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_BathGases_RPBMR",
                                                      OpenSMOKE::VECTOR_STRING,
                                                      "TODO",
                                                      false,
                                                      "none",
                                                      "@ListOfTarget_RPBMR_Reactions",
                                                      "none"));
    // FORD
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOf_reactions_FORD",
                                                      OpenSMOKE::VECTOR_INT,
                                                      "TODO",
                                                      false,
                                                      "none",
                                                      "none",
                                                      "none"));
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOf_species_FORD",
                                                      OpenSMOKE::VECTOR_STRING,
                                                      "TODO",
                                                      false,
                                                      "none",
                                                      "none",
                                                      "none"));
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_FORD",
                                                      OpenSMOKE::VECTOR_DOUBLE,
                                                      "TODO",
                                                      false,
                                                      "none",
                                                      "none",
                                                      "none"));
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_FORD",
                                                      OpenSMOKE::VECTOR_DOUBLE,
                                                      "TODO",
                                                      false,
                                                      "none",
                                                      "none",
                                                      "none"));
    // RORD
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOf_reactions_RORD",
                                                      OpenSMOKE::VECTOR_INT,
                                                      "TODO",
                                                      false,
                                                      "none",
                                                      "none",
                                                      "none"));
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOf_species_RORD",
                                                      OpenSMOKE::VECTOR_STRING,
                                                      "TODO",
                                                      false,
                                                      "none",
                                                      "none",
                                                      "none"));
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_RORD",
                                                      OpenSMOKE::VECTOR_DOUBLE,
                                                      "TODO",
                                                      false,
                                                      "none",
                                                      "none",
                                                      "none"));
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_RORD",
                                                      OpenSMOKE::VECTOR_DOUBLE,
                                                      "TODO",
                                                      false,
                                                      "none",
                                                      "none",
                                                      "none"));

    AddKeyWord(
        OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ReactionsClassesDefinitions",
                                               OpenSMOKE::SINGLE_PATH,
                                               "Path to the file containing the definitions of the reaction/s classes.",
                                               false));
  }
};
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
