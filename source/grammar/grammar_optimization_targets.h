/*-----------------------------------------------------------------------*\
|     ____            _  ______ __  __  ____  _  ________                 |
|    / __ \       _  (_)/  ___ |  \/  |/ __ \| |/ /  ____|                |
|   | |  | |_ __ | |_ _|  (___ | \  / | |  | | ' /| |__    _     _        |
|   | |  | | '_ \|  _| |\___  \| |\/| | |  | |  < |  __| _| |_ _| |_      |
|   | |__| | |_) | |_| |____)  | |  | | |__| | . \| |___|_   _|_   _|     |
|    \____/| .__/\___|_|______/|_|  |_|\____/|_|\_\______||_|   |_|       |
|          | |                                                            |
|          |_|                                                            |
|                                                                         |
|            Authors: Magnus Fürst <magnus.furst@ulb.ac.be>               |
|                     Andrea Bertolino <andrea.bertolino@ulb.be>          |
|                     Timoteo Dinelli <timoteo.dinelli@polimi.it>         |
|-------------------------------------------------------------------------|
|   License                                                               |
|                                                                         |
|   This file is part of OptiSMOKE.                                       |
|   Copyright (C) 2020 by Magnus Fürst and Andrea Bertolino               |
|                                                                         |
|   OptiSMOKE is free software: you can redistribute it and/or modify     |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OptiSMOKE is distributed in the hope that it will be useful,          |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OptiSMOKE. If not, see <http://www.gnu.org/licenses/>.     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"
#include "dictionary/OpenSMOKE_DictionaryManager.h"

namespace OptiSMOKE
{
class grammar_optimization_targets : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
{
  protected:
    virtual void DefineRules()
    {
        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfBatchReactor", OpenSMOKE::SINGLE_INT,
                                                          "Number of batch reactor datasets", false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfPlugFlowReactor", OpenSMOKE::SINGLE_INT,
                                                          "Number of plug flow reactor datasets", false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfPerfectlyStirredReactor", OpenSMOKE::SINGLE_INT,
                                                          "Number of perfectly stirred reactor datasets", false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfPremixedLaminarFlame", OpenSMOKE::SINGLE_INT,
                                                          "Number of laminar flame datasets", false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfCounterFlowFlame", OpenSMOKE::SINGLE_INT,
                                                          "Number of counterflow flame datasets", false));

        // list of constraints gli posso dare k max e k min con un file di testo
        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfConstraints", OpenSMOKE::SINGLE_STRING,
                                                          "File including the path for list of input files", false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@ListOfTargetUncertaintyFactors", OpenSMOKE::VECTOR_INT,
            "List of reaction indices (starting from 1) for which uncertainty factors are defined", false, "none",
            "@ListOfUncertaintyFactors", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfUncertaintyFactors", OpenSMOKE::VECTOR_DOUBLE,
                                                          "List of uncertainty factors", false, "none",
                                                          "@ListOfTargetUncertaintyFactors", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@ListOfTargetUncertaintyFactors_inf", OpenSMOKE::VECTOR_INT,
            "List of reaction indices (starting from 1) for which uncertainty factors are defined", false, "none",
            "@ListOfUncertaintyFactors_inf", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfUncertaintyFactors_inf", OpenSMOKE::VECTOR_DOUBLE,
                                                          "List of uncertainty factors", false, "none",
                                                          "@ListOfTargetUncertaintyFactors_inf", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@ListOfTarget_lnA", OpenSMOKE::VECTOR_INT,
            "List of target reactions for frequency factors (indices starting from 1)", false, "none", "none", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@ListOfTarget_Beta", OpenSMOKE::VECTOR_INT,
            "List of target reactions for temperature exponents (indices starting from 1)", false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@ListOfTarget_E_over_R", OpenSMOKE::VECTOR_INT,
            "List of target reactions for activation temperatures (indices starting from 1)", false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@ListOfTarget_lnA_inf", OpenSMOKE::VECTOR_INT,
            "List of target reactions for frequency factors (inf) (indices starting from 1)", false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@ListOfTarget_Beta_inf", OpenSMOKE::VECTOR_INT,
            "List of target reactions for temperature exponents (inf) (indices starting from 1)", false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@ListOfTarget_E_over_R_inf", OpenSMOKE::VECTOR_INT,
            "List of target reactions for activation temperatures (inf) (indices starting from 1)", false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_ThirdBody_Reactions", OpenSMOKE::VECTOR_INT,
                                                          "List of target third body reactions", false, "none",
                                                          "@ListOfTarget_ThirdBody_Species", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_ThirdBody_Species", OpenSMOKE::VECTOR_STRING,
                                                          "List of target third body species", false, "none",
                                                          "@ListOfTarget_ThirdBody_Reactions", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_ThirdBody_Eff", OpenSMOKE::VECTOR_STRING,
                                                          "List of maximum values for the third body efficiencies",
                                                          false, "none", "@ListOfTarget_ThirdBody_Reactions", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_ThirdBody_Eff", OpenSMOKE::VECTOR_STRING,
                                                          "List of minimum values for the third body efficiencies",
                                                          false, "none", "@ListOfTarget_ThirdBody_Reactions", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_classic_PLOG_Reactions", OpenSMOKE::VECTOR_INT,
                                                          "List of target classic pressure logarithmic reactions",
                                                          false, "none", "@ListOfUncertaintyFactors_classic_PLOG",
                                                          "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfUncertaintyFactors_classic_PLOG",
                                                          OpenSMOKE::VECTOR_DOUBLE,
                                                          "List of uncertainty factors for classic PLOG", false, "none",
                                                          "@ListOfTarget_classic_PLOG_Reactions", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_RPBMR_Reactions", OpenSMOKE::VECTOR_INT,
                                                          "TODO", false, "none", "@ListOfUncertaintyFactors_RPBMR",
                                                          "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfUncertaintyFactors_RPBMR", OpenSMOKE::VECTOR_DOUBLE,
                                                          "TODO", false, "none", "@ListOfTarget_RPBMR_Reactions",
                                                          "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_BathGases_RPBMR", OpenSMOKE::VECTOR_STRING,
                                                          "TODO", false, "none", "@ListOfTarget_RPBMR_Reactions",
                                                          "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@ReactionsClassesDefinitions", OpenSMOKE::SINGLE_PATH,
            "Path to the file containing the definitions of the reaction/s classes.", false));
    }
};
} // namespace OptiSMOKE
