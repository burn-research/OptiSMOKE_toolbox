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

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OptiSMOKE
{
	class grammar_optimization_targets : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
            // To be removed
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TypeOfReactor",
                OpenSMOKE::SINGLE_STRING,
                "Type of ideal reactor",
                false,
                "none",
                "none",
                "@NumberOfBatchReactor @NumberOfPlugFlowReactor @NumberOfPerfectlyStirredReactor @NumberOfPremixedLaminarFlame @NumberOfCounterFlowFlame"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfBatchReactor",
                OpenSMOKE::SINGLE_INT,
                "Number of batch reactor datasets",
                false,
                "none",
                "none",
                "@TypeOfReactor"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfPlugFlowReactor",
                OpenSMOKE::SINGLE_INT,
                "Number of plug flow reactor datasets",
                false,
                "none",
                "none",
                "@TypeOfReactor"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfPerfectlyStirredReactor",
                OpenSMOKE::SINGLE_INT,
                "Number of perfectly stirred reactor datasets",
                false,
                "none",
                "none",
                "@TypeOfReactor"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfPremixedLaminarFlame",
                OpenSMOKE::SINGLE_INT,
                "Number of laminar flame datasets",
                false,
                "none",
                "none",
                "@TypeOfReactor"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfCounterFlowFlame",
                OpenSMOKE::SINGLE_INT,
                "Number of counterflow flame datasets",
                false,
                "none",
                "none",
                "@TypeOfReactor"));

            // list of constraints gli posso dare k max e k min con un file di testo  
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfConstraints",
                OpenSMOKE::SINGLE_STRING,
                "File including the path for list of input files",
                false,
                "none",
                "none",
                "none"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTargetUncertaintyFactors",
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

            // limite a pressione infinito TROE se metto inf prende il limite di pressione infinita 
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTargetUncertaintyFactors_inf",
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

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_A",
                OpenSMOKE::VECTOR_INT,
                "List of target reactions for frequency factors (indices starting from 1)",
                false,
                "none",
                "none",
                "@ListOfTarget_lnA"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_lnA",
                OpenSMOKE::VECTOR_INT,
		        "List of target reactions for frequency factors (indices starting from 1)",
		        false,
                "none",
                "none",
                "@ListOfTarget_A"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_Beta",
		        OpenSMOKE::VECTOR_INT,
		        "List of target reactions for temperature exponents (indices starting from 1)",
		        false,
                "none",
                "none",
                "none"));
            
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_E",
                OpenSMOKE::VECTOR_INT,
                "List of target reactions for activation temperatures (indices starting from 1)",
                false,
                "none",
                "none",
                "@ListOfTarget_E_over_R"));

	        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_E_over_R",
		        OpenSMOKE::VECTOR_INT,
		        "List of target reactions for activation temperatures (indices starting from 1)",
		        false,
                "none",
                "none",
                "@ListOfTarget_E"));
            
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_A_inf",
                OpenSMOKE::VECTOR_INT,
                "List of target reactions for frequency factors (inf) (indices starting from 1)",
                false,
                "none",
                "none",
                "@ListOfTarget_lnA_inf"));
            
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_lnA_inf",
                OpenSMOKE::VECTOR_INT,
                "List of target reactions for frequency factors (inf) (indices starting from 1)",
                false,
                "none",
                "none",
                "@ListOfTarget_A_inf"));
            
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_Beta_inf",
                OpenSMOKE::VECTOR_INT,
                "List of target reactions for temperature exponents (inf) (indices starting from 1)",
                false,
                "none",
                "none",
                "none"));
            
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_E_inf",
                OpenSMOKE::VECTOR_INT,
                "List of target reactions for activation temperatures (inf) (indices starting from 1)",
                false,
                "none",
                "none",
                "@ListOfTarget_E_over_R_inf"));
            
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_E_over_R_inf",
                OpenSMOKE::VECTOR_INT,
                "List of target reactions for activation temperatures (inf) (indices starting from 1)",
                false,
                "none",
                "none",
                "@ListOfTarget_E_inf"));

		    // H+O2 TROE limite di alta e bassa pressione con thirdbody (questa è per ottimizzare le efficienze di terzo corpo)
	        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_ThirdBody_Reactions",
                OpenSMOKE::VECTOR_INT,
		        "List of target third body reactions",
		        false,
                "none",
                "@ListOfTarget_ThirdBody_Species",
                "none"));

		    // qua nelle thirdbody gli dici quale specie voglio ottimizzare
	        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_ThirdBody_Species",
		        OpenSMOKE::VECTOR_STRING,
		        "List of target third body species",
		        false,
                "none",
                "@ListOfTarget_ThirdBody_Reactions",
                "none"));

            // Extended PLOG reactions for THIRD BODY ESTIMATION
	        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_ExtPLOG_Reactions_TB",
		        OpenSMOKE::VECTOR_INT,
		        "Contains a list of indices of the extended PLOG reactions for the optimization of Third Bodies",
		        false,
                "none",
                "none",
                "none"));

	        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_ExtPLOG_Species",
		        OpenSMOKE::VECTOR_STRING,
		        "List of target species in extendend pressure logarithmic reactions",
		        false,
                "none",
                "@ListOfTarget_ExtPLOG_Reactions_TB",
                "none"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMin_TBeff_ExtPLOG",
                OpenSMOKE::VECTOR_DOUBLE,
                "List of minimum absolute values for the third body efficiencies of colliders",
                false,
                "none",
                "@ListOfTarget_ExtPLOG_Reactions_TB",
                "none"));
        
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMax_TBeff_ExtPLOG",
                OpenSMOKE::VECTOR_DOUBLE,
                "List of maximum absolute values for the third body efficiencies of colliders",
                false,
                "none",
                "@ListOfTarget_ExtPLOG_Reactions_TB",
                "none"));

	        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_ExtPLOG_Reactions",
                OpenSMOKE::VECTOR_INT,
		        "List of target extended pressure logarithmic reactions",
		        false,
                "none",
                "@ListOfUncertaintyFactors_ExtPLOG",
                "none"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfUncertaintyFactors_ExtPLOG",
                OpenSMOKE::VECTOR_DOUBLE,
                "List of uncertainty factors for extended pressure logarithmic reactions",
                false,
                "none",
                "@ListOfTarget_ExtPLOG_Reactions",
                "none"));
            // Extended PLOG reactions for THIRD BODY ESTIMATION

            // Extended PLOG reactions
	        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_EPLR",
		        OpenSMOKE::VECTOR_INT,
		        "Contains a list of indices of the extended PLOG reaction for their optimization",
		        false,
                "none",
                "none",
                "none"));

	        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfTarget_BathGases_EPLR",
		        OpenSMOKE::VECTOR_STRING,
		        "Contains a list of species of interest for the optimization of their corresponding rate",
		        false,
                "none",
                "@ListOfTarget_EPLR",
                "none"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfUncertaintyFactors_EPLR",
                OpenSMOKE::VECTOR_DOUBLE,
                "List of uncertainty factors for EPLR",
                false,
                "none",
                "@ListOfTarget_EPLR",
                "none"));
            // Extended PLOG reactions

            // PLOG
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
            // PLOG

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ReactionsClassesDefinitions",
                OpenSMOKE::SINGLE_PATH,
                "Path to the file containing the definitions of the reaction/s classes.",
                false,
                "none",
                "none",
                "none"));

        }
    };
} // namespace OptiSMOKE