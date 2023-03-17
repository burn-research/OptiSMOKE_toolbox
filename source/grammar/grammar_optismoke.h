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

#ifndef GRAMMAR_OPTISMOKEPP_H
#define GRAMMAR_OPTISMOKEPP_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OptiSMOKE
{
	class grammar_optismoke : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

	virtual void DefineRules()
	{

                AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Debug",
                        OpenSMOKE::SINGLE_BOOL,
                        "Specify if debugging print outs should be used (default: false)",
                        false) );

                AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DebugSim",
                        OpenSMOKE::SINGLE_BOOL,
                        "Specify if debugging print outs should be used (default: false)",
                        false) );

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DakotaOptions", 
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary with options for Dakota", 
                        false,
                        "none",
                        "@OptimizationLibrary",
                        "@NLOPTOptions"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NLOPTOptions", 
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary with options for NLOPT", 
                        false,
                        "none",
                        "@OptimizationLibrary",
                        "@DakotaOptions"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CurveMatchingOptions", 
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary with options for setting Curve Matching as the objective function", 
                        false,
                        "none",
                        "none",
                        "none"));

                // TODO revision based on dakota
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OptimizationSetup",
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary to set the options for the optimization",
                        true,
                        "none",
                        "none",
                        "none"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OptimizationTarget",
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary to set the targets for the optimization",
                        true,
                        "none",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsFolder",
                        OpenSMOKE::SINGLE_PATH,
                        "Name of the folder containing the kinetic scheme (XML Version)",
                        true,
                        "@KineticsPreProcessor",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NominalKineticsFolder",
                        OpenSMOKE::SINGLE_PATH,
                        "Name of the folder containing the nominal kinetic scheme (XML Version)",
                        true,
                        "@NominalKineticsPreProcessor",
                        "none",
                        "none"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor",
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary containing the list of kinetic files to be interpreted",
                        true,
                        "@KineticsFolder",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NominalKineticsPreProcessor",
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary containing the list of nominal kinetic files to be interpreted",
                        true,
                        "@NominalKineticsFolder",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NameOfOptimizedKineticsFolder",
                        OpenSMOKE::SINGLE_PATH,
                        "Name of the folder where the optimized kinetics is printed at the end in CHEMKIN and XML format.",
                        false,
                        "none",
                        "none",
                        "none"));

                // To be removed        
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfOpenSMOKEInputFiles",
                        OpenSMOKE::VECTOR_STRING,
                        "List of OpenSMOKE input files",
                        true,
                        "@PathOpenSMOKEInputFiles @PathDatasetInputFiles",
                        "none",
                        "@PathOpenSMOKEInputFiles @PathDatasetInputFiles"));

                // To be removed
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfExperimentalDataFiles",
                        OpenSMOKE::VECTOR_STRING,
                        "List of Experimental data files",
                        true,
                        "@PathExperimentalDataFiles",
                        "none",
                        "@PathExperimentalDataFiles"));

                // To be removed
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PathOpenSMOKEInputFiles",
                        OpenSMOKE::SINGLE_PATH,
                        "File including the path for list of input files",
                        true,
                        "@ListOfOpenSMOKEInputFiles @PathDatasetInputFiles",
                        "none",
                        "@ListOfOpenSMOKEInputFiles @PathDatasetInputFiles"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PathDatasetInputFiles",
                        OpenSMOKE::SINGLE_PATH,
                        "File including the path for list of dataset specific input files",
                        true,
                        "@ListOfOpenSMOKEInputFiles @PathOpenSMOKEInputFiles",
                        "none",
                        "@ListOfOpenSMOKEInputFiles @PathOpenSMOKEInputFiles"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PathExperimentalDataFiles",
                        OpenSMOKE::SINGLE_PATH,
                        "File including the path for list of Experimental data files",
                        true,
                        "@ListOfExperimentalDataFiles",
                        "none",
                        "@ListOfExperimentalDataFiles"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OptimizationLibrary",
                        OpenSMOKE::SINGLE_STRING,
                        "Name of the library containing the optimization routines",
                        true,
                        "none",
                        "none",
                        "none"));
                
                // List of MAX ABS and MIN ABS and MAX REL and MIN REL vanno tenute e messe in optimization target
                // List initial li facciamo volare potenti e ce ne fottiamo
                // ICA e PCA gli dedicherei un dizionario a parte
                /*
                // ICA for direct reactions // //
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ICA_DirectReactionsIndices",
                        OpenSMOKE::VECTOR_INT,
                                "List of target reactions to be treated with ICA (indices starting from 1)",
                                false,
                        "none",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ICA_AbsMin_DirectReactions",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum values for ICA variables in direct reactions",
                        false,
                        "none",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ICA_AbsInit_DirectReactions",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum values for ICA variables in direct reactions",
                        false,
                        "none",
                        "none",
                        "none"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ICA_AbsMax_DirectReactions",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum values for ICA variables in direct reactions",
                        false,
                        "none",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ICA_files_MixingMatrices",
                        OpenSMOKE::SINGLE_STRING,
                        "File including the path for list of ICA mixing matrices files",
                        false,
                        "none",
                        "none",
                        "none"));
                // ICA for direct reactions // //

                // AB_  PCA for direct reactions // //
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_DirectReactions",
                        OpenSMOKE::VECTOR_INT,
                                "List of target reactions to be treated with PCA (indices starting from 1)",
                                false,
                        "none",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_MinAbs_DirectReactions",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum values for PCA variables in direct reactions",
                        false,
                        "none",
                        "none",
                        "none"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_MaxAbs_DirectReactions",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum values for PCA variables in direct reactions",
                        false,
                        "none",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_files_DirectReactions",
                        OpenSMOKE::SINGLE_STRING,
                        "File including the path for list of PCA input files",
                        false,
                        "none",
                        "none",
                        "none"));
                // AB_  PCA for direct reactions // //

                // AB_  PCA for P_inf reactions // //
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_LimitPressureReactions",
                        OpenSMOKE::VECTOR_INT,
                                "List of target reactions to be treated with PCA (indices starting from 1)",
                                false,
                        "none",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_MinAbs_LimitPressureReactions",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum values for PCA variables in direct reactions",
                        false,
                        "none",
                        "none",
                        "none"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_MaxAbs_LimitPressureReactions",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum values for PCA variables in direct reactions",
                        false,
                        "none",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_files_LimitPressureReactions",
                        OpenSMOKE::SINGLE_STRING,
                        "File including the path for list of PCA input files",
                        false,
                        "none",
                        "none",
                        "none"));
                // AB_  PCA for P_inf reactions // //

                // AB_  PCA for PLOG reactions // //
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_PLOGReactions",
                        OpenSMOKE::VECTOR_INT,
                                "List of target reactions to be treated with PCA (indices starting from 1)",
                                false,
                        "none",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_MinAbs_PLOGReactions",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum values for PCA variables in direct reactions",
                        false,
                        "none",
                        "none",
                        "none"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_MaxAbs_PLOGReactions",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum values for PCA variables in direct reactions",
                        false,
                        "none",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PCA_files_PLOGReactions",
                        OpenSMOKE::SINGLE_STRING,
                        "File including the path for list of PCA input files",
                        false,
                        "none",
                        "none",
                        "none"));
                // AB_  PCA for PLOG reactions // //

                // PARLARE CON AB SE ListOfInitial li facciamo volare o no.
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfInitial_A",
                        OpenSMOKE::VECTOR_STRING,
                        "List of initial values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_A",
                        "@ListOfInitial_lnA"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfInitial_lnA",
                        OpenSMOKE::VECTOR_STRING,
                        "List of initial values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_lnA",
                        "@ListOfInitial_A"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfInitial_Beta",
                        OpenSMOKE::VECTOR_STRING,
                        "List of initial values for temperature exponent",
                        false,
                        "none",
                        "@ListOfTarget_Beta",
                        "none"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfInitial_E",
                        OpenSMOKE::VECTOR_STRING,
                        "List of initial values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E",
                        "@ListOfInitial_E_over_R"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfInitial_E_over_R",
                        OpenSMOKE::VECTOR_STRING,
                        "List of initial values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_over_R",
                        "@ListOfInitial_E"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfInitial_ThirdBody_Eff",
                        OpenSMOKE::VECTOR_STRING,
                        "List of initial values for third body efficiencies",
                        false,
                        "none",
                        "@ListOfTarget_ThirdBody_Reactions",
                        "none"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_A",
                        OpenSMOKE::VECTOR_STRING,
                        "List of minimum absolute values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_A",
                        "@ListOfMinRel_A"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_A",
                        OpenSMOKE::VECTOR_STRING,
                        "List of maximum absolute values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_A",
                        "@ListOfMaxRel_A"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_A",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum relative values for frequency factors",
                        false,
                        "none",
                        "@ListOfTarget_A",
                        "@ListOfMinAbs_A"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_A",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum relative values for frequency factors",
                        false,
                        "none",
                        "@ListOfTarget_A",
                        "@ListOfMaxAbs_A"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_lnA",
                        OpenSMOKE::VECTOR_STRING,
                        "List of minimum absolute values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_lnA",
                        "@ListOfMinRel_lnA"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_lnA",
                        OpenSMOKE::VECTOR_STRING,
                        "List of maximum absolute values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_lnA",
                        "@ListOfMaxRel_lnA"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_lnA",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum relative values for frequency factors",
                        false,
                        "none",
                        "@ListOfTarget_lnA",
                        "@ListOfMinAbs_lnA"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_lnA",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum relative values for frequency factors",
                        false,
                        "none",
                        "@ListOfTarget_lnA",
                        "@ListOfMaxAbs_lnA"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_Beta",
                        OpenSMOKE::VECTOR_STRING,
                        "List of minimum absolute values for temperature exponent",
                        false,
                        "none",
                        "@ListOfTarget_Beta",
                        "@ListOfMinRel_Beta"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_Beta",
                        OpenSMOKE::VECTOR_STRING,
                        "List of maximum absolute values for temperature exponent",
                        false,
                        "none",
                        "@ListOfTarget_Beta",
                        "@ListOfMaxRel_Beta"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_Beta",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum relative values for temperature exponent",
                        false,
                        "none",
                        "@ListOfTarget_Beta",
                        "@ListOfMinAbs_Beta"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_Beta",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum relative values for temperature exponent",
                        false,
                        "none",
                        "@ListOfTarget_Beta",
                        "@ListOfMaxAbs_Beta"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_E",
                        OpenSMOKE::VECTOR_STRING,
                        "List of minimum absolute values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E",
                        "@ListOfMinRel_E"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_E",
                        OpenSMOKE::VECTOR_STRING,
                        "List of maximum absolute values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E",
                        "@ListOfMaxRel_E"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_E",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum relative values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E",
                        "@ListOfMinAbs_E"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_E",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum relative values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E",
                        "@ListOfMaxAbs_E"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_E_over_R",
                                OpenSMOKE::VECTOR_STRING,
                                "List of minimum absolute values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_over_R",
                        "@ListOfMinRel_E_over_R"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_E_over_R",
                                OpenSMOKE::VECTOR_STRING,
                                "List of maximum absolute values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_over_R",
                        "@ListOfMaxRel_E_over_R"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_E_over_R",
                                OpenSMOKE::VECTOR_DOUBLE,
                                "List of minimum relative values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_over_R",
                        "@ListOfMinAbs_E_over_R"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_E_over_R",
                                OpenSMOKE::VECTOR_DOUBLE,
                                "List of maximum relative values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_over_R",
                        "@ListOfMaxAbs_E_over_R"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfInitial_A_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of initial values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_A_inf",
                        "@ListOfInitial_lnA_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfInitial_lnA_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of initial values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_lnA_inf",
                        "@ListOfInitial_A_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfInitial_Beta_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of initial values for temperature exponent",
                        false,
                        "none",
                        "@ListOfTarget_Beta_inf",
                        "none"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfInitial_E_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of initial values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_inf",
                        "@ListOfInitial_E_over_R_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfInitial_E_over_R_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of initial values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_over_R_inf",
                        "@ListOfInitial_E_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_A_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of minimum absolute values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_A_inf",
                        "@ListOfMinRel_A_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_A_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of maximum absolute values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_A_inf",
                        "@ListOfMaxRel_A_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_A_inf",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum relative values for frequency factors",
                        false,
                        "none",
                        "@ListOfTarget_A_inf",
                        "@ListOfMinAbs_A_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_A_inf",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum relative values for frequency factors",
                        false,
                        "none",
                        "@ListOfTarget_A_inf",
                        "@ListOfMaxAbs_A_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_lnA_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of minimum absolute values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_lnA_inf",
                        "@ListOfMinRel_lnA_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_lnA_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of maximum absolute values for frequency factors (units: kmol, m, s)",
                        false,
                        "none",
                        "@ListOfTarget_lnA_inf",
                        "@ListOfMaxRel_lnA_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_lnA_inf",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum relative values for frequency factors",
                        false,
                        "none",
                        "@ListOfTarget_lnA_inf",
                        "@ListOfMinAbs_lnA_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_lnA_inf",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum relative values for frequency factors",
                        false,
                        "none",
                        "@ListOfTarget_lnA_inf",
                        "@ListOfMaxAbs_lnA_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_Beta_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of minimum absolute values for temperature exponent",
                        false,
                        "none",
                        "@ListOfTarget_Beta_inf",
                        "@ListOfMinRel_Beta_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_Beta_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of maximum absolute values for temperature exponent",
                        false,
                        "none",
                        "@ListOfTarget_Beta_inf",
                        "@ListOfMaxRel_Beta_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_Beta_inf",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum relative values for temperature exponent",
                        false,
                        "none",
                        "@ListOfTarget_Beta_inf",
                        "@ListOfMinAbs_Beta_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_Beta_inf",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum relative values for temperature exponent",
                        false,
                        "none",
                        "@ListOfTarget_Beta_inf",
                        "@ListOfMaxAbs_Beta_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_E_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of minimum absolute values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_inf",
                        "@ListOfMinRel_E_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_E_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of maximum absolute values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_inf",
                        "@ListOfMaxRel_E_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_E_inf",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum relative values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_inf",
                        "@ListOfMinAbs_E_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_E_inf",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum relative values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_inf",
                        "@ListOfMaxAbs_E_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_E_over_R_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of minimum absolute values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_over_R_inf",
                        "@ListOfMinRel_E_over_R_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_E_over_R_inf",
                        OpenSMOKE::VECTOR_STRING,
                        "List of maximum absolute values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_over_R_inf",
                        "@ListOfMaxRel_E_over_R_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_E_over_R_inf",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of minimum relative values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_over_R_inf",
                        "@ListOfMinAbs_E_over_R_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_E_over_R_inf",
                        OpenSMOKE::VECTOR_DOUBLE,
                        "List of maximum relative values for activation temperature",
                        false,
                        "none",
                        "@ListOfTarget_E_over_R_inf",
                        "@ListOfMaxAbs_E_over_R_inf"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinAbs_ThirdBody_Eff",
                                OpenSMOKE::VECTOR_STRING,
                                "List of minimum absolute values fo third body efficiency",
                                false,
                        "none",
                        "@ListOfTarget_ThirdBody_Reactions",
                        "@ListOfMinRel_ThirdBody_Eff"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxAbs_ThirdBody_Eff",
                                OpenSMOKE::VECTOR_STRING,
                                "List of maximum absolute values fo third body efficiency",
                                false,
                        "none",
                        "@ListOfTarget_ThirdBody_Reactions",
                        "@ListOfMaxRel_ThirdBody_Eff"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMinRel_ThirdBody_Eff",
                                OpenSMOKE::VECTOR_DOUBLE,
                                "List of minimum relative values fo third body efficiency",
                                false,
                        "none",
                        "@ListOfTarget_ThirdBody_Reactions",
                        "@ListOfMinAbs_ThirdBody_Eff"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfMaxRel_ThirdBody_Eff",
                                OpenSMOKE::VECTOR_DOUBLE,
                                "List of maximum relative values fo third body efficiency",
                                false,
                        "none",
                        "@ListOfTarget_ThirdBody_Reactions",
                        "@ListOfMaxAbs_ThirdBody_Eff"));
                */
	}
    };
}

#endif // GRAMMAR_OPTISMOKEPP_H
