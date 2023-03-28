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
|					  Timoteo Dinelli <timoteo.dinelli@polimi.it>	      |
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

namespace OptiSMOKE
{
    options_optimization_target::options_optimization_target() 
    {
        numberOfBatchReactor_ = 0;
        numberOfPlugFlowReactor_ = 0;
        numberOfPerfectlyStirredReactor_ = 0;
        numberOfPremixedLaminarFlame_ = 0;
        numberOfCounterFlowFlame_ = 0;

		numberOfParameters_ = 0;
    }
    
    options_optimization_target::~options_optimization_target() {}

    void options_optimization_target::SetupFromDictionary
                (OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                std::string dictionary_name)
    {
        dictionary_manager(dictionary_name).SetGrammar(optimization_target_grammar_);

        if (dictionary_manager(dictionary_name).CheckOption("@NumberOfBatchReactor"))
			dictionary_manager(dictionary_name).ReadInt("@NumberOfBatchReactor", 
													numberOfBatchReactor_);

        if (dictionary_manager(dictionary_name).CheckOption("@NumberOfPlugFlowReactor"))
		dictionary_manager(dictionary_name).ReadInt("@NumberOfPlugFlowReactor", 
													numberOfPlugFlowReactor_);

        if (dictionary_manager(dictionary_name).CheckOption("@NumberOfPerfectlyStirredReactor"))
			dictionary_manager(dictionary_name).ReadInt("@NumberOfPerfectlyStirredReactor", 
														numberOfPerfectlyStirredReactor_);

        if (dictionary_manager(dictionary_name).CheckOption("@NumberOfPremixedLaminarFlame"))
			dictionary_manager(dictionary_name).ReadInt("@NumberOfPremixedLaminarFlame", 
														numberOfPremixedLaminarFlame_);

        if(dictionary_manager(dictionary_name).CheckOption("@NumberOfCounterFlowFlame"))
            dictionary_manager(dictionary_name).ReadInt("@NumberOfCounterFlowFlame", 
														numberOfCounterFlowFlame_);

        // EPLR - Which reactions for Optimization of A, n, Ea
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_EPLR"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_EPLR", list_of_target_EPLR_);
		
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_BathGases_EPLR"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_BathGases_EPLR", 
															list_of_bath_gases_EPLR_);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_EPLR"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors_EPLR", 
															list_of_uncertainty_factors_EPLR_);

		numberOfParameters_ += list_of_target_EPLR_.size() * 3;
		
		// Extended PLOG - Which reactions for Optimization of A, n, Ea
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ExtPLOG_Reactions"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_ExtPLOG_Reactions", 
															list_of_target_extplog_);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_ExtPLOG"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors_ExtPLOG", 
															list_of_uncertainty_factors_extplog_);

		numberOfParameters_ += list_of_target_extplog_.size() * 3;
		
		// Extended PLOG - Which reactions and which Third Bodies to Optimize
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ExtPLOG_Reactions_TB"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_ExtPLOG_Reactions_TB", 
															list_of_target_extended_plog_reactions_);

		numberOfParameters_ += list_of_target_extended_plog_reactions_.size();

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ExtPLOG_Species"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_ExtPLOG_Species", 
															list_of_target_extended_plog_species_);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfMin_TBeff_ExtPLOG"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfMin_TBeff_ExtPLOG", 
															list_of_min_tb_extplog_);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfMax_TBeff_ExtPLOG"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfMax_TBeff_ExtPLOG", 
															list_of_max_tb_extplog_);

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
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_E_over_R_inf", 
															list_of_target_E_over_R_inf_);

		numberOfParameters_ += list_of_target_E_over_R_inf_.size();

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ThirdBody_Reactions"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_ThirdBody_Reactions", 
															list_of_target_thirdbody_reactions_);

		numberOfParameters_ += list_of_target_thirdbody_reactions_.size();
		
		// Classic PLOG - which
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_classic_PLOG_Reactions"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_classic_PLOG_Reactions", 
															list_of_target_classic_plog_reactions_);

		numberOfParameters_ += list_of_target_classic_plog_reactions_.size() * 3;

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_classic_PLOG"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors_classic_PLOG", 
															list_of_uncertainty_factors_classic_plog_);

		// List of third body species
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ThirdBody_Species"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_ThirdBody_Species", 
															list_of_target_thirdbody_species_);
		
		// List of target reactions for uncertainty factors
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTargetUncertaintyFactors"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTargetUncertaintyFactors", 
															list_of_target_uncertainty_factors_);
	
		// List of uncertainty factors
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors", 
															list_of_uncertainty_factors_);
	
		// List of target inf reactions for uncertainty factors
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTargetUncertaintyFactors_inf"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTargetUncertaintyFactors_inf", 
															list_of_target_uncertainty_factors_inf_);
	
		// List of inf uncertainty factors
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_inf"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors_inf", 
															list_of_uncertainty_factors_inf_);

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

		if(dictionary_manager(dictionary_name).CheckOption("@ReactionsClassesDefinitions"))
		{
            dictionary_manager(dictionary_name).ReadPath("@ReactionsClassesDefinitions", reactions_classes_definition_);
			if(!fs::exists(reactions_classes_definition_))
				OptiSMOKE::FatalErrorMessage("The @ReactionsClassesDefinitions path does not exists!");

			ReadReactionClassesDefinition(reactions_classes_definition_);
        }
    }

    void options_optimization_target::ReadReactionClassesDefinition(fs::path classes_definition)
	{
	}

} // namespace OptiSMOKE
