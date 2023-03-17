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

#ifndef OPTIONS_OPTIMIZATION_TARGET_H
#define OPTIONS_OPTIMIZATION_TARGET_H

namespace OptiSMOKE
{
    class options_optimization_target
    {
    
    public:

        options_optimization_target();
        ~options_optimization_target();

        void SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                                std::string dictionary_name);

        inline const int& number_of_batch_reactor() const {return numberOfBatchReactor_;};
        inline const int& number_of_plug_flow_reactor() const {return numberOfPlugFlowReactor_;};
        inline const int& number_of_perfectly_stirred_reactor() const {return numberOfPerfectlyStirredReactor_;};
        inline const int& number_of_premixed_laminar_flame() const {return numberOfPremixedLaminarFlame_;};
        inline const int& number_of_counter_flow_flame() const {return numberOfCounterFlowFlame_;};
        inline const int& number_of_parameters() const {return numberOfParameters_;};

	    inline const std::vector<int>& list_of_target_EPLR() const {return list_of_target_EPLR_;};
	    inline const std::vector<std::string>& list_of_bath_gases_EPLR() const {return list_of_bath_gases_EPLR_;};
	    inline const std::vector<double>& list_of_uncertainty_factors_EPLR() const {return list_of_uncertainty_factors_EPLR_;};
            	
	    inline const std::vector<int>& list_of_target_extplog() const {return list_of_target_extplog_;};
	    inline const std::vector<double>& list_of_uncertainty_factors_extplog() const {return list_of_uncertainty_factors_extplog_;};

        inline const std::vector<int>& list_of_target_extended_plog_reactions() const {return list_of_target_extended_plog_reactions_;};
	    inline const std::vector<std::string>& list_of_target_extended_plog_species() const {return list_of_target_extended_plog_species_;};
            
	    inline const std::vector<int>& list_of_target_lnA() const {return list_of_target_lnA_;};
	    inline const std::vector<int>& list_of_target_Beta() const {return list_of_target_Beta_;};
	    inline const std::vector<int>& list_of_target_E_over_R() const {return list_of_target_E_over_R_;};

	    inline const std::vector<int>& list_of_target_lnA_inf() const {return list_of_target_lnA_inf_;};
	    inline const std::vector<int>& list_of_target_Beta_inf() const {return list_of_target_Beta_inf_;};
	    inline const std::vector<int>& list_of_target_E_over_R_inf() const {return list_of_target_E_over_R_inf_;};

	    inline const std::vector<int>& list_of_target_thirdbody_reactions() const {return list_of_target_thirdbody_reactions_;};
	    inline const std::vector<std::string>& list_of_target_thirdbody_species() const {return list_of_target_thirdbody_species_;};

	    inline const std::vector<int>& list_of_target_classic_plog_reactions() const {return list_of_target_classic_plog_reactions_;};
	    inline const std::vector<double>& list_of_uncertainty_factors_classic_plog() const {return list_of_uncertainty_factors_classic_plog_;};

        inline const std::vector<double>& list_of_min_tb_extplog() const {return list_of_min_tb_extplog_;};
        inline const std::vector<double>& list_of_max_tb_extplog() const {return list_of_max_tb_extplog_;};
            
        inline const std::vector<int>& list_of_target_uncertainty_factors() const {return list_of_target_uncertainty_factors_;};
	    inline const std::vector<double>& list_of_uncertainty_factors() const {return list_of_uncertainty_factors_;};

	    inline const std::vector<int>& list_of_target_uncertainty_factors_inf() const {return list_of_target_uncertainty_factors_inf_;};
	    inline const std::vector<double>& list_of_uncertainty_factors_inf() const {return list_of_uncertainty_factors_inf_;};
        
        inline const std::vector<double>& list_of_min_rel_lnA() const {return list_of_min_rel_lnA_;};
	    inline const std::vector<double>& list_of_max_rel_lnA() const {return list_of_max_rel_lnA_;};

        inline const std::vector<double>& list_of_min_rel_Beta() const {return list_of_min_rel_Beta_;};
	    inline const std::vector<double>& list_of_max_rel_Beta() const {return list_of_max_rel_Beta_;};

        inline const std::vector<double>& list_of_min_rel_E_over_R() const {return list_of_min_rel_E_over_R_;};
	    inline const std::vector<double>& list_of_max_rel_E_over_R() const {return list_of_max_rel_E_over_R_;};

        inline const std::vector<double>& list_of_min_rel_lnA_inf() const {return list_of_min_rel_lnA_inf_;};
	    inline const std::vector<double>& list_of_max_rel_lnA_inf() const {return list_of_max_rel_lnA_inf_;};

        inline const std::vector<double>& list_of_min_rel_Beta_inf() const {return list_of_min_rel_Beta_inf_;};
	    inline const std::vector<double>& list_of_max_rel_Beta_inf() const {return list_of_max_rel_Beta_inf_;};

        inline const std::vector<double>& list_of_min_rel_E_over_R_inf() const {return list_of_min_rel_E_over_R_inf_;};
	    inline const std::vector<double>& list_of_max_rel_E_over_R_inf() const {return list_of_max_rel_E_over_R_inf_;};

        inline const std::vector<double>& list_of_min_rel_thirdbody_eff() const {return list_of_min_rel_thirdbody_eff_;};
	    inline const std::vector<double>& list_of_max_rel_thirdbody_eff() const {return list_of_max_rel_thirdbody_eff_;};

        inline const std::vector<std::string>& list_of_min_abs_thirdbody_eff() const {return list_of_min_abs_thirdbody_eff_;};
	    inline const std::vector<std::string>& list_of_max_abs_thirdbody_eff() const {return list_of_max_abs_thirdbody_eff_;};
        
        // PENSIAMO ALLA CLASSE PIÙ GENERICA DELLA TERRA DOVREI
        // USARE LA STRUTTURA DI SOPRA ALLO STESSO MODO NO? 
        struct reactions_classes_
        {
            std::string class_name;
            std::vector<int> list_of_target_reactions;
            std::vector<double> list_of_uncertainty_factor;
        };
    
    private:

        grammar_optimization_targets optimization_target_grammar_;

        void ReadReactionClassesDefinition(fs::path reaction_classes);

        int numberOfBatchReactor_;
        int numberOfPlugFlowReactor_;
        int numberOfPerfectlyStirredReactor_;
        int numberOfPremixedLaminarFlame_;
        int numberOfCounterFlowFlame_;
        int numberOfParameters_;

        bool iReactionClasses; // Per ora uno non deve sbagliare ma sarebbe meglio mettere un check
                               // nel senso che in teoria è vero che uno puo ottimizzare sia classi 
                               // che altro assieme

        fs::path reactions_classes_definition_;
        std::vector<reactions_classes_> classes_structure_;

        // EPLR
	    std::vector<int> list_of_target_EPLR_;
	    std::vector<std::string> list_of_bath_gases_EPLR_;
	    std::vector<double> list_of_uncertainty_factors_EPLR_;
            	
        // Indices of extended plog for optimization
	    std::vector<int> list_of_target_extplog_;
	    std::vector<double> list_of_uncertainty_factors_extplog_;

        std::vector<int> list_of_target_extended_plog_reactions_;
	    std::vector<std::string> list_of_target_extended_plog_species_;
            
        // DIRECT REACTIONS
	    std::vector<int> list_of_target_lnA_;
	    std::vector<int> list_of_target_Beta_;
	    std::vector<int> list_of_target_E_over_R_;

	    // INF REACTIONS
	    std::vector<int> list_of_target_lnA_inf_;
	    std::vector<int> list_of_target_Beta_inf_;
	    std::vector<int> list_of_target_E_over_R_inf_;

	    // INF REACTIONS
	    std::vector<int> list_of_target_thirdbody_reactions_;
	    std::vector<std::string> list_of_target_thirdbody_species_;

	    // PLOG -- AB
	    std::vector<int> list_of_target_classic_plog_reactions_;
	    // List of uncertainties for PLOGS -- AB
	    std::vector<double> list_of_uncertainty_factors_classic_plog_;

        std::vector<double> list_of_min_tb_extplog_;
        std::vector<double> list_of_max_tb_extplog_;
            
        std::vector<int> list_of_target_uncertainty_factors_;
	    std::vector<double> list_of_uncertainty_factors_;

	    std::vector<int> list_of_target_uncertainty_factors_inf_;
	    std::vector<double> list_of_uncertainty_factors_inf_;

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


    };
    
} // namespace OptiSMOKE

#include "options_optimization_target.hpp"
#endif // OPTIONS_OPTIMIZATION_TARGET_H