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
        
    private:

        grammar_optimization_targets optimization_target_grammar_;

        void ReadReactionClassesDefinition(fs::path reaction_classes);

        int NumberOfBatchReactor;
        int NumberOfPlugFlowReactor;
        int NumberOfPerfectlyStirredReactor;
        int NumberOfPremixedLaminarFlame;
        int NumberOfCounterFlowFlame;

        bool iReactionClasses; // Per ora uno non deve sbagliare ma sarebbe meglio mettere un check
                               // nel senso che in teoria è vero che uno puo ottimizzare sia classi 
                               // che altro assieme

        fs::path reactions_classes_definition_;
        
        struct reactions_
        {
            // EPLR
	        std::vector<int> list_of_target_EPLR;
	        std::vector<std::string> list_of_bath_gases_EPLR;
	        std::vector<double> list_of_uncertainty_factors_EPLR;
            	
            // Indices of extended plog for optimization
	        std::vector<int> list_of_target_extplog;
	        std::vector<double> list_of_uncertainty_factors_extplog;

            std::vector<int> list_of_target_extended_plog_reactions;
	        std::vector<std::string> list_of_target_extended_plog_species;
            
            // DIRECT REACTIONS
	        std::vector<int> list_of_target_lnA;
	        std::vector<int> list_of_target_Beta;
	        std::vector<int> list_of_target_E_over_R;

	        // INF REACTIONS
	        std::vector<int> list_of_target_lnA_inf;
	        std::vector<int> list_of_target_Beta_inf;
	        std::vector<int> list_of_target_E_over_R_inf;

	        // INF REACTIONS
	        std::vector<int> list_of_target_thirdbody_reactions;
	        std::vector<std::string> list_of_target_thirdbody_species;

	        // PLOG -- AB
	        std::vector<int> list_of_target_classic_plog_reactions;
	        // AB // List of uncertainties for plogs
	        std::vector<double> list_of_uncertainty_factors_classic_plog;

            std::vector<double> list_of_min_tb_extplog;
            std::vector<double> list_of_max_tb_extplog;
            
            std::vector<int> list_of_target_uncertainty_factors;
	        std::vector<double> list_of_uncertainty_factors;

	        std::vector<int> list_of_target_uncertainty_factors_inf;
	        std::vector<double> list_of_uncertainty_factors_inf;
        };

        // PENSIAMO ALLA CLASSE PIÙ GENERICA DELLA TERRA DOVREI
        // USARE LA STRUTTURA DI SOPRA ALLO STESSO MODO NO? 
        struct reactions_classes_
        {
            std::string class_name;
            std::vector<int> list_of_target_reactions;
            std::vector<double> list_of_uncertainty_factor;
        };

        std::vector<reactions_classes_> classes_structure_;
        reactions_ target_optimization_;

    public:

        options_optimization_target();
        ~options_optimization_target();

        void SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                                std::string dictionary_name);

        inline const int& number_of_batch_reactor() {return NumberOfBatchReactor;};
        inline const int& number_of_plug_flow_reactor() {return NumberOfPlugFlowReactor;};
        inline const int& number_of_perfectly_stirred_reactor() {return NumberOfPerfectlyStirredReactor;};
        inline const int& number_of_premixed_laminar_flame_reactor() {return NumberOfPremixedLaminarFlame;};
        inline const int& number_of_counter_flow_flame_reactor() {return NumberOfCounterFlowFlame;};

        inline const reactions_& optimization_target() const {return target_optimization_;};
        inline const std::vector<reactions_classes_>& classes_structure() const {return classes_structure_;};

    };
    
} // namespace OptiSMOKE

#include "options_optimization_target.hpp"
#endif // OPTIONS_OPTIMIZATION_TARGET_H