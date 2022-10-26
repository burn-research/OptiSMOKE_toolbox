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
        NumberOfBatchReactor = 0;
        NumberOfPlugFlowReactor = 0;
        NumberOfPerfectlyStirredReactor = 0;
        NumberOfPremixedLaminarFlame = 0;
        NumberOfCounterFlowFlame = 0;
    }
    
    options_optimization_target::~options_optimization_target() {}

    void options_optimization_target::SetupFromDictionary
                (OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                std::string dictionary_name)
    {
        dictionary_manager(dictionary_name).SetGrammar(optimization_target_grammar_);

        if (dictionary_manager(dictionary_name).CheckOption("@NumberOfBatchReactor"))
			dictionary_manager(dictionary_name).ReadInt("@NumberOfBatchReactor", NumberOfBatchReactor);

        if (dictionary_manager(dictionary_name).CheckOption("@NumberOfPlugFlowReactor"))
			dictionary_manager(dictionary_name).ReadInt("@NumberOfPlugFlowReactor", NumberOfPlugFlowReactor);

        if (dictionary_manager(dictionary_name).CheckOption("@NumberOfPerfectlyStirredReactor"))
			dictionary_manager(dictionary_name).ReadInt("@NumberOfPerfectlyStirredReactor", 
														NumberOfPerfectlyStirredReactor);

        if (dictionary_manager(dictionary_name).CheckOption("@NumberOfPremixedLaminarFlame"))
			dictionary_manager(dictionary_name).ReadInt("@NumberOfPremixedLaminarFlame", 
														NumberOfPremixedLaminarFlame);

        if(dictionary_manager(dictionary_name).CheckOption("@NumberOfCounterFlowFlame"))
            dictionary_manager(dictionary_name).ReadInt("@NumberOfCounterFlowFlame", NumberOfCounterFlowFlame);

        // EPLR - Which reactions for Optimization of A, n, Ea
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_EPLR"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_EPLR", target_optimization_.list_of_target_EPLR);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_BathGases_EPLR"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_BathGases_EPLR", 
															target_optimization_.list_of_bath_gases_EPLR);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_EPLR"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors_EPLR", 
															target_optimization_.list_of_uncertainty_factors_EPLR);
		// EPLR - Which reactions for Optimization of A, n, Ea
		
		// Extended PLOG - Which reactions for Optimization of A, n, Ea
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ExtPLOG_Reactions"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_ExtPLOG_Reactions", 
															target_optimization_.list_of_target_extplog);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_ExtPLOG"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors_ExtPLOG", 
															target_optimization_.list_of_uncertainty_factors_extplog);

		// Extended PLOG - Which reactions and which Third Bodies to Optimize
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ExtPLOG_Reactions_TB"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_ExtPLOG_Reactions_TB", 
															target_optimization_.list_of_target_extended_plog_reactions);
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ExtPLOG_Species"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_ExtPLOG_Species", 
															target_optimization_.list_of_target_extended_plog_species);
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfMin_TBeff_ExtPLOG"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfMin_TBeff_ExtPLOG", 
															target_optimization_.list_of_min_tb_extplog);
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfMax_TBeff_ExtPLOG"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfMax_TBeff_ExtPLOG", 
															target_optimization_.list_of_max_tb_extplog);

		// Direct reactions - specific parameters to be optimized	
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_lnA"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_lnA", target_optimization_.list_of_target_lnA);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_Beta"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_Beta", target_optimization_.list_of_target_Beta);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_E_over_R"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_E_over_R", target_optimization_.list_of_target_E_over_R);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_lnA_inf"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_lnA_inf", target_optimization_.list_of_target_lnA_inf);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_Beta_inf"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_Beta_inf", target_optimization_.list_of_target_Beta_inf);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_E_over_R_inf"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_E_over_R_inf", 
															target_optimization_.list_of_target_E_over_R_inf);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ThirdBody_Reactions"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_ThirdBody_Reactions", 
															target_optimization_.list_of_target_thirdbody_reactions);
		
		// Classic PLOG - which
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_classic_PLOG_Reactions"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_classic_PLOG_Reactions", 
															target_optimization_.list_of_target_classic_plog_reactions);

		if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_classic_PLOG"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors_classic_PLOG", 
															target_optimization_.list_of_uncertainty_factors_classic_plog);

		// List of third body species
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTarget_ThirdBody_Species"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTarget_ThirdBody_Species", 
															target_optimization_.list_of_target_thirdbody_species);
	
		// List of target reactions for uncertainty factors
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTargetUncertaintyFactors"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTargetUncertaintyFactors", 
															target_optimization_.list_of_target_uncertainty_factors);
	
		// List of uncertainty factors
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors", 
															target_optimization_.list_of_uncertainty_factors);
	
		// List of target inf reactions for uncertainty factors
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfTargetUncertaintyFactors_inf"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfTargetUncertaintyFactors_inf", 
															target_optimization_.list_of_target_uncertainty_factors_inf);
	
		// List of inf uncertainty factors
		if (dictionary_manager(dictionary_name).CheckOption("@ListOfUncertaintyFactors_inf"))
			dictionary_manager(dictionary_name).ReadOption("@ListOfUncertaintyFactors_inf", 
															target_optimization_.list_of_uncertainty_factors_inf);

        if(dictionary_manager(dictionary_name).CheckOption("@ReactionsClassesDefinitions"))
        {
            dictionary_manager(dictionary_name).ReadPath("@ReactionsClassesDefinitions", reactions_classes_definition_);
			ReadReactionClassesDefinition(reactions_classes_definition_);
        }
    }

    void options_optimization_target::ReadReactionClassesDefinition(fs::path classes_definition)
	{
		fs::ifstream fileHandler(classes_definition.c_str());
    	std::string line;
    	int numberOfLines = 0;
    	std::vector<std::string> content;

   		while (getline(fileHandler, line)) {
        	numberOfLines++;
        	content.push_back(line);
    	}

    	int NumberOfReactionClasses = numberOfLines / 5;

		classes_structure_.resize(NumberOfReactionClasses);

    	std::vector<std::string> reactionClassName;
    	std::vector<std::string> str_reaction_index; // tmp vector dove mi salvo le stringhe 
    	std::vector<std::string> str_unc; // tmp vector dove mi salvo le incertezze
    	std::vector <std::string> str_QOI; // tmp vector dove mi salvo le stringhe delle QOI

    	/*
        	Devo prendere content e:
            	- La prima riga di ogni blocco corrisponde al nome della classe
            	- la seconda riga agli indici delle reazione
            	- la terza riga bo scaling factor
            	- la quarta riga uncertainty
            	- la quinta riga quello che voglio ottimizzare ovvero che parametri dell'arrehnius 
			Per ora solo le dirette dopo faccio la funzione per scegliere il tipo
    	*/

   		// Con questo ciclo for mi sono separato le cose e ora ci posso lavorare
   		// però ho gia tutto quello che mi serve!
    	for(int i = 0; i<NumberOfReactionClasses; i++){
			classes_structure_[i].class_name = content[0 + i * 5];
        	reactionClassName.push_back(content[0 + i * 5]);
        	str_reaction_index.push_back(content[1 + i * 5]);
        	str_unc.push_back(content[3 + i * 5]);
        	str_QOI.push_back(content[4 + i * 5]);
    	}
    
    	for(int i = 0; i < NumberOfReactionClasses; i++){
        	std::vector<std::string> tmp_idx; // reaction index
        	std::vector<int> tmp_int_idx;
        	std::vector<std::string> tmp_unc; // uncertainty factors
        	std::vector<double> tmp_double_unc;
        	std::vector<std::string> tmp_QOI; // Quantity of interest

        	boost::split(tmp_idx, str_reaction_index[i], boost::is_any_of(" "));
        	boost::split(tmp_unc, str_unc[i], boost::is_any_of(" "));
        	boost::split(tmp_QOI, str_QOI[i], boost::is_any_of(" "));
        
        	for(int j = 0; j < tmp_idx.size(); j++){
            	tmp_int_idx.push_back(std::stoi(tmp_idx[j]));
        	}
        	for(int k = 0; k < tmp_unc.size(); k++){
            	tmp_double_unc.push_back(std::stod(tmp_unc[k]));
       		}
        	// matrixOfReactionIndex.push_back(tmp_int_idx);
			classes_structure_[i].list_of_target_reactions = tmp_int_idx;
        	// matrixOfUnceratintyFactors.push_back(tmp_double_unc);
			classes_structure_[i].list_of_uncertainty_factor = tmp_double_unc;
        	// matrixOfQOI.push_back(tmp_QOI);

        	tmp_idx.clear();
        	tmp_int_idx.clear();
        	tmp_unc.clear();
        	tmp_double_unc.clear();
        	tmp_QOI.clear();
   		}
		
    	for(int i = 0; i<NumberOfReactionClasses; i++)
		{
			target_optimization_.list_of_target_lnA.push_back(classes_structure_[i].list_of_target_reactions[0]);
			target_optimization_.list_of_target_Beta.push_back(classes_structure_[i].list_of_target_reactions[0]);
			target_optimization_.list_of_target_E_over_R.push_back(classes_structure_[i].list_of_target_reactions[0]);
			target_optimization_.list_of_target_uncertainty_factors.push_back(classes_structure_[i].list_of_target_reactions[0]);

			// list_of_initial_lnA.push_back(boost::lexical_cast<std::string>(std::log(kineticsMapXML->A(matrixOfReactionIndex[i][0]-1))));
			// list_of_initial_Beta.push_back(boost::lexical_cast<std::string>(kineticsMapXML->Beta(matrixOfReactionIndex[i][0]-1)));
			// list_of_initial_E_over_R.push_back(boost::lexical_cast<std::string>(kineticsMapXML->E_over_R(matrixOfReactionIndex[i][0]-1)));
        	
		}
		
		// Da capire se sta roba ha senso
		for(int i = 0; i<NumberOfReactionClasses; i++)
		{
			target_optimization_.list_of_uncertainty_factors.push_back(classes_structure_[i].list_of_uncertainty_factor[0]);
		}

		std::cout << "Reading the reaction classes definition in: " << classes_definition.c_str() << std::endl;
    	std::cout << " * Number of  reaction classess defined: " << numberOfLines/5 << std::endl;
		for(int i = 0; i < NumberOfReactionClasses; i++)
		{
			std::cout << " * Reaction class: " << classes_structure_[i].class_name << std::endl;
			for(int j = 0; j < classes_structure_[i].list_of_target_reactions.size(); j++)
			{
				std::cout << "   " <<  classes_structure_[i].list_of_target_reactions[j] 
				<< "   " << classes_structure_[i].list_of_uncertainty_factor[j] << std::endl;
			}
		}
	}

} // namespace OptiSMOKE
