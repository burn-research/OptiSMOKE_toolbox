namespace OptiSMOKE{

	SimulationsInterface::SimulationsInterface(const OptiSMOKE::InputManager& data) : data_(data)
	{
		// Resize reactors objects
		
		n_batch = data_.optimization_target().number_of_batch_reactor();
		n_pfr = data_.optimization_target().number_of_plug_flow_reactor();
		n_psr = data_.optimization_target().number_of_perfectly_stirred_reactor();
		n_premixed = data_.optimization_target().number_of_premixed_laminar_flame();
		n_counterflow = data_.optimization_target().number_of_counter_flow_flame();

		// Resize simulations results
		simulations_results_.resize(data_.expdata_x().size());
		for(unsigned int i = 0; i < data_.expdata_x().size(); i++){
			simulations_results_[i].resize(data_.expdata_x()[i].size());
			for(unsigned int j = 0; j < data_.expdata_x()[i].size(); j++){
				simulations_results_[i][j].resize(data_.expdata_x()[i][j].size());
			}
		}
	}

	SimulationsInterface::~SimulationsInterface(){}

	void SimulationsInterface::Setup()
	{

		//-------------------------------------------------------//
		//               Allocations of reactors object          //
		//-------------------------------------------------------//
		unsigned int offset = 0;
		if(n_batch != 0){
			batch_reactors.resize(n_batch);
			for(unsigned int i = 0; i < n_batch; i++)
				batch_reactors[i] = new OptiSMOKE::BatchReactor[data_.input_paths()[i].size()];
			
			offset += n_batch;
		}

		if(n_pfr != 0){
			plugflow_reactors.resize(n_pfr);
			for(unsigned int i = 0; i < n_pfr; i++)
				plugflow_reactors[i] = new OptiSMOKE::PlugFlowReactor[data_.input_paths()[i+offset].size()];
			
			offset += n_pfr; 
		}

		//-------------------------------------------------------//
		//               Reactions constraints                   //
		//-------------------------------------------------------//

		// Create constraints for reactions
		// Initializes vector of uncertainty factors for Direct reactions
		// Compute max and min reaction rates according to temperature dependent 
		// or constant uncertainty factors for direct reactions
		std::vector<double> f_factors = data_.optimization_target().list_of_uncertainty_factors();

		// Computing boundaries for all the reactions
		int rows = data_.optimization_target().list_of_target_uncertainty_factors().size();
		int columns = T_span.size();
		std::vector<std::vector<double>> k_0(rows, std::vector<double>(columns));
		k_upper.resize(rows, std::vector<double>(columns));
		k_lower.resize(rows, std::vector<double>(columns));

		// if (ObjectInput2.udc_bool == true){
		// Not yet added
		// } 
		// else{
		for (int j = 0; j < rows; j++){
			int reaction_index = data_.optimization_target().list_of_target_uncertainty_factors()[j];
			double A_0 = data_.nominalkineticsMapXML_->A(reaction_index-1);
			double Beta_0 = data_.nominalkineticsMapXML_->Beta(reaction_index-1);
			double E_over_R_0 = data_.nominalkineticsMapXML_->E_over_R(reaction_index-1);
			for (int i=0; i < columns; i++){
				k_0[j][i] = A_0 * std::pow(T_span[i], Beta_0) * std::exp((-1*E_over_R_0)/T_span[i]);
				k_upper[j][i] = k_0[j][i] * std::pow(10, f_factors[j]) * ((double)data_.optimization_setup().sigma_k_distribution()/2);
				k_lower[j][i] = k_0[j][i] * std::pow(10,-1*f_factors[j]) / ((double)data_.optimization_setup().sigma_k_distribution()/2);
			}

			// std::cout << "Reaction " << reaction_index << ": " << std::endl;
			// for (int i=0; i < T_span.size(); i++) 
            // {    
            //     std::cout<< "The value of f for the lower constraint is: " << k_lower[j][i] <<std::endl;
            // }    

            // for (int i=0; i < T_span.size(); i++) 
            // {    
            //     std::cout<< "The value of f for the upper constraint is: " << k_upper[j][i] <<std::endl;
            // }   
		}
		// }

		// Compute max and min reaction rates according to temperature dependent 
		// or constant uncertainty factors for fall off reactions
		std::vector<double> f_factors_inf = data_.optimization_target().list_of_uncertainty_factors_inf();
		int rows_inf = data_.optimization_target().list_of_target_uncertainty_factors_inf().size();
		int columns_inf = T_span.size();
		std::vector<unsigned int> indices_of_falloff_reactions = data_.nominalkineticsMapXML_->IndicesOfFalloffReactions();
		std::vector<std::vector<double>> k_0_inf(rows_inf, std::vector<double>(columns_inf));
		k_upper_inf.resize(rows_inf, std::vector<double>(columns_inf));
		k_lower_inf.resize(rows_inf, std::vector<double>(columns_inf));

		for (int j=0; j < rows_inf; j++){
			int falloff_reaction_index = data_.optimization_target().list_of_target_uncertainty_factors_inf()[j];
			int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(), indices_of_falloff_reactions.end(), falloff_reaction_index) - indices_of_falloff_reactions.begin();
			double A_0_inf = data_.nominalkineticsMapXML()->A_falloff_inf(pos_FallOff_Reaction);
			double Beta_0_inf = data_.nominalkineticsMapXML()->Beta_falloff_inf(pos_FallOff_Reaction);
			double E_over_R_0_inf = data_.nominalkineticsMapXML()->E_over_R_falloff_inf(pos_FallOff_Reaction);
			for (int i=0; i < columns_inf; i++){
				k_0_inf[j][i] = A_0_inf * std::pow(T_span[i],Beta_0_inf) * std::exp((-1*E_over_R_0_inf)/T_span[i]);
				k_upper_inf[j][i] = k_0_inf[j][i] * std::pow(10,f_factors_inf[j])* (double)data_.optimization_setup().sigma_k_distribution()/2;
				k_lower_inf[j][i] = k_0_inf[j][i] * std::pow(10,-1*f_factors_inf[j])/ ((double)data_.optimization_setup().sigma_k_distribution()/2);
			}
		}

		// Compute max and min reaction rates according to temperature dependent 
		// or constant uncertainty factors for for PLOG reactions
		std::vector<std::vector<std::vector<double>>> k_classic_plog;
		std::vector<unsigned int> indices_of_classic_plog = data_.nominalkineticsMapXML()->IndicesOfPLOGReactions();
		int rows_cp = data_.optimization_target().list_of_uncertainty_factors_classic_plog().size();
		//resize first dimensions according to the number of classic plog we are interested in for optimisation
    	k_classic_plog.resize(rows_cp);		
		k_upper_classic_plog.resize(rows_cp);
		k_lower_classic_plog.resize(rows_cp);
	
		for (unsigned int j=0; j < rows_cp; j++){

			//finding the position of the target reaction in correspondent OS++ object
			int pos_classic_plog_reaction = std::find(
				indices_of_classic_plog.begin(), 
				indices_of_classic_plog.end(), 
				data_.optimization_target().list_of_target_classic_plog_reactions()[j]
			)-indices_of_classic_plog.begin();
			
			for(int k=0; k<data_.nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA().size(); k++){
				if(data_.nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA()[k].size() != 1){
					std::cout << " * WARNING the plog reactions requested for the optimization is a duplicate!" << std::endl;
					std::cout << "   * Reaction: " << data_.optimization_target().list_of_target_classic_plog_reactions()[j] << std::endl;
				}
			}	
			
			k_classic_plog[j].resize(data_.nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA().size());
			k_upper_classic_plog[j].resize(data_.nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA().size());
			k_lower_classic_plog[j].resize(data_.nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA().size());

			for (int k=0; k < data_.nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA().size(); k++){
				double A_CP = std::exp(data_.nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA()[k][0]);
				double n_CP = data_.nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).Beta()[k][0];
				double E_over_R_CP = data_.nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).E_over_R()[k][0];

				k_classic_plog[j][k].resize(T_span.size());
				k_upper_classic_plog[j][k].resize(T_span.size());
				k_lower_classic_plog[j][k].resize(T_span.size());

				for (int i=0; i < T_span.size(); i++){
					k_classic_plog[j][k][i] = A_CP * std::pow(T_span[i],n_CP) * std::exp((-1*E_over_R_CP)/T_span[i]);
					k_upper_classic_plog[j][k][i] = k_classic_plog[j][k][i] * std::pow(10, data_.optimization_target().list_of_uncertainty_factors_classic_plog()[j]) * ((double)data_.optimization_setup().sigma_k_distribution()/2);
					k_lower_classic_plog[j][k][i] = k_classic_plog[j][k][i] * std::pow(10, -1*data_.optimization_target().list_of_uncertainty_factors_classic_plog()[j])/ ((double)data_.optimization_setup().sigma_k_distribution()/2);
				}
			}		
		}
	}

	void SimulationsInterface::run()
	{
		OpenSMOKE::KineticsMap_CHEMKIN* kinetics = data_.kineticsMapXML_;
		OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermo = data_.thermodynamicsMapXML_;

		// Loop over all datasets
		// Here data_.path_experimental_data_files().size() this is 
		// misleading however keep in mind that only the size matters

		for(unsigned int i = 0; i < data_.path_experimental_data_files().size(); i++)
		{
			std::string qoi = data_.QoI()[i];
			std::string qoi_target = data_.QoI_target()[i];
			std::string solver = data_.solver_name()[i];

			if (solver == "BatchReactor"){
				if(qoi == "IDT"){
					for(unsigned int j = 0; j < data_.input_paths()[i].size(); j++){
						// For the moment stay simple but keep in mind that this every time 
						// re-read the OS input file
						batch_reactors[i][j].Setup(data_.input_paths()[i][j], thermo, kinetics);
						batch_reactors[i][j].Solve();
						simulations_results_[i][0][j] = batch_reactors[i][j].GetIgnitionDelayTime(qoi_target) * std::pow(10, 6);
					}
				}
				else if (qoi == "Composition"){
					OptiSMOKE::FatalErrorMessage("Compositions profile measurements in batch reactors not yet implemented!");
				}
				else{
					OptiSMOKE::FatalErrorMessage("Unknown QoI: " + qoi);
				}
			}

			if(solver == "PlugFlowReactor"){
				if(qoi == "Composition"){
					for(unsigned int j = 0; j < data_.input_paths()[i].size(); j++){
						plugflow_reactors[i][j].Setup(data_.input_paths()[i][j], thermo, kinetics);
						plugflow_reactors[i][j].Solve();
						for(unsigned int k = 0; k < species_list; k++){
							simulations_results_[i][k][j] = plugflow_reactors[i][j].GetMolefraction();
						}
					}
				}
				else{
					OptiSMOKE::FatalErrorMessage("Unknown QoI: " + qoi);
				}
			}

			if(solver == "PerfectlyStirredReactor"){
				OptiSMOKE::FatalErrorMessage(solver + " not supported yet!");
			}
			
			if (solver == "PremixedLaminarFlame1D"){
				OptiSMOKE::FatalErrorMessage(solver + " not supported yet!");
			}
            
			if (solver == "CounterFlowFlame1D"){
				OptiSMOKE::FatalErrorMessage(solver + " not supported yet!");
			}
		}
	}

	double SimulationsInterface::ComputeObjectiveFunction()
	{
		double objective_function = 0;

		if (data_.optimization_setup().objective_function_type() == "CurveMatching") {
			std::vector<double> CM_indexes;
			std::vector<std::vector<double>> CM_score;
			CM_score.resize(simulations_results_.size());
			
			double tStart_CM = OpenSMOKE::OpenSMOKEGetCpuTime();
			for(unsigned int i = 0; i < simulations_results_.size(); i++){
				CM_score[i].resize(simulations_results_[i].size());
				for(unsigned int j = 0; j < simulations_results_[i].size(); j++){
					if(data_.QoI()[i] == "IDT"){
						std::vector<double> tmp;
						std::vector<double> tmp2;
						for(int k = 0; k < data_.expdata_y()[i][j].size(); k++){
							tmp.push_back(std::log(data_.expdata_y()[i][j][k]));
							tmp2.push_back(std::log(simulations_results_[i][j][k]));
						}
						CM_indexes = curveMatching(data_.curvematching_options().number_of_bootstrap(),
									data_.expdata_x()[i][j],
									tmp,
									data_.expdata_x()[i][j],
									tmp2,
									data_.uncertainty()[i][j]);
					}
					else{
						CM_indexes = curveMatching(data_.curvematching_options().number_of_bootstrap(),
									data_.expdata_x()[i][j],
									data_.expdata_y()[i][j],
									data_.expdata_x()[i][j],
									simulations_results_[i][j],
									data_.uncertainty()[i][j]);
					}
					CM_score[i][j] = (CM_indexes[0] + CM_indexes[1] + CM_indexes[2] + CM_indexes[3])/4;
				}	
			}
			
			double final_index = 0; 
            for (int i =0; i < CM_score.size(); i++){    
                double i_th_index = 0; 
                for (int j=0; j < CM_score[i].size(); j++){
                	i_th_index = i_th_index  + CM_score[i][j]/CM_score[i].size(); 
                }
                std::cout << "    * The Curve Matching score of ";
				std::cout << data_.dataset_names()[i] << " is: " << i_th_index << std::endl; 
                final_index = final_index+i_th_index/CM_score.size();
            }

			double tEnd_CM = OpenSMOKE::OpenSMOKEGetCpuTime();
            std::cout << " * Time to compute the Curve Matching using ";
			std::cout << data_.curvematching_options().number_of_bootstrap();
			std::cout << " bootstrap variations: " << tEnd_CM - tStart_CM << std::endl;
			std::cout << " * The final Curve Matching score is ";
			std::cout << final_index << std::endl;
            objective_function = 1 - final_index;
		}
		else if (data_.optimization_setup().objective_function_type() == "L1-norm"){
			OptiSMOKE::FatalErrorMessage("L1-norm not yet implemented!");
		} 
		else if (data_.optimization_setup().objective_function_type() == "L2-norm"){
			OptiSMOKE::FatalErrorMessage("L2-norm not yet implemented!");
		}
		else {
			OptiSMOKE::FatalErrorMessage("Unknown type of objective function. Available are: L1-norm | L2-norm | CurveMatching");		
		}
		return objective_function;	
	}

	bool SimulationsInterface::CheckKineticConstasts()
	{
		// OptiSMOKE mantainer of the future in order to make things more effcient consider
		// to not substituting kinetics and then creating constraints and then checking 
		// just check and substitute or not
		for (int j=0; j < data_.optimization_target().list_of_target_uncertainty_factors().size(); j++){
			std::vector<double> k_check;
			k_check.resize(T_span.size());
			for (int i=0; i < T_span.size(); i++){
				// Calculating the k value for the new set of kinetic parameters for a temperature span of 300-3000 K
				unsigned int reaction_index = data_.optimization_target().list_of_target_uncertainty_factors()[j]-1;
				double A_j = data_.kineticsMapXML()->A(reaction_index);
				double Beta_j = data_.kineticsMapXML()->Beta(reaction_index);
				double E_over_R_j = data_.kineticsMapXML()->E_over_R(reaction_index);

				k_check[i] = A_j * std::pow(T_span[i], Beta_j) * std::exp((-1 * E_over_R_j)/T_span[i]);

				// If at one temperature the k value is either below the lower bound or above the upper bound, 
				// forcefully set the objective function value to 1e7 and print out for which reaction the violation occured 
				if ((k_check[i]<=k_lower[j][i]) || (k_check[i]>=k_upper[j][i])){
					std::cout << " * Violation for reaction: ";
					std::cout << reaction_index + 1 << std::endl;
					return true;
				}
			}
		}

		std::vector<unsigned int> indices_of_falloff_reactions = data_.kineticsMapXML()->IndicesOfFalloffReactions();
		for (int j=0; j < data_.optimization_target().list_of_target_uncertainty_factors_inf().size(); j++){
			std::vector<double> k_check_inf;
			k_check_inf.resize(T_span.size());
			
			// Finding position of the fall off reaction
			int pos_FallOff_Reaction = std::find(
				indices_of_falloff_reactions.begin(),
				indices_of_falloff_reactions.end(),
				data_.optimization_target().list_of_target_uncertainty_factors_inf()[j]
			) - indices_of_falloff_reactions.begin();
			
			for (int i=0; i < T_span.size(); i++){
				double A_falloff_inf_j = data_.kineticsMapXML()->A_falloff_inf(pos_FallOff_Reaction);
				double Beta_falloff_inf_j = data_.kineticsMapXML()->Beta_falloff_inf(pos_FallOff_Reaction);
				double E_over_R_falloff_inf_j = data_.kineticsMapXML()->E_over_R_falloff_inf(pos_FallOff_Reaction);
				k_check_inf[i] = A_falloff_inf_j * std::pow(T_span[i], Beta_falloff_inf_j) * std::exp((-1 * E_over_R_falloff_inf_j)/T_span[i]);

				if ((k_check_inf[i] <= k_lower_inf[j][i]) || (k_check_inf[i] >= k_upper_inf[j][i])){
					std::cout << " * Violation for reaction: ";
					std::cout << data_.optimization_target().list_of_target_uncertainty_factors_inf()[j];
					std::cout << " (inf) " << std::endl;
					return true;
				}
			}
		}

		std::vector<unsigned int> indices_of_classic_plog = data_.nominalkineticsMapXML()->IndicesOfPLOGReactions();
		for (int j=0; j < data_.optimization_target().list_of_uncertainty_factors_classic_plog().size(); j++){

			int pos_classic_plog_reaction = std::find(
				indices_of_classic_plog.begin(),
				indices_of_classic_plog.end(),
				data_.optimization_target().list_of_target_classic_plog_reactions()[j]
			) - indices_of_classic_plog.begin();

			for (int k=0; k < data_.kineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA().size(); k++){
				double A_CP_trial = std::exp(data_.kineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA()[k][0]);
				double n_CP_trial = data_.kineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).Beta()[k][0];
				double E_over_R_CP_trial = data_.kineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).E_over_R()[k][0];

				std::vector<double> k_check_CP;
				k_check_CP.resize(T_span.size());

				for (int i=0; i < T_span.size(); i++){
					k_check_CP[i] = A_CP_trial * std::pow(T_span[i], n_CP_trial) * std::exp((-1*E_over_R_CP_trial)/T_span[i]);

					if ((k_check_CP[i] <= k_lower_classic_plog[j][k][i]) || (k_check_CP[i] >= k_upper_classic_plog[j][k][i])){
						std::cout << " * Violation for PLOG reaction: ";
						std::cout << data_.optimization_target().list_of_target_classic_plog_reactions()[j] << std::endl;
						return true;
					}
				}
			}
		}
		return false;
	}

	void SimulationsInterface::SubstituteKineticParameters(const Dakota::RealVector& c_vars){
		unsigned int count = 0;

		// lnA
		if(data_.optimization_target().list_of_target_lnA().size() != 0){
			for(unsigned int i = 0; i < data_.optimization_target().list_of_target_lnA().size(); i++){
				ChangeDirectParamaters("lnA", data_.optimization_target().list_of_target_lnA()[i], c_vars[count]);
				count += 1;
			}
		}

		// lnA_inf
		if(data_.optimization_target().list_of_target_lnA_inf().size() != 0){
			for(unsigned int i = 0; i < data_.optimization_target().list_of_target_lnA_inf().size(); i++){
				ChangeFallOffParamaters("lnA", data_.optimization_target().list_of_target_lnA_inf()[i], c_vars[count]);
				count += 1;
			}
		}
		
		// Beta
		if(data_.optimization_target().list_of_target_Beta().size() != 0){
			for(unsigned int i = 0; i < data_.optimization_target().list_of_target_Beta().size(); i++){
				ChangeDirectParamaters("Beta", data_.optimization_target().list_of_target_Beta()[i], c_vars[count]);
				count += 1;
			}
		}
		
		// Beta_inf
		if(data_.optimization_target().list_of_target_Beta_inf().size() != 0){
			for(unsigned int i = 0; i < data_.optimization_target().list_of_target_Beta_inf().size(); i++){
				ChangeFallOffParamaters("Beta", data_.optimization_target().list_of_target_Beta_inf()[i], c_vars[count]);
				count += 1;
			}
		}
		
		// E_over_R
		if(data_.optimization_target().list_of_target_E_over_R().size() != 0){
			for(unsigned int i = 0; i < data_.optimization_target().list_of_target_E_over_R().size(); i++){
				ChangeDirectParamaters("E_over_R", data_.optimization_target().list_of_target_E_over_R()[i], c_vars[count]);
				count += 1;
			}
		}
		
		// E_over_R_inf
		if(data_.optimization_target().list_of_target_E_over_R_inf().size() != 0){
			for(unsigned int i = 0; i < data_.optimization_target().list_of_target_E_over_R_inf().size(); i++){
				ChangeFallOffParamaters("E_over_R", data_.optimization_target().list_of_target_E_over_R_inf()[i], c_vars[count]);
				count += 1;
			}			
		}
		
		// 3B eff
		if(data_.optimization_target().list_of_target_thirdbody_reactions().size() != 0){
			for (unsigned int i = 0; i <  data_.optimization_target().list_of_target_thirdbody_reactions().size(); i++){
                ChangeThirdBodyEfficiencies(data_.optimization_target().list_of_target_thirdbody_reactions()[i], 
											data_.optimization_target().list_of_target_thirdbody_species()[i], 
											c_vars[count]);
                count += 1;
        	}	
		}
		
		// Classic PLOG
		if(data_.optimization_target().list_of_target_classic_plog_reactions().size() != 0){
			// TODO
		}
	}

	void SimulationsInterface::ChangeDirectParamaters(std::string type, int index, double parameter){
		
		if(type == "lnA")
			data_.kineticsMapXML_->Set_A(index-1, std::exp(parameter));
		if(type == "Beta")
			data_.kineticsMapXML_->Set_Beta(index-1, parameter);
		if(type == "E_over_R")
			data_.kineticsMapXML_->Set_E_over_R(index-1, parameter);
	}

	void SimulationsInterface::ChangeFallOffParamaters(std::string type, int index, double parameter){
		std::vector<unsigned int> indices_of_falloff_reactions = data_.kineticsMapXML()->IndicesOfFalloffReactions();
		int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(), 
			indices_of_falloff_reactions.end(), 
			index
		) - indices_of_falloff_reactions.begin();
		
		if(type == "lnA")
			data_.kineticsMapXML_->Set_A_falloff_inf(pos_FallOff_Reaction, std::exp(parameter));
		if(type == "Beta")
			data_.kineticsMapXML_->Set_Beta_falloff_inf(pos_FallOff_Reaction, parameter);
		if(type == "E_over_R")
			data_.kineticsMapXML_->Set_E_over_R_falloff_inf(pos_FallOff_Reaction, parameter);
	}

	void SimulationsInterface::ChangeThirdBodyEfficiencies(unsigned int i, std::string name, double parameter){
		// Finding position of the fall off reaction
        // Finding the index of the third body species
        int iSpecies = data_.thermodynamicsMapXML()->IndexOfSpecies(name);

        // Finding position of the third body species to be changed
        // Changing the value of the thirdbody species
        data_.kineticsMapXML_->Set_ThirdBody(i-1, iSpecies-1, parameter);
	}

} // namespace OptiSMOKE