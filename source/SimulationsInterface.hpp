
namespace OptiSMOKE{

	SimulationsInterface::SimulationsInterface(OptiSMOKE::InputManager* data)
	{
		// Set options
		data_ = data;

		// Resize reactors objects
		
		n_batch = data_->optimization_target().number_of_batch_reactor();
		n_pfr = data_->optimization_target().number_of_plug_flow_reactor();
		n_psr = data_->optimization_target().number_of_perfectly_stirred_reactor();
		n_premixed = data_->optimization_target().number_of_premixed_laminar_flame();
		n_counterflow = data_->optimization_target().number_of_counter_flow_flame();

		batch_reactors.resize(n_batch);

		// Resize simulations results
		simulations_results_.resize(data_->expdata_x().size());
		for(unsigned int i = 0; i < data_->expdata_x().size(); i++){
			simulations_results_[i].resize(data_->expdata_x()[i].size());
		}
	}

	SimulationsInterface::~SimulationsInterface(){}

	void SimulationsInterface::Setup(){

		// Allocate Reactors Object
		if(n_batch != 0){
			batch_reactors.resize(n_batch);
			for(unsigned int i = 0; i < n_batch; i++)
				batch_reactors[i] = new OptiSMOKE::BatchReactor[data_->input_paths()[i].size()];
		}


		//-------------------------------------------------------//
		//               Reactions constraints                   //
		//-------------------------------------------------------//

		// Create constraints for reactions
		// Initializes vector of uncertainty factors for Direct reactions
		// Compute max and min reaction rates according to temperature dependent 
		// or constant uncertainty factors for direct reactions
		std::vector<double> T_span = {300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500};
		std::vector<double> f_factors = data_->optimization_target().list_of_uncertainty_factors();

		// Computing boundaries for all the reactions
		int rows = data_->optimization_target().list_of_target_uncertainty_factors().size();
		int columns = T_span.size();
		std::vector<std::vector<double>> k_0(rows, std::vector<double>(columns));
		k_upper.resize(rows, std::vector<double>(columns));
		k_lower.resize(rows, std::vector<double>(columns));

		// if (ObjectInput2.udc_bool == true){
		// Not yet added
		// } 
		// else{
		for (int j = 0; j < rows; j++){
			int reaction_index = data_->optimization_target().list_of_target_uncertainty_factors()[j];
			double A_0 = data_->nominalkineticsMapXML()->A(reaction_index-1);
			double Beta_0 = data_->nominalkineticsMapXML()->Beta(reaction_index-1);
			double E_over_R_0 = data_->nominalkineticsMapXML()->E_over_R(reaction_index-1);
			for (int i=0; i < columns; i++){
				k_0[j][i] = A_0 * std::pow(T_span[i], Beta_0) * std::exp((-1*E_over_R_0)/T_span[i]);
				k_upper[j][i] = k_0[j][i] * std::pow(10, f_factors[j]) * ((double)data_->optimization_setup().sigma_k_distribution()/2);
				k_lower[j][i] = k_0[j][i] * std::pow(10,-1*f_factors[j]) / ((double)data_->optimization_setup().sigma_k_distribution()/2);
			}
		}
		// }

		// Compute max and min reaction rates according to temperature dependent 
		// or constant uncertainty factors for fall off reactions
		std::vector<double> f_factors_inf = data_->optimization_target().list_of_uncertainty_factors_inf();
		int rows_inf = data_->optimization_target().list_of_target_uncertainty_factors_inf().size();
		int columns_inf = T_span.size();
		std::vector<unsigned int> indices_of_falloff_reactions = data_->nominalkineticsMapXML()->IndicesOfFalloffReactions();
		std::vector<std::vector<double>> k_0_inf(rows_inf, std::vector<double>(columns_inf));
		k_upper_inf.resize(rows_inf, std::vector<double>(columns_inf));
		k_lower_inf.resize(rows_inf, std::vector<double>(columns_inf));

		for (int j=0; j < rows_inf; j++){
			int falloff_reaction_index = data_->optimization_target().list_of_target_uncertainty_factors_inf()[j];
			int pos_FallOff_Reaction = std::find(indices_of_falloff_reactions.begin(), indices_of_falloff_reactions.end(), falloff_reaction_index) - indices_of_falloff_reactions.begin();
			double A_0_inf = data_->nominalkineticsMapXML()->A_falloff_inf(pos_FallOff_Reaction);
			double Beta_0_inf = data_->nominalkineticsMapXML()->Beta_falloff_inf(pos_FallOff_Reaction);
			double E_over_R_0_inf = data_->nominalkineticsMapXML()->E_over_R_falloff_inf(pos_FallOff_Reaction);
			for (int i=0; i < columns_inf; i++){
				k_0_inf[j][i] = A_0_inf * std::pow(T_span[i],Beta_0_inf) * std::exp((-1*E_over_R_0_inf)/T_span[i]);
				k_upper_inf[j][i] = k_0_inf[j][i] * std::pow(10,f_factors_inf[j])* (double)data_->optimization_setup().sigma_k_distribution()/2;
				k_lower_inf[j][i] = k_0_inf[j][i] * std::pow(10,-1*f_factors_inf[j])/ ((double)data_->optimization_setup().sigma_k_distribution()/2);
			}
		}

		// Compute max and min reaction rates according to temperature dependent 
		// or constant uncertainty factors for for PLOG reactions
		std::vector<std::vector<std::vector<double>>> k_classic_plog;
		std::vector<unsigned int> indices_of_classic_plog = data_->nominalkineticsMapXML()->IndicesOfPLOGReactions();
		int rows_cp = data_->optimization_target().list_of_uncertainty_factors_classic_plog().size();
		//resize first dimensions according to the number of classic plog we are interested in for optimisation
    	k_classic_plog.resize(rows_cp);		
		k_upper_classic_plog.resize(rows_cp);
		k_lower_classic_plog.resize(rows_cp);
	
		for (unsigned int j=0; j < rows_cp; j++){

			//finding the position of the target reaction in correspondent OS++ object
			int pos_classic_plog_reaction = std::find(
				indices_of_classic_plog.begin(), 
				indices_of_classic_plog.end(), 
				data_->optimization_target().list_of_target_classic_plog_reactions()[j]
			)-indices_of_classic_plog.begin();
			
			for(int k=0; k<data_->nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA().size(); k++){
				if(data_->nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA()[k].size() != 1){
					std::cout << " * WARNING one of the plog reactions requested for the optimization is a duplicate!" << std::endl;
				}
			}	
			
			k_classic_plog[j].resize(data_->nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA().size());
			k_upper_classic_plog[j].resize(data_->nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA().size());
			k_lower_classic_plog[j].resize(data_->nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA().size());

			for (int k=0; k < data_->nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA().size(); k++){
				double A_CP = std::exp(data_->nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).lnA()[k][0]);
				double n_CP = data_->nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).Beta()[k][0];
				double E_over_R_CP = data_->nominalkineticsMapXML()->pressurelog_reactions(pos_classic_plog_reaction).E_over_R()[k][0];

				k_classic_plog[j][k].resize(T_span.size());
				k_upper_classic_plog[j][k].resize(T_span.size());
				k_lower_classic_plog[j][k].resize(T_span.size());

				for (int i=0; i < T_span.size(); i++){
					k_classic_plog[j][k][i] = A_CP * std::pow(T_span[i],n_CP) * std::exp((-1*E_over_R_CP)/T_span[i]);
					k_upper_classic_plog[j][k][i] = k_classic_plog[j][k][i] * std::pow(10, data_->optimization_target().list_of_uncertainty_factors_classic_plog()[j]) * ((double)data_->optimization_setup().sigma_k_distribution()/2);
					k_lower_classic_plog[j][k][i] = k_classic_plog[j][k][i] * std::pow(10, -1*data_->optimization_target().list_of_uncertainty_factors_classic_plog()[j])/ ((double)data_->optimization_setup().sigma_k_distribution()/2);
				}
			}		
		}

		//-------------------------------------------------------//
		//             Spline for experimental data              //
		//-------------------------------------------------------//

		if(data_->optimization_setup().objective_function_type() == "CurveMatching"){

			splines_exp.resize(data_->expdata_x().size());
			if(data_->curvematching_options().use_bootstrap() == false){
				for (int a=0; a < data_->expdata_x().size(); a++){	
					splines_exp[a].resize(data_->expdata_x()[a].size());
					for (int b = 0; b < data_->expdata_x()[a].size(); b++){
						splines_exp[a][b].resize(1);
						if(data_->QoI()[a] == "IDT"){
							// Initialize 1 vector to get the ordinates
							std::vector<double> temporary_vector;			
							temporary_vector.resize(data_->expdata_y()[a][b].size());
							for (int z=0; z < data_->expdata_y()[a][b].size(); z++){
								temporary_vector[z] = std::log(data_->expdata_y()[a][b][z]);
							}
							splines_exp[a][b][0].solve(data_->expdata_x()[a][b], temporary_vector, 0, 0, false);
						} 
						else {
							splines_exp[a][b][0].solve(data_->expdata_x()[a][b], data_->expdata_y()[a][b], 0, 0, false);	
						}
						splines_exp[a][b][0].removeAsymptotes();	
					}
				}
			}
			else{
				// Bootstrap is active

				
			}
		}
		else if(data_->optimization_setup().objective_function_type() == "L2-norm"){

		}
		else if(data_->optimization_setup().objective_function_type() == "L1-norm"){
		}
		else{
			OptiSMOKE::FatalErrorMessage("Unknown type for the Objective Function, available are: CurveMatching | L2-norm | L1-norm");
		}
	}

	void SimulationsInterface::run(){

		OpenSMOKE::KineticsMap_CHEMKIN* kinetics = data_->kineticsMapXML();
		OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermo = data_->thermodynamicsMapXML();

		// Loop over all datasets
		// Here data_->path_experimental_data_files().size() this is 
		// misleading however keep in mind that only the size matters
		for(unsigned int i = 0; i < data_->path_experimental_data_files().size(); i++){
			std::string qoi = data_->QoI()[i];
			std::string qoi_target = data_->QoI_target()[i];
			std::string solver = data_->solver_name()[i];
			
			if (solver == "BatchReactor"){
				if(qoi == "IDT"){
					for(unsigned int j = 0; j < data_->input_paths()[i].size(); j++){
						// For the moment stay simple but keep in mind that this every time 
						// re-read the OS input file
						batch_reactors[i][j].Setup(data_->input_paths()[i][j], thermo, kinetics);
						batch_reactors[i][j].Solve();
						simulations_results_[i][0].push_back(batch_reactors[i][j].GetIgnitionDelayTime(qoi_target)*1e6);
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
				OptiSMOKE::FatalErrorMessage(solver + " not yet supported!");
			}

			if(solver == "PerfectlyStirredReactor"){
				OptiSMOKE::FatalErrorMessage(solver + " not yet supported!");
			}
			
			if (solver == "PremixedLaminarFlame1D"){
				OptiSMOKE::FatalErrorMessage(solver + " not yet supported!");
			}
            
			if (solver == "CounterFlowFlame1D"){
				OptiSMOKE::FatalErrorMessage(solver + " not yet supported!");
			}
		}
	}

	void SimulationsInterface::BootstrappingData(std::vector<std::vector<std::vector<double>>>){
	
		bootstrap_exp.resize(data_->expdata_y().size());

    	for (int a = 0; a < data_->expdata_y().size(); a++) {
			bootstrap_exp[a].resize(data_->expdata_y()[a].size()); 
			for (int b=0; b < data_->expdata_y()[a].size(); b++) {
				bootstrap_exp[a][b].resize(data_->curvematching_options().number_of_bootstrap()); 
				for (int c=0; c < data_->curvematching_options().number_of_bootstrap(); c++) {
					bootstrap_exp[a][b][c].resize(data_->expdata_y()[a][b].size()); 
					for (int d=0; d < data_->expdata_y()[a][b].size(); d++) {
						bootstrap_exp[a][b][c][d] = data_->expdata_y()[a][b][d];
					}
				}
			}		
		}
		
		// Replace BS variations with samples from a 
		// gaussian distribution

		for (int i=0; i < bootstrap_exp.size(); i++) {
			for (int c=0; c < bootstrap_exp[i].size(); c++){
				for (int a=0; a < bootstrap_exp[i][c][0].size(); a++) {
					// create an object of class std::default_random_engine, 
					// which generates pseudo-random numbers
					std::default_random_engine generator;
					// initializes the seed of the random_engine_generator
					generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
					std::normal_distribution<double> distribution(data_->expdata_y()[i][c][a], ObjectInput2.standard_deviations[i][c][a]);	
					// loop over the number of Bootstrap variations
					for (int b=1; b < data_->curvematching_options().number_of_bootstrap(); b++){
						// generate a random number from the distribution
						double number = distribution(generator);
						// if negative replace with 0
						if (number < 0 && data_->curvematching_options().possible_negative_ordinates() == false){
							number = 0;
						}
						// if logScale is false, returns the number
						// otherwise, return the log of the number to bootstrap.
						if (logScale == false){
							// if the data is 0, then always put 0, otherwise, replace the nominal value with the sampled number
							if (Exp_data[i][c][1][a]==0){
								bootstrapExp[i][c][b][a] = 0;
							} else{
								bootstrapExp[i][c][b][a] = number;
							}
						}
						else{
							// IN PRINCIPLES IT WILL NEVER GO HERE
							bootstrapExp[i][c][b][a] = log(number);
						}	
					}
				}
			}
		}
	}

} // namespace OptiSMOKE