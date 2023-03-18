
namespace OptiSMOKE{

	SimulationsInterface::SimulationsInterface(OptiSMOKE::InputManager* data)
	{
		data_ = data;
		n_batch = data_->optimization_target().number_of_batch_reactor();
		n_pfr = data_->optimization_target().number_of_plug_flow_reactor();
		n_psr = data_->optimization_target().number_of_perfectly_stirred_reactor();
		n_premixed = data_->optimization_target().number_of_premixed_laminar_flame();
		n_counterflow = data_->optimization_target().number_of_counter_flow_flame();

		batch_reactors.resize(n_batch);
	}

	SimulationsInterface::~SimulationsInterface(){}

	void SimulationsInterface::Setup(){

		if(n_batch != 0){
			batch_reactors.resize(n_batch);
			for(unsigned int i = 0; i < n_batch; i++)
				batch_reactors[i] = new OptiSMOKE::BatchReactor[data_->input_paths()[i].size()];
		}
	}

	void SimulationsInterface::run(){

		OpenSMOKE::KineticsMap_CHEMKIN* kinetics = data_->kineticsMapXML();
		OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermo = data_->thermodynamicsMapXML();

		for(unsigned int i = 0; i < n_batch; i++){
			std::string qoi = data_->QoI()[i];
			std::string qoi_target = data_->QoI_target()[i];
			if(qoi == "IDT"){
				for(unsigned int j = 0; j < data_->input_paths()[i].size(); j++){
					// For the moment stay simple but keep in mind that this every time 
					// re-read the OS input file
					batch_reactors[i][j].Setup(data_->input_paths()[i][j], thermo, kinetics);
					batch_reactors[i][j].Solve();
					std::cout << batch_reactors[i][j].GetIgnitionDelayTime(qoi_target) << std::endl;
				}
			}
		}
	}

} // namespace OptiSMOKE