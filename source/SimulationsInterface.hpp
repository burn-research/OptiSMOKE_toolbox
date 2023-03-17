
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

	void SimulationsInterface::run(){

		OpenSMOKE::KineticsMap_CHEMKIN* kinetics = data_->kineticsMapXML();
		OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermo = data_->thermodynamicsMapXML();

		// batch_reactors = new OptiSMOKE::BatchReactor[n_batch];
		for(unsigned int i = 0; i < n_batch; i++){
			batch_reactors[i] = new OptiSMOKE::BatchReactor[data_->input_paths()[i].size()];
			for(unsigned int j = 0; j < data_->input_paths()[i].size(); j++)
				batch_reactors[i][j].Setup(data_->input_paths()[i][j], thermo, kinetics);

		}

	}

} // namespace OptiSMOKE