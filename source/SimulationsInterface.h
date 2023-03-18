#ifndef SIMULATIONS_INTERFACE_H
#define SIMULATIONS_INTERFACE_H

namespace OptiSMOKE {

	class SimulationsInterface
	{
  	public:
		SimulationsInterface(OptiSMOKE::InputManager* data);

		~SimulationsInterface();
  	
		void run();

		void Setup();
		
  	private:
		OptiSMOKE::InputManager* data_;
	
		std::vector<OptiSMOKE::BatchReactor*> batch_reactors;
	
		unsigned int n_batch;
		unsigned int n_pfr;
		unsigned int n_psr;
		unsigned int n_premixed;
		unsigned int n_counterflow;
	};
} // namespace OptiSMOKE

#include "SimulationsInterface.hpp"
#endif // SIMULATIONS_INTERFACE_H