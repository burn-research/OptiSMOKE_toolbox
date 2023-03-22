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

		double ComputeObjectiveFunction();
	
	private:
		OptiSMOKE::InputManager* data_;
	
		std::vector<OptiSMOKE::BatchReactor*> batch_reactors;
	
		unsigned int n_batch;
		unsigned int n_pfr;
		unsigned int n_psr;
		unsigned int n_premixed;
		unsigned int n_counterflow;

		std::vector<std::vector<double>> k_upper;
		std::vector<std::vector<double>> k_lower;

		std::vector<std::vector<double>> k_upper_inf;
		std::vector<std::vector<double>> k_lower_inf;

		std::vector<std::vector<std::vector<double>>> k_upper_classic_plog;
		std::vector<std::vector<std::vector<double>>> k_lower_classic_plog;

		std::vector<std::vector<std::vector<double>>> simulations_results_;
	
		std::vector<std::vector<std::vector<std::vector<double>>>> bootstrap_exp;
		
		void BootstrappingData(std::vector<std::vector<std::vector<double>>> ciao);

	};
} // namespace OptiSMOKE

#include "SimulationsInterface.hpp"
#endif // SIMULATIONS_INTERFACE_H