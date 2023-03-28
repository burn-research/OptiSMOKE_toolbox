#ifndef SIMULATIONS_INTERFACE_H
#define SIMULATIONS_INTERFACE_H

namespace OptiSMOKE {

	class SimulationsInterface
	{
  	public:
		SimulationsInterface(const OptiSMOKE::InputManager& data);

		~SimulationsInterface();
  	
		void run();

		void Setup();

		double ComputeObjectiveFunction();

		bool CheckKineticConstasts();
			
		void SubstituteKineticParameters(const Dakota::RealVector& c_vars);
	
	private:
		const OptiSMOKE::InputManager& data_;
	
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

		std::vector<double> T_span = {300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, \
		1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500};

		void ChangeDirectParamaters(std::string type, int index, double parameter);

		void ChangeFallOffParamaters(std::string type, int index, double parameter);

		void ChangeThirdBodyEfficiencies(unsigned int i, std::string name, double parameter);
	
	};
} // namespace OptiSMOKE

#include "SimulationsInterface.hpp"
#endif // SIMULATIONS_INTERFACE_H