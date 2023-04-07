#ifndef OPTIMIZEDKINETICS_H
#define OPTIMIZEDKINETICS_H

namespace OptiSMOKE{
	
	class OptimizedKinetics
	{
	public:
		OptimizedKinetics(const OptiSMOKE::InputManager& data,
						const OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
						OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML);
	
		~OptimizedKinetics ();
	
		void SetChemkinName (const fs::path& path);
	
		void WriteOptimizedMechanism ();

	private:
		
		//References
		const OptiSMOKE::InputManager& data_;
		const OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;

		int NS_, NR_;
		fs::path chemkin_path_;
		bool isChemkinNameSet_;

		std::ofstream fChemKin_;
		std::ofstream fSpecies_;

		std::vector<std::string> species_list_;
		std::vector<bool> iSpeciesList_;

		PreProcessorKinetics_CHEMKIN* preprocessor_kinetics_;

		void MemoryAllocation ();
		void CheckElementPresence ();
		void SelectReactions();
	};
} // namespace OptiSMOKE

#include "OptimizedKinetics.hpp"

#endif // OPTIMIZEDKINETICS_H