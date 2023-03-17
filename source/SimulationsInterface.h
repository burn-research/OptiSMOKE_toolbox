#ifndef SIMULATIONS_INTERFACE_H
#define SIMULATIONS_INTERFACE_H

namespace OptiSMOKE {

class SimulationsInterface
{
  public:
	SimulationsInterface(const OptiSMOKE::InputManager* data);
	~SimulationsInterface();
  protected:
  private:
	const OptiSMOKE::InputManager* data_;
};
} // namespace OptiSMOKE

#include "SimulationsInterface.hpp"
#endif // SIMULATIONS_INTERFACE_H