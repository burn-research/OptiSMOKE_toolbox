#ifndef SERIAL_DAKOTA_INTERFACE_H
#define SERIAL_DAKOTA_INTERFACE_H

#include "SimulationsInterface.h"

namespace SIM {

class SerialDakotaInterface: public Dakota::DirectApplicInterface
{
  public:

    // constructor
    SerialDakotaInterface(const Dakota::ProblemDescDB& problem_db, const OptiSMOKE::InputManager& data);
  
    // destructor
    ~SerialDakotaInterface();

  protected:

    // execute an analysis code portion of a direct evaluation invocation
    int derived_map_ac(const Dakota::String& ac_name);

    void derived_map_asynch(const Dakota::ParamResponsePair& pair);

    // evaluate the batch of jobs contained in prp_queue
    void wait_local_evaluations(Dakota::PRPQueue& prp_queue);
  
    // invokes wait_local_evaluations() (no special nowait support)
    void test_local_evaluations(Dakota::PRPQueue& prp_queue);

    // no-op hides default run-time error checks at DirectApplicInterface level
    void set_communicators_checks(int max_eval_concurrency);

  private:

    // Rosenbrock plug-in test function
    // Keep this for future ideas and testing maybe?
    // int rosenbrock(const Dakota::RealVector& c_vars, short asv,
    //              Dakota::Real& fn_val, Dakota::RealVector& fn_grad,
    //             Dakota::RealSymMatrix& fn_hess);
    //

    int simulations_interface(const Dakota::RealVector& c_vars, short asv,Dakota::Real& fn_val);
         
    const OptiSMOKE::InputManager& data_;

    OptiSMOKE::SimulationsInterface* sim_iface_;

    OptiSMOKE::OptimizedKinetics* opti_kinetics_;

    int eval_nr;

    double prev_fn_val;

	  bool violated_uncertainty;
};

  // Constructor
  inline SerialDakotaInterface::SerialDakotaInterface (
    const Dakota::ProblemDescDB& problem_db, const OptiSMOKE::InputManager& data
  ) : Dakota::DirectApplicInterface(problem_db), data_(data)
  {
    eval_nr = 0;
    violated_uncertainty = false;
    
    sim_iface_ = new OptiSMOKE::SimulationsInterface(data_);
    sim_iface_->Setup();

    opti_kinetics_ = new OptiSMOKE::OptimizedKinetics(data_, 
      data_.thermodynamicsMapXML_, 
      data_.kineticsMapXML_);
  }

  // Destructor
  inline SerialDakotaInterface::~SerialDakotaInterface(){ }


  inline void SerialDakotaInterface::derived_map_asynch(const Dakota::ParamResponsePair& pair)
  {
    // no-op (just hides base class error throw). Jobs are run exclusively within
    // wait_local_evaluations(), prior to there existing true batch processing
    // facilities.
  }


  /** For use by ApplicationInterface::serve_evaluations_asynch(), which can
    provide a batch processing capability within message passing schedulers
    (called using chain IteratorScheduler::run_iterator() --> Model::serve()
    --> ApplicationInterface::serve_evaluations()
    --> ApplicationInterface::serve_evaluations_asynch()). */
  inline void SerialDakotaInterface::test_local_evaluations(Dakota::PRPQueue& prp_queue)
  { 
    wait_local_evaluations(prp_queue);
  }

  // Hide default run-time error checks at DirectApplicInterface level
  inline void SerialDakotaInterface::set_communicators_checks(int max_eval_concurrency){ }

} // namespace SIM

#include "SerialDakotaInterface.hpp"

#endif // SERIAL_DAKOTA_INTERFACE_H
