/*-----------------------------------------------------------------------*\
|     ____            _  ______ __  __  ____  _  ________                 |
|    / __ \       _  (_)/  ___ |  \/  |/ __ \| |/ /  ____|                |
|   | |  | |_ __ | |_ _|  (___ | \  / | |  | | ' /| |__    _     _        |
|   | |  | | '_ \|  _| |\___  \| |\/| | |  | |  < |  __| _| |_ _| |_      |
|   | |__| | |_) | |_| |____)  | |  | | |__| | . \| |___|_   _|_   _|     |
|    \____/| .__/\___|_|______/|_|  |_|\____/|_|\_\______||_|   |_|       |
|          | |                                                            |
|          |_|                                                            |
|                                                                         |
|                                                                         |
|            Authors: Magnus Fürst <magnus.furst@ulb.ac.be>               |
|                     Andrea Bertolino <andrea.bertolino@ulb.be>          |
|                     Timoteo Dinelli <timoteo.dinelli@polimi.it>         |
|                                                                         |
|                                                                         |
|-------------------------------------------------------------------------|
|   License                                                               |
|                                                                         |
|   This file is part of OptiSMOKE.                                       |
|   Copyright (C) 2020 by Magnus Fürst and Andrea Bertolino               |
|                                                                         |
|   OptiSMOKE is free software: you can redistribute it and/or modify     |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OptiSMOKE is distributed in the hope that it will be useful,          |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OptiSMOKE. If not, see <http://www.gnu.org/licenses/>.     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "OptiSMOKEpp.h"

int main(int argc, char *argv[]) {

  #ifdef HAVE_AMPL
    // Switch to 53-bit rounding if appropriate, to eliminate some
    // cross-platform differences.
    fpinit_ASL();
  #endif

    // whether running in parallel
    bool parallel = Dakota::MPIManager::detect_parallel_launch(argc, argv);

    // Define MPI_DEBUG in dakota_global_defs.cpp to cause a hold here
    Dakota::mpi_debug_hold();

  #ifdef DAKOTA_HAVE_MPI
    if (parallel)
      MPI_Init(&argc, &argv); // initialize MPI
  #endif                      // DAKOTA_HAVE_MPI

    const char *dakota_input_string = NULL;

    OptiSMOKE::OptiSMOKE_logo("OptiSMOKE++",
                              "M. Furst, A. Bertolino, T. Dinelli");

    input.SetInputOptions(argc, argv);
    input.ReadDictionary();
    input.ReadExperimentalDataFiles();

    if (input.optimization_library() == "dakota") {
      input.DakotaInputString();
      dakota_input_string = input.dakota_input_string().c_str();
      run_dakota_parse(dakota_input_string,
                       input.dakota_options().echo_dakota_string());

      // Note: Dakota objects created in above function calls need to go
      // out of scope prior to MPI_Finalize so that MPI code in
      // destructors works properly in library mode.

  #ifdef DAKOTA_HAVE_MPI
      if (parallel)
        MPI_Finalize(); // finalize MPI
  #endif                // DAKOTA_HAVE_MPI

      return 0;
    }
    else if (input.optimization_library() == "nlopt") {

      violated_uncertainty = false;
      numberOfGradientEvaluations = 0;
      numberOfFunctionEvaluations = 0;

      sim_iface_ = new OptiSMOKE::SimulationsInterface(input);
      opti_kinetics_ = new OptiSMOKE::OptimizedKinetics(
          input, input.thermodynamicsMapXML_, input.kineticsMapXML_);

      sim_iface_->Setup();
      opti_kinetics_->SetChemkinName(input.output_folder() /
                                     input.optimized_kinetics_folder() /
                                     "OptimalMechanism.CKI");

      input.SetUpNLOPT();

      nlopt::algorithm algo = static_cast<nlopt::algorithm>(input.nlopt_options().algo_int());
      nlopt::opt opt(algo, input.optimization_target().number_of_parameters());

      opt.set_lower_bounds(input.lb());
      opt.set_upper_bounds(input.ub());
      opt.set_min_objective(NLOptFunction, NULL);
      opt.set_maxeval(input.nlopt_options().max_function_evaluations());
      opt.set_ftol_abs(1e-8);
      opt.set_ftol_rel(1e-6);

      std::vector<double> initial_values = input.initial_values();
      double minf;
      try {
        nlopt::result result = opt.optimize(initial_values, minf);
      } catch (std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
      }
    }
    else if (input.optimization_library() == "optimlib") {
      OptiSMOKE::FatalErrorMessage("OptimLIB not yet implemented!");
    } else {
      OptiSMOKE::FatalErrorMessage(
          "Available libraries for the optimization are: DAKOTA | NLOPT");
    }
  }

  void run_dakota_parse(const char *dakota_input_string,
                        bool echo_dakota_string) {

    // Parse input and construct Dakota LibraryEnvironment,
    // performing input data checks
    Dakota::ProgramOptions opts;
    opts.input_string(dakota_input_string);
    opts.echo_input(echo_dakota_string);

    // Defaults constructs the MPIManager, which assumes COMM_WORLD
    Dakota::LibraryEnvironment env(opts);

    if (env.mpi_manager().world_rank() == 0)
      Cout << "Library mode 1: run_dakota_parse()\n";

    // plug the client's interface (function evaluator) into the Dakota
    // environment; in serial case, demonstrate the simpler plugin method
    if (env.mpi_manager().mpirun_flag())
      // parallel_interface_plugin(env);
      OptiSMOKE::ErrorMessage("run_dakota_parse",
                              "Parallel interface not implemented yet!");
    else
      opensmoke_interface_plugin(env);

    // Execute the environment
    env.execute();
  }

  void opensmoke_interface_plugin(Dakota::LibraryEnvironment &env) {

    std::string model_type(""); // demo: empty string will match any model
    std::string interf_type("direct");
    std::string an_driver("opensmoke_plugin");

    Dakota::ProblemDescDB &problem_db = env.problem_description_db();

    std::shared_ptr<Dakota::Interface> serial_iface =
        std::make_shared<SIM::SerialDakotaInterface>(problem_db, input);

    bool plugged_in =
        env.plugin_interface(model_type, interf_type, an_driver,
        serial_iface);

    if (!plugged_in) {
      Cerr << "Error: no serial interface plugin performed. Check  "
           << "compatibility between parallel\n       configuration and "
           << "selected analysis_driver." << std::endl;
      Dakota::abort_handler(-1);
    }
  }

  double NLOptFunction(const std::vector<double> &x, std::vector<double> & grad, void *my_func_data)
{
    std::cout << "-------------------------------------------------------"
              << std::endl;
    std::cout << " Begin evaluation " << numberOfFunctionEvaluations + 1
              << std::endl;
    std::cout << "-------------------------------------------------------"
              << std::endl;

    for (unsigned int i = 0; i < x.size(); i++)
      std::cout << "\t" << std::scientific << std::setw(35) << std::left
                << input.param_str()[i] << std::scientific << std::setw(35)
                << std::left << std::setprecision(5) << x[i] << std::endl;

    const double f = OptFunction(x, numberOfFunctionEvaluations);
    if (!grad.empty()) {
      // FOR THE MOMENT NO GRADIENT METHODS
      numberOfGradientEvaluations++;
    }

    numberOfFunctionEvaluations++;
    return f;
  }

  double OptFunction(const std::vector<double> &b, unsigned int eval_nr) {
    double fn_val;
    sim_iface_->SubstituteKineticParameters(b);
    if (input.optimization_setup().penalty_function())
      violated_uncertainty = sim_iface_->CheckKineticConstasts();

    if (violated_uncertainty) {
      if (input.optimization_setup().objective_function_type() ==
      "CurveMatching")
        fn_val = 1;
      else
        fn_val = 10000000;
    } else {
      sim_iface_->run();
      fn_val = sim_iface_->ComputeObjectiveFunction();

      if (eval_nr == 0) {
        prev_fn_val = fn_val;
        sim_iface_->PrepareASCIIFile(fOut, input.parametric_file_name(),
                                     input.param_str());
      }
      if (prev_fn_val > fn_val) {
        prev_fn_val = fn_val;
        opti_kinetics_->WriteOptimizedMechanism();
        std::cout << " * Wrote optimized mechanism" << std::endl;
      }
    }
    sim_iface_->PrintASCIIFile(fOut, eval_nr, b, fn_val);
    return fn_val;
}
