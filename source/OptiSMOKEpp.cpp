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

/// Default Dakota input string for parallel case (text_book)
static const char parallel_input[] = 
  "	method,"
  "		optpp_q_newton"
  "		  max_iterations = 50"
  "		  convergence_tolerance = 1e-4"
  "	variables,"
  "		continuous_design = 2"
  "		  descriptors 'x1' 'x2'"
  "	interface,"
  "		direct"
  "		  analysis_driver = 'plugin_text_book'"
  "	responses,"
  "		num_objective_functions = 1"
  "		num_nonlinear_inequality_constraints = 2"
  "		analytic_gradients"
  "		no_hessians";
  
int main(int argc, char* argv[]){
    
    // whether running in parallel
    //bool parallel = Dakota::MPIManager::detect_parallel_launch(argc, argv);
    
	#ifdef OPTISMOKE_USE_MPI
    MPI_Init(&argc, &argv);
    #endif

    int rank, nprocs;
    #ifdef OPTISMOKE_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    #else
    rank = 0;
    nprocs = 0;
    #endif
    
    if(rank == 0)
        OptiSMOKE::OptiSMOKE_logo("OptiSMOKE++", "M. Furst, A. Bertolino, T. Dinelli");
    
    #ifdef OPTISMOKE_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    boost::posix_time::ptime t_program_start = boost::posix_time::second_clock::local_time();
    
    input.SetInputOptions(argc, argv);
    input.ReadDictionary();
    input.ReadExperimentalDataFiles();
    input.loadDistributions();
    
    // Dunno if this barrier is needed
    #ifdef OPTISMOKE_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    if(input.optimization_library() == "dakota"){
    
        const char *dakota_input_string = NULL;
        input.DakotaInputString();
        dakota_input_string = input.dakota_input_string().c_str();

        run_dakota_parse(dakota_input_string, input.dakota_options().echo_dakota_string(), rank);
    }
    #ifdef OPTISMOKE_USE_NLOPT
    else if (input.optimization_library() == "nlopt"){
        
        violated_uncertainty = false;
        numberOfGradientEvaluations = 0;
        numberOfFunctionEvaluations = 0;
        
        sim_iface_ = new OptiSMOKE::SimulationsInterface(input);
        opti_kinetics_ = new OptiSMOKE::OptimizedKinetics(input, input.thermodynamicsMapXML_, input.kineticsMapXML_);

        sim_iface_->Setup();
        opti_kinetics_->SetChemkinName(input.output_folder() / input.optimized_kinetics_folder() / "OptimalMechanism.CKI");

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
        
        #ifdef OPTISMOKE_USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        #endif
        
        double minf;
        try{
            nlopt::result result = opt.optimize(initial_values, minf);
        }
        catch(std::exception &e) {
            if(input.isMaster())
                std::cout << "nlopt failed: " << e.what() << std::endl;
        }
    }
    #endif // OPTISMOKE_USE_NLOPT
    else if (input.optimization_library() == "optimlib"){
        if(input.isMaster())
            OptiSMOKE::FatalErrorMessage("OptimLIB not yet implemented!");
    }
    else{
        if(input.isMaster())
            OptiSMOKE::FatalErrorMessage("Available libraries for the optimization are: DAKOTA | NLOPT");
    }
    
    #ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    boost::posix_time::ptime t_program_end = boost::posix_time::second_clock::local_time();
    boost::posix_time::time_duration tdiff = t_program_end - t_program_start;

    if (rank == 0) {
        std::cout << "Time elapsed: ";
        std::cout << std::fixed << tdiff.total_seconds() << " seconds." << std::endl;
    }
    
    #ifdef OPTISMOKE_USE_MPI
    MPI_Finalize(); // finalize MPI
    #endif // OPTISMOKE_USE_MPI
    
    return 0;
}

void run_dakota_parse(const char* dakota_input_string, bool echo_dakota_string, int rank){

    // Parse input and construct Dakota LibraryEnvironment, 
    // performing input data checks
    Dakota::ProgramOptions opts(0);
    opts.input_string(parallel_input);
    opts.echo_input(echo_dakota_string);
    // opts.world_rank();
    // opts.write_restart_file();

    // Defaults constructs the MPIManager, which assumes COMM_WORLD
    Dakota::LibraryEnvironment env(opts);

    if (rank == 0)
        std::cout << "\nDakota Library mode 1: run_dakota_parse()\n";
    
    // plug the client's interface (function evaluator) into the Dakota
    // environment; in serial case, demonstrate the simpler plugin method
    // #if OPTISMOKE_USE_MPI
    // if(rank == 0)
    //     OptiSMOKE::FatalErrorMessage("Parallell interface not available yet!");
    if(env.mpi_manager().mpirun_flag())
		if(rank == 0)
			OptiSMOKE::FatalErrorMessage("Parallell interface for dakota not available yet!");
		//parallel_interface_plugin(env);
    else
        opensmoke_interface_plugin(env);    


    // Execute the environment
    env.execute();
}

void opensmoke_interface_plugin(Dakota::LibraryEnvironment& env){    
    
    std::string model_type(""); // demo: empty string will match any model type
    std::string interf_type("direct");
    std::string an_driver("opensmoke_plugin");

    Dakota::ProblemDescDB& problem_db = env.problem_description_db();

    std::shared_ptr<Dakota::Interface> serial_iface = 
        std::make_shared<SIM::SerialDakotaInterface>(problem_db, input);
    
    bool plugged_in = env.plugin_interface(model_type, interf_type, an_driver, serial_iface);
            
    if (!plugged_in) {
        Cerr << "Error: no serial interface plugin performed. Check  "
            << "compatibility between parallel\n       configuration and "
            << "selected analysis_driver." << std::endl;
        Dakota::abort_handler(-1);
    }
}

#ifdef OPTISMOKE_USE_MPI
/*void parallel_interface_plugin(Dakota::LibraryEnvironment& env)
{
    // std::string model_type(""); // demo: empty string will match any model type
    // std::string interf_type("direct");
    // std::string an_driver("opensmoke_plugin");
    // Dakota::ProblemDescDB& problem_db = env.problem_description_db();
    // std::cout << "Diocane" << std::endl;
     
    // get the list of all models matching the specified model, interface, driver:
    Dakota::ModelList filt_models = env.filtered_model_list("simulation", "direct", "plugin_text_book");

    if (filt_models.empty()) {
        Cerr << "Error: no parallel interface plugin performed.  Check "
	        << "compatibility between parallel\n       configuration and "
	        << "selected analysis_driver." << std::endl;
        Dakota::abort_handler(-1);
    }

    Dakota::ProblemDescDB& problem_db = env.problem_description_db();
    Dakota::ModelLIter ml_iter;
    size_t model_index = problem_db.get_db_model_node(); // for restoration
    for (ml_iter = filt_models.begin(); ml_iter != filt_models.end(); ++ml_iter) {
        // set DB nodes to input specification for this Model
        problem_db.set_db_model_nodes(ml_iter->model_id());

        Dakota::Interface& model_interface = ml_iter->derived_interface();

        // Parallel case: plug in derived Interface object with an analysisComm.
        // Note: retrieval and passing of analysisComm is necessary only if
        // parallel operations will be performed in the derived constructor.

        // retrieve the currently active analysisComm from the Model.  In the most
        // general case, need an array of Comms to cover all Model configurations.
        const MPI_Comm& analysis_comm = ml_iter->analysis_comm();

        // don't increment ref count since no other envelope shares this letter
        model_interface.assign_rep(std::make_shared<SIM::ParallelDirectApplicInterface>(problem_db, analysis_comm));
    }
    problem_db.set_db_model_nodes(model_index);
}*/
#endif // OPTISMOKE_USE_MPI

#ifdef OPTISMOKE_USE_NLOPT
double NLOptFunction(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if(input.isMaster()){
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << " Begin evaluation " << numberOfFunctionEvaluations + 1 << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;

        for(unsigned int i = 0; i < x.size(); i++)
            std::cout << "\t" << std::scientific << std::setw(35) << std::left
                << input.param_str()[i] << std::scientific << std::setw(35) 
                << std::left << std::setprecision(5) << x[i] << std::endl;
    }

    #ifdef OPTISMOKE_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    const double f = OptFunction(x, numberOfFunctionEvaluations);
	
    #ifdef OPTISMOKE_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    if (!grad.empty())
	{
        // FOR THE MOMENT NO GRADIENT METHODS
		numberOfGradientEvaluations++;
	}
	numberOfFunctionEvaluations++;
	return f;
}

double OptFunction(const std::vector<double>& b, unsigned int eval_nr)
{
    double fn_val;
    sim_iface_->SubstituteKineticParameters(b);
	if(input.optimization_setup().penalty_function())
		violated_uncertainty = sim_iface_->CheckKineticConstasts();
    
    #ifdef OPTISMOKE_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
	
    if(violated_uncertainty){
		if (input.optimization_setup().objective_function_type() == "CurveMatching")
			fn_val = 1;
		else
			fn_val = 10000000;
	}
    else{
        
        sim_iface_->run();
        fn_val = sim_iface_->ComputeObjectiveFunction();
        
        #ifdef OPTISMOKE_USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        #endif
        
        if(eval_nr == 0){
            prev_fn_val = fn_val;
            if(input.isMaster())
                sim_iface_->PrepareASCIIFile(fOut, input.parametric_file_name(), input.param_str());
        }
        if(prev_fn_val > fn_val) {
            prev_fn_val = fn_val;
            if(input.isMaster()){
                opti_kinetics_->WriteOptimizedMechanism();
                std::cout << " * Wrote optimized mechanism" << std::endl;
            }
        }
    }    
    if(input.isMaster())
        sim_iface_->PrintASCIIFile(fOut, eval_nr, b, fn_val);
    
    #ifdef OPTISMOKE_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    return fn_val;
}
#endif // OPTISMOKE_USE_NLOPT
