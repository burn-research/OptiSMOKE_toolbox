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

int main(int argc, char* argv[]){

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
    #endif // DAKOTA_HAVE_MPI
    
    const char *dakota_input_string = NULL;

    OptiSMOKE::OptiSMOKE_logo("OptiSMOKE++", "M. Furst, A. Bertolino, T. Dinelli");
    
    input.SetInputOptions(argc, argv);
    input.ReadDictionary();
    input.ReadExperimentalDataFiles();

    if(input.optimization_library() == "dakota"){
        input.DakotaInputString();
        dakota_input_string = input.dakota_input_string().c_str();
        run_dakota_parse(dakota_input_string, input.dakota_options().echo_dakota_string());

        // Note: Dakota objects created in above function calls need to go
        // out of scope prior to MPI_Finalize so that MPI code in
        // destructors works properly in library mode.

        #ifdef DAKOTA_HAVE_MPI
        if (parallel)
            MPI_Finalize(); // finalize MPI
        #endif // DAKOTA_HAVE_MPI
    
        return 0;
    }
    else if (input.optimization_library() == "nlopt"){
        
        violated_uncertainty = false;
        numberOfGradientEvaluations = 0;
        numberOfFunctionEvaluations = 0;

        sim_iface_.Setup();
        opti_kinetics_.SetChemkinName(input.optimized_kinetics_folder() / "OptimalMechanism.CKI");

        nlopt_opt opt;
        opt = nlopt_create(NLOPT_LD_TNEWTON, input.optimization_target().number_of_parameters());

        /*double* lb = new double[number_parameters];// lower bounds
		double* ub = new double[number_parameters];// upper bounds
		double* x = new double[number_parameters];// first guess
		for (unsigned int i = 0; i < number_parameters; i++){
			lb[i] = bMin(i);
			ub[i] = bMax(i);
			x[i] = b0(i);
		}

		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
		nlopt_set_min_objective(opt, NLOptFunction, NULL);
		nlopt_set_maxeval(opt, max_eval);

		double fOpt;
		if (nlopt_optimize(opt, x, &fOpt) < 0)
		{
			std::cout << "NLopt failed!" << std::endl;
			getchar();
            exit(-1);
		}
		else
		{
			for (unsigned int i = 0; i < number_parameters; i++)
				bOpt(i) = x[i];
		}*/

    }
    else if (input.optimization_library() == "optimlib"){
        OptiSMOKE::FatalErrorMessage("OptimLIB not yet implemented!");
    }
    else{
        OptiSMOKE::FatalErrorMessage("Available libraries for the optimization are: DAKOTA | NLOPT");
    }
}

void run_dakota_parse(const char* dakota_input_string, bool echo_dakota_string){

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
        OptiSMOKE::ErrorMessage("run_dakota_parse", "Parallel interface not implemented yet!");
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


//  #if OPTISMOKE_USE_NLOPT
double NLOptFunction(unsigned n, const double *x, double *grad, void *my_func_data)
{
    int number_parameters = input.optimization_target().number_of_parameters();
	Eigen::VectorXd b(number_parameters);
	for (unsigned int i = 0; i < number_parameters; i++)
		b(i) = x[i];

	const double f = OptFunction(b);

	if (grad)
	{
		std::cout << "    * Gradient evaluation..." << std::endl;

		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double ETA3 = std::pow(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE, 1./3.);
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_DOUBLE);
		
		// Dimensions of parameters vector
		Eigen::VectorXd b_dimensions(input.optimization_target().number_of_parameters());
		b_dimensions.setConstant(1.);

		// Choosing between forward and centrale difference
		double eta = ETA2;
		if (central_difference == true)
			eta = ETA3;
		
		// Estimate the increment for forward approximation
		Eigen::VectorXd deltab(number_parameters);
		for (unsigned int i = 0; i < number_parameters; i++)
		{
			if (b(i) < 0.)
				deltab(i) = -eta * std::max(std::fabs(b(i)), std::fabs(b_dimensions(i)));
			else
				deltab(i) =  eta * std::max(std::fabs(b(i)), std::fabs(b_dimensions(i)));
			
			if (deltab(i) == 0.)
				deltab(i) = ZERO_DER;
		}

		// Forward gradient
		if (central_difference == false)
		{
			Eigen::VectorXd b_plus = b;
			for (unsigned int j = 0; j < number_parameters; j++)
			{
				b_plus(j) = b(j) + deltab(j);
				const double f_plus = OptFunction(b_plus);
				
				grad[j] = (f_plus - f) / deltab(j);

				b_plus(j) = b(j);
			}
		}

		// Central gradient
		if (central_difference == true)
		{
			Eigen::VectorXd b_star = b;
			for (unsigned int j = 0; j < number_parameters; j++)
			{
				b_star(j) = b(j) + deltab(j);
				const double f_plus = OptFunction(b_star);
				b_star(j) = b(j) - deltab(j);
				const double f_minus = OptFunction(b_star);

				grad[j] = (f_plus - f_minus) / (2.*deltab(j));

				b_star(j) = b(j);
			}
		}

		numberOfGradientEvaluations++;
	}
	
	numberOfFunctionEvaluations++;
	return f;
}

double OptFunction(const Eigen::VectorXd &b, unsigned int eval_nr)
{
    double fn_val;
    sim_iface_.SubstituteKineticParameters(b);

	if(input.optimization_setup().penalty_function())
		violated_uncertainty = sim_iface_.CheckKineticConstasts();
        
	if(violated_uncertainty){
		if (input.optimization_setup().objective_function_type() == "CurveMatching")
			fn_val = 1;
		else
			fn_val = 10000000;
	}
    else{
        sim_iface_.run();
        fn_val = sim_iface_.ComputeObjectiveFunction();
        
        if(eval_nr == 1)
            prev_fn_val = fn_val;
        
        if(prev_fn_val > fn_val) {
            prev_fn_val = fn_val;
            opti_kinetics_.WriteOptimizedMechanism();
            std::cout << " * Wrote optimized mechanism" << std::endl;
        }
    }
    return fn_val;
}
// #endif