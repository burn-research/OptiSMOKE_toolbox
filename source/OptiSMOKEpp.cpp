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

    input->SetInputOptions(argc, argv);
    input->ReadDictionary();
    input->ReadExperimentalDataFiles();

    if(input->optimization_library() == "dakota"){
        input->DakotaInputString();
        dakota_input_string = input->dakota_input_string().c_str();
        run_dakota_parse(dakota_input_string, input->dakota_options().echo_dakota_string());

        // Note: Dakota objects created in above function calls need to go
        // out of scope prior to MPI_Finalize so that MPI code in
        // destructors works properly in library mode.

        #ifdef DAKOTA_HAVE_MPI
        if (parallel)
            MPI_Finalize(); // finalize MPI
        #endif // DAKOTA_HAVE_MPI
    
        return 0;
    }
    else if (input->optimization_library() == "nlopt"){
        OptiSMOKE::FatalErrorMessage("Not yet come up with the implementation!");
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
