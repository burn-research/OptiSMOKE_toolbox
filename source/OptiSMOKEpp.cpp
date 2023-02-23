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
|            Authors: Magnus Fürst <magnus.furst@ulb.ac.be>               |
|                     Andrea Bertolino <andrea.bertolino@ulb.be>          |
|                     Timoteo Dinelli <timoteo.dinelli@polimi.it>         |
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

int main(int argc, char* argv[])
{

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

    // Allow MPI to extract its command line arguments first in detect above,
    // then detect "-mixed" and dakota_input_file
    OptiSMOKE::OptiSMOKE_logo("OptiSMOKE++", "M. Furst, A. Bertolino, T. Dinelli");
    
    bool mixed_input = false;
    const char *plugin_input_file = NULL;
    plugin_input_file = argv[1];

    // Reading OptiSMOKE Input file
    OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
    OptiSMOKE::InputManager input(dictionaries);

    input.SetInputOptions(argc, argv);
    input.ReadDictionary();

    run_dakota_parse(plugin_input_file); // mode 1: parse

    // Note: Dakota objects created in above function calls need to go
    // out of scope prior to MPI_Finalize so that MPI code in
    // destructors works properly in library mode.

    #ifdef DAKOTA_HAVE_MPI
    if (parallel)
        MPI_Finalize(); // finalize MPI
    #endif // DAKOTA_HAVE_MPI

    return 0;
}

void run_dakota_parse(const char* plugin_input_file)
{

    // Parse input and construct Dakota LibraryEnvironment, performing
    // input data checks
    Dakota::ProgramOptions opts;
    // opts.input_string(ObjectInput1.dakota_options_string);

    // Defaults constructs the MPIManager, which assumes COMM_WORLD
    
    Dakota::LibraryEnvironment env(opts);

    if (env.mpi_manager().world_rank() == 0){
        Cout << "Library mode 1: run_dakota_parse()\n";
    }
    
    // plug the client's interface (function evaluator) into the Dakota
    // environment; in serial case, demonstrate the simpler plugin method
    /*if (env.mpi_manager().mpirun_flag())
    parallel_interface_plugin(env);
    else */
    //opensmoke_interface_plugin(env,plugin_input_file);

    // Execute the environment
    //env.execute();

}

/*
void opensmoke_interface_plugin(Dakota::LibraryEnvironment& env,const char* plugin_input_file)//const char opensmoke_input_file)
{
    std::string model_type(""); // demo: empty string will match any model type
    std::string interf_type("direct");
    std::string an_driver("plugin_opensmoke");
        
    Dakota::ProblemDescDB& problem_db = env.problem_description_db();
        
    // here it actually uses the DirectApplicInterface
    
    // This is the old way of instace connection with Dakota
    // Dakota::Interface* serial_iface = new SIM::OpenSMOKEDirectApplicInterface(problem_db,plugin_input_file);
    
    // This is the new way with the shared pointer
    std::shared_ptr<Dakota::Interface> serial_iface =  
        std::make_shared<SIM::OpenSMOKEDirectApplicInterface>(problem_db, plugin_input_file);
    bool plugged_in = env.plugin_interface(model_type, interf_type, an_driver, serial_iface);
        
    if (!plugged_in) {
        Cerr << "Error: no serial interface plugin performed.  Check "
        << "compatibility between parallel\n       configuration and "
        << "selected analysis_driver." << std::endl;
        Dakota::abort_handler(-1);
    }
}*/
