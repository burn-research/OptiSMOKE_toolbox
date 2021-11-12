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

#include "Grammar_files/Grammar_OptiSMOKEpp.h"
#include "Grammar_files/Grammar_DakotaOptions.h"
#include "Grammar_files/Grammar_CurveMatchingOptions.h"
#include "Read_Input/Read_Input.h"
#include "ParallelLibrary.hpp"
#include "ProblemDescDB.hpp"
#include "LibraryEnvironment.hpp"
#include "DakotaModel.hpp"
#include "DakotaInterface.hpp"
#include "Dakota_Plugin.h"

//#include "CURVE_MATCHING/BasisFunction.h"
//#include "CURVE_MATCHING/Indexes.h"
//#include "CURVE_MATCHING/Spline.h"
//#include "CURVE_MATCHING/Utilities.h"

#ifdef HAVE_AMPL 
/** Floating-point initialization from AMPL: switch to 53-bit rounding
    if appropriate, to eliminate some cross-platform differences. */
extern "C" void fpinit_ASL(); 
#endif 

//#ifndef DAKOTA_HAVE_MPI
//#define MPI_COMM_WORLD 0
//#endif // not DAKOTA_HAVE_MPI


/// Run a Dakota LibraryEnvironment, mode 1: parsing an input file
void run_dakota_parse(const char* plugin_input_file);

void opensmoke_interface_plugin(Dakota::LibraryEnvironment& env,const char* plugin_input_file);

int main(int argc, char* argv[])
{

#ifdef HAVE_AMPL
  // Switch to 53-bit rounding if appropriate, to eliminate some
  // cross-platform differences.
  fpinit_ASL();	
#endif

  // whether running in parallel
//  bool parallel = Dakota::MPIManager::detect_parallel_launch(argc, argv);
  bool parallel = false;
  // Define MPI_DEBUG in dakota_global_defs.cpp to cause a hold here
  Dakota::mpi_debug_hold();

//#ifdef DAKOTA_HAVE_MPI
//  if (parallel)
//    MPI_Init(&argc, &argv); // initialize MPI
//#endif // DAKOTA_HAVE_MPI

  // Allow MPI to extract its command line arguments first in detect above,
  // then detect "-mixed" and dakota_input_file
  bool mixed_input = false;
  const char *plugin_input_file = NULL;
  plugin_input_file = argv[1];

  run_dakota_parse(plugin_input_file); // mode 1: parse

  // Note: Dakota objects created in above function calls need to go
  // out of scope prior to MPI_Finalize so that MPI code in
  // destructors works properly in library mode.

//#ifdef DAKOTA_HAVE_MPI
//  if (parallel)
//    MPI_Finalize(); // finalize MPI
//#endif // DAKOTA_HAVE_MPI

  return 0;
}

/** Simplest library case: this function parses from an input file to define the
    ProblemDescDB data. */
void run_dakota_parse(const char* plugin_input_file)// const char opensmoke_input_file)	
{
	// static const char dakota_options;
  // initialized and object of class Read_Input
	OpenSMOKE::Read_Input ObjectInput1;

  // read everything into this object
	ObjectInput1.ReadInfo(plugin_input_file,false);

	//Print out OptiSMOKE++ logo
	std::cout<<"//----------------------------------------------------------------------//"<<std::endl;
	std::cout<<"//     ____            _  ______ __  __  ____  _  ________              //"<<std::endl;
	std::cout<<"//    / __ \\       _  (_)/  ___ |  \\/  |/ __ \\| |/ /  ____|             //"<<std::endl;
	std::cout<<"//   | |  | |_ __ | |_ _ | (___ | \\  / | |  | | ' /| |__    _     _     //"<<std::endl;
	std::cout<<"//   | |  | | '_ \\|  _| |\\___  \\| |\\/| | |  | |  < |  __| _| |_ _| |_   //"<<std::endl;
	std::cout<<"//   | |__| | |_) | |_| |____)  | |  | | |__| | . \\| |___|_   _|_   _|  //"<<std::endl;
	std::cout<<"//    \\____/| .__/\\___|_|______/|_|  |_|\\____/|_|\\_\\______||_|   |_|    //"<<std::endl;
	std::cout<<"//          | |                                                         //"<<std::endl;
	std::cout<<"//          |_|                                                         //"<<std::endl;
	std::cout<<"//                                                                      //"<<std::endl;
	std::cout<<"//            Author: Magnus Fürst <magnus.furst@ulb.ac.be>             //"<<std::endl;
	std::cout<<"//                    Andrea Bertolino <andrea.bertolino@ulb.be>        //"<<std::endl;
	std::cout<<"//----------------------------------------------------------------------//"<<std::endl;
	std::cout<<""<<std::endl;
	std::cout<<"Input file: "<<plugin_input_file<<std::endl;
	std::cout<<""<<std::endl;
	// Prepare Dakota input string
	ObjectInput1.DakotaInputString();

  // Parse input and construct Dakota LibraryEnvironment, performing
  // input data checks
  Dakota::ProgramOptions opts;
  opts.input_string(ObjectInput1.dakota_options_string);

  // Defaults constructs the MPIManager, which assumes COMM_WORLD
  Dakota::LibraryEnvironment env(opts);

  if (env.mpi_manager().world_rank() == 0)
    Cout << "Library mode 1: run_dakota_parse()\n";

  // plug the client's interface (function evaluator) into the Dakota
  // environment; in serial case, demonstrate the simpler plugin method
  /*if (env.mpi_manager().mpirun_flag())
    parallel_interface_plugin(env);
  else*/
    opensmoke_interface_plugin(env,plugin_input_file);

  // Execute the environment
  env.execute();
  // Print out new kinetics.CKI file with the best parameter values
  // retrieve the final parameter values
  /*const Dakota::Variables& vars = env.variables_results(); 
  Dakota::RealVector best_parameters;
  best_parameters = vars.all_continuous_variables();
  std::vector<double> best_param_vec;
  for(int i = 0; i < best_parameters.length(); i++)
  {
	best_param_vec.push_back(best_parameters[i]);
  }
  if (ObjectInput1.CKI_File_Read)
  {
  	ObjectInput1.PrintFinalMechanism(best_param_vec);
  } else
  {
  	std::cout<<"Not able to write optimized parameters to a new CKI file, as neither @KineticsPreProcessor nor @NominalKineticsPreProcessor is specified!"<<std::endl;
  }*/

}


void opensmoke_interface_plugin(Dakota::LibraryEnvironment& env,const char* plugin_input_file)//const char opensmoke_input_file)
{
    std::string model_type(""); // demo: empty string will match any model type
    std::string interf_type("direct");
    std::string an_driver("plugin_opensmoke");
    Dakota::ProblemDescDB& problem_db = env.problem_description_db();
    // here it actually uses the DirectApplicInterface
    Dakota::Interface* serial_iface = new SIM::OpenSMOKEDirectApplicInterface(problem_db,plugin_input_file);
    bool plugged_in =
    env.plugin_interface(model_type, interf_type, an_driver, serial_iface);
    
    if (!plugged_in) {
        Cerr << "Error: no serial interface plugin performed.  Check "
        << "compatibility between parallel\n       configuration and "
        << "selected analysis_driver." << std::endl;
        Dakota::abort_handler(-1);
    }
}



