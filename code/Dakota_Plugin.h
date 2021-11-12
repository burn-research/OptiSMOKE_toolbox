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
|   Copyright (C) 2020 by Magnus Fürst                                    |
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

// ciao
#ifndef DAKOTA_PLUGIN_H
#define DAKOTA_PLUGIN_H

// Include OpenMP Header file
//#if defined(_OPENMP)
#include <omp.h>
//#endif

#include "DakotaResponse.hpp"
#include "DakotaConstraints.hpp"
#include "ParamResponsePair.hpp"
#include "ProblemDescDB.hpp"
#include "ParallelLibrary.hpp"

#include "idealreactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"
#include "utilities/cema/OnTheFlyCEMA.h"
#include "utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h"
#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"
#include "Grammar_files/Grammar_PlugFlowReactor.h"
#include "idealreactors/plugflow/PlugFlowReactor"
#include "PFR/PlugFlowReactor_Plugin.h"
#include "Grammar_files/Grammar_BatchReactor.h"
#include "idealreactors/batch/BatchReactor"
#include "Batch/BatchReactor_Plugin.h"
#include "Grammar_files/Grammar_PerfectlyStirredReactor.h"
#include "idealreactors/psr/PerfectlyStirredReactor"
#include "PSR/PerfectlyStirredReactor_Plugin.h"
#include "DirectApplicInterface.hpp"
// 1D grid
#include "utilities/grids/adaptive/Grid1D.h"
#include "utilities/grids/adaptive/Grammar_Grid1D.h"
#include "utilities/grids/adaptive/Adapter_Grid1D.h"	

#include <memory>
// Hybrid Method of Moments
#include "utilities/soot/hmom/HMOM.h"
#include "idealreactors/utilities/Grammar_LewisNumbers.h"
#include "Flame1D/utilities/Utilities.h"
#include "Flame1D/OpenSMOKE_PremixedLaminarFlame1D.h"
#include "Grammar_files/Grammar_PremixedLaminarFlame1D.h"
#include "Flame1D/Flame1D_Plugin.h"
//OpenSMOKE::OpenSMOKE_PremixedLaminarFlame1D* flame_premixed;
#include "Flame1D/OpenSMOKE_PremixedLaminarFlame1D.hpp"
#include "Flame1D/Flame1D_Plugin.hpp"
#include "CURVE_MATCHING/Indexes.h"

namespace OpenSMOKE{
	class PlugFlowReactor_Plugin;
	class BatchReactor_Plugin;
	class PerfectlyStirredReactor_Plugin;
  	class Flame1D_Plugin;
}
/// A sample namespace for derived classes that use assign_rep() to 
/// plug facilities into DAKOTA.

/** A typical use of plug-ins with assign_rep() is to publish a 
    simulation interface for use in library mode  See \ref DakLibrary
    for more information. */

namespace SIM {

/// Sample derived interface class for testing serial simulator
/// plug-ins using assign_rep().

/** The plug-in SerialDirectApplicInterface resides in namespace SIM
    and uses a copy of rosenbrock() to perform serial parameter to
    response mappings.  It is used to demonstrate plugging in a serial
    direct analysis driver into Dakota in library mode.  Test input
    files can then use an analysis_driver of "plugin_rosenbrock". */

//SIM::OpenSMOKEDirectApplicInterface.laminar_flames

class OpenSMOKEDirectApplicInterface: public Dakota::DirectApplicInterface
{
public:

  //
  //- Heading: Constructor and destructor
  //
	
	// Check that the kinetic parameters are changed correctly
	bool Debug_ChangeParam;
	

	// Object containing input file information
	OpenSMOKE::Read_Input ObjectInput2;

	// Reactor objects
  // This is one of the most important part, where i had objects of the class referred to each canonical reactor;
	std::vector<OpenSMOKE::PlugFlowReactor_Plugin*> plugflow_reactors;
	std::vector<OpenSMOKE::BatchReactor_Plugin*> batch_reactors;
	std::vector<OpenSMOKE::PerfectlyStirredReactor_Plugin*> psr_reactors;
  	//std::vector<std::vector<std::shared_ptr<OpenSMOKE::Flame1D_Plugin>>> laminar_flames;
	//std::shared_ptr<OpenSMOKE::Flame1D_Plugin> laminar_flames;
	
	// Functions for calculating the objective function (Least Squared Residuals and Normalized)
	double Calculate_ObjFunc(std::vector<std::vector<std::vector<double>>> Sim_values);

	// Function for controlling that proposed parameters respect uncertainty limits of the rate coefficients
	bool Check_k(Dakota::Real& fn_val);

	// Functions for changing kinetic parameters
  void ChangelnA(int i, double lnA_new);
  void ChangeBeta(int i, double Beta_new);
  void ChangeE_over_R(int i, double E_over_R_new);
  void ChangelnA_inf(int i, double lnA_inf_new);
  void ChangeBeta_inf(int i, double Beta_inf_new);
  void ChangeE_over_R_inf(int i, double E_over_R_inf_new);
	void ChangeThirdBody_Eff(int i, std::string Species, double ThirdBody_Eff_new);

  void ChangelnA_EPLR(int i, double coefficient_new, std::string bath_gas);
  void Change_ER_EPLR(int i, double coefficient_new, std::string bath_gas);
  void Change_Beta_EPLR(int i, double coefficient_new, std::string bath_gas);
  
  void ChangelnA_classic_PLOG(int i, double lnA_PLOG_coefficient_new);
  void Change_ER_classic_PLOG(int i, double lnA_PLOG_coefficient_new);
  void Change_Beta_classic_PLOG(int i, double lnA_PLOG_coefficient_new);

  void ChangelnA_ExtPLOG(int i, double lnA_PLOG_coefficient_new, std::vector<int> list_of_reactions, std::vector<std::string> list_of_species);
  void Change_ER_ExtPLOG(int i, double lnA_PLOG_coefficient_new, std::vector<int> list_of_reactions, std::vector<std::string> list_of_species);
  void Change_Beta_ExtPLOG(int i, double lnA_PLOG_coefficient_new, std::vector<int> list_of_reactions, std::vector<std::string> list_of_species);
  void Change_ExtPLOG_TB(int i, std::string Species, double new_value);
  int Extract_SpeciesPos_ExtPLOG(int extPLOG_pos, std::string species);

  // AB // Function to reconstruct sampled parameters according to eigenvectors, scaling and centering factors obtained in a pre-processing off-line step
  std::vector <double> reconstuction_pca(int i, std::vector<double> &z_parameters, std::vector<std::vector<std::vector<double>>> &eigenvectors, std::vector<std::vector<double>> &scaling, std::vector<std::vector<double>> &centering);
  
  // Function to mix the source signals from ICA transformation 
  std::vector <double> ica_mixing(int i, std::vector<double> &source_parameters, std::vector<std::vector<double>> &mixing_matrix);
  
	// Setup function
	void Setup_Plugin();

	// Number of input files for opensmoke simulations
	int nInputFiles;

	// Boolean for checking if uncertainty limits of the rate coefficient has been violated
	bool violated_uncertainty;


	// Nominal and limits of the rate constant
	std::vector<std::vector<double>> k_upper;
	std::vector<std::vector<double>> k_lower;
	std::vector<std::vector<double>> k_upper_inf;
	std::vector<std::vector<double>> k_lower_inf;
  //Constraints for PLOGs
	std::vector<std::vector<std::vector<double>>> k_UB_classic_plog;
  std::vector<std::vector<std::vector<double>>> k_LB_classic_plog;
  //Constraints for ExtPLOGs
	std::vector<std::vector<std::vector<double>>> k_UB_ext_plog;
  std::vector<std::vector<std::vector<double>>> k_LB_ext_plog;
  //Constraints for EPLR
	std::vector<std::vector<std::vector<double>>> k_UB_EPLR;
  std::vector<std::vector<std::vector<double>>> k_LB_EPLR;
 
	// Current evaluation number
	int eval_nr;
	double prev_fn_val;
	// String with name of input file
 	const char* plugin_input_file_;

  /// constructor
  OpenSMOKEDirectApplicInterface(const Dakota::ProblemDescDB& problem_db,const char* plugin_input_file);
  /// destructor
  ~OpenSMOKEDirectApplicInterface();

  bool logScale = false;
  std::vector<Spline> splines_Sim;
  std::vector<std::vector<std::vector<Indexes>>> indexes;
  std::vector<std::vector<std::vector<Spline>>> splinesExp;
  std::vector<std::vector<std::vector<std::vector<double>>>> bootstrapExp;

  void BootStrapping_exp_data(std::vector<std::vector<std::vector<std::vector<double>>>> Exp_data);
    
protected:

  //
  //- Heading: Virtual function redefinitions
  //

  // execute the input filter portion of a direct evaluation invocation
  //int derived_map_if(const Dakota::String& if_name);

  /// execute an analysis code portion of a direct evaluation invocation
  int derived_map_ac(const Dakota::String& ac_name);

  // execute the output filter portion of a direct evaluation invocation
  //int derived_map_of(const Dakota::String& of_name);

  /// no-op hides base error; job batching occurs within
  /// wait_local_evaluations()
  void derived_map_asynch(const Dakota::ParamResponsePair& pair);

  /// evaluate the batch of jobs contained in prp_queue
  void wait_local_evaluations(Dakota::PRPQueue& prp_queue);
  /// invokes wait_local_evaluations() (no special nowait support)
  void test_local_evaluations(Dakota::PRPQueue& prp_queue);

  /// no-op hides default run-time error checks at DirectApplicInterface level
  void set_communicators_checks(int max_eval_concurrency);

private:

  //
  //- Heading: Convenience functions
  //

  /// Rosenbrock plug-in test function
  /*int rosenbrock(const Dakota::RealVector& c_vars, short asv,
		 Dakota::Real& fn_val, Dakota::RealVector& fn_grad,
		 Dakota::RealSymMatrix& fn_hess);*/
    
  /// OpenSMOKE plug-in
    int opensmoke_interface(const Dakota::RealVector& c_vars, short asv,Dakota::Real& fn_val);//,Dakota::RealVector& constr_val);

  //
  //- Heading: Data
  //

};

// Constructor
// here it runs the setup_plugin function.
inline OpenSMOKEDirectApplicInterface::
OpenSMOKEDirectApplicInterface(const Dakota::ProblemDescDB& problem_db,const char* plugin_input_file):
  Dakota::DirectApplicInterface(problem_db), plugin_input_file_(plugin_input_file)
{Setup_Plugin(); }


inline OpenSMOKEDirectApplicInterface::~OpenSMOKEDirectApplicInterface()
{ /* Virtual destructor handles referenceCount at Interface level. */ }


inline void OpenSMOKEDirectApplicInterface::
derived_map_asynch(const Dakota::ParamResponsePair& pair)
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
inline void OpenSMOKEDirectApplicInterface::
test_local_evaluations(Dakota::PRPQueue& prp_queue)
{ wait_local_evaluations(prp_queue); }


// Hide default run-time error checks at DirectApplicInterface level
inline void OpenSMOKEDirectApplicInterface::
set_communicators_checks(int max_eval_concurrency)
{ }

} // namespace Dakota
#include "Dakota_Plugin.hpp"
#endif
