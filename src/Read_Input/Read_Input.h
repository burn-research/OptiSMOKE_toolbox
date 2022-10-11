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

/* include guards they check if Read_Input_H is defined, if not, they do. */
#ifndef Read_Input_H
#define Read_Input_H
// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"
// Thermodynamics
#include "kernel/thermo/Species.h"
#include "kernel/thermo/ThermoPolicy_CHEMKIN.h"
#include "kernel/thermo/ThermoReader.h"
#include "kernel/thermo/ThermoReaderPolicy_CHEMKIN.h"
// Kinetics
#include "kernel/kinetics/ReactionPolicy_CHEMKIN.h"
// Preprocessing
#include "preprocessing/PreProcessorSpecies.h"
#include "preprocessing/PreProcessorKinetics.h"
#include "preprocessing/PreProcessorKineticsPolicy_CHEMKIN.h"
#include "preprocessing/PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport.h"
#include "preprocessing/PreProcessorSpeciesPolicy_CHEMKIN_WithTransport.h"
// CHEMKIN maps
#include "maps/Maps_CHEMKIN"
// Reactor utilities
#include "idealreactors/utilities/Utilities"
// CurveMatching random header
#include <random>
#include <chrono>
// CURVE_MATCHING classes
//#include "../CURVE_MATCHING/Indexes.h"
//#include "../CURVE_MATCHING/GlobalVariables.h"
//#include "../CURVE_MATCHING/BasisFunction.h"
//#include "../CURVE_MATCHING/Indexes.h"
//#include "../CURVE_MATCHING/Spline.h"
//#include "../CURVE_MATCHING/Utilities.h"

//#include <boost/algorithm/string.hpp>
//#include "preprocessing/PreProcessor_CHEMKIN"
//#include <experimental/filesystem>
namespace OpenSMOKE
{
	class ThermodynamicsMap_CHEMKIN;	
	class KineticsMap_CHEMKIN;
	class TransporPropertiestMap_CHEMKIN;	
	//class PreProcessorKinetics_CHEMKIN;
	class Read_Input
	{
	// Global parameters
	public:

	bool Debug_print_out;
	bool Debug_Sim = false;
	

	bool Read_transport = true;

	double default_sigma=2;
	std::string boundaries_method;	
	//std::string type_of_reactor;
	int N_of_plugflow_datasets;
	int N_of_batch_datasets;
	int N_of_psr_datasets;
	int N_of_laminar_flame_datasets;
	int N_of_direct_measurement_dataset = 0;

	std::string distribution;
	// for constraints!
	bool udc_bool = false;
	std::vector<std::vector<double>> ud_constraints;

	// AB // created these 3 variables to store contents of the header file of experimental data files
	// Now i will add a third structure (v of v of doubles) with the standard deviation assuming 2 sigma;
	std::vector<std::string> QoI;
	std::vector<std::string> type_of_reactor;
	std::vector<std::string> type_KTP;
	std::vector<double> value_KTP;

	// UPDATED
	std::vector<std::vector<std::string>> what_2_calc;
	std::vector<std::string> Sigma_vector;
	// UPDATED
	std::vector<std::vector<std::vector<double>>> standard_deviations;

	std::vector<std::string> tau_calc_type;
	std::vector<std::string> species_of_interest;
	std::string Objective_Function;
	int  SigmaExpDistribution;
	int AcceptedSigmaKDistribution;
	std::vector<std::vector<std::string>> list_of_opensmoke_input_files;
	std::vector<std::string> experimental_data_files;
	std::vector<std::vector<std::string>> Order_of_tau_calc;
	std::vector<std::vector<std::vector<std::vector<double>>>> Exp_data;
	//std::vector<std::vector<std::vector<double>>> bootstrapExp;
	int number_of_parameters;
	int number_of_threads;
	std::string file_extension;
	bool penalty_function;	
	bool UseBootStrap;
	int numberOfBootstrapVariations;
	bool print_indexes;
	bool print_splines;
	bool print_bootstrap;
	std::vector<unsigned int> indices_of_falloff_reactions;
	std::vector<std::vector<unsigned int>> falloff_indices_of_thirdbody_species;

	// EPLR
	std::vector<int> list_of_target_EPLR;
	std::vector<std::string>  list_of_bath_gases_EPLR;
	std::vector<double>       list_of_uncertainty_factors_EPLR;
	// 
	std::vector<std::string> list_of_nominal_lnA_EPLR;
	std::vector<std::string> list_of_min_lnA_EPLR;
	std::vector<std::string> list_of_max_lnA_EPLR;
	//
	std::vector<std::string> list_of_nominal_ER_EPLR;
	std::vector<std::string> list_of_min_ER_EPLR;
	std::vector<std::string> list_of_max_ER_EPLR;
	// 
	std::vector<std::string> list_of_nominal_Beta_EPLR;
	std::vector<std::string> list_of_min_Beta_EPLR;
	std::vector<std::string> list_of_max_Beta_EPLR;
	// EPLR
	
	// EXTENDED PLOG
	// Indices of extended plog
	std::vector<unsigned int> indices_of_extendedplogs;
	std::vector<unsigned int> nominal_indices_of_extendedplogs;
	// Indices of extended plog for optimization
	std::vector<int> 	list_of_target_extplog;
	std::vector<double> list_of_uncertainty_factors_extplog;

	// extended PLOG for Third Bodies
	std::vector<unsigned int> indices_of_extendedplogs_opt;
	std::vector<unsigned int> nominal_indices_of_extendedplogs_opt;
	std::vector<int> list_of_target_extended_plog_reactions;
	std::vector<std::string> list_of_target_extended_plog_species;
	std::vector<double> list_of_min_tb_extplog;
	std::vector<double> list_of_max_tb_extplog;
	// STRINGS of nominal, min and max value have to be defined here, as they are used in different scopes //

	std::vector<std::string> list_of_nominal_TB_ExtPLOG;
	std::vector<std::string> list_of_min_TB_ExtPLOG;
	std::vector<std::string> list_of_max_TB_ExtPLOG;
	// 
	std::vector<std::string> list_of_nominal_lnA_ext_plog_coefficients;
	std::vector<std::string> list_of_min_lnA_ext_plog_coefficients;
	std::vector<std::string> list_of_max_lnA_ext_plog_coefficients;
	//
	std::vector<std::string> list_of_nominal_ER_ext_plog_coefficients;
	std::vector<std::string> list_of_min_ER_ext_plog_coefficients;
	std::vector<std::string> list_of_max_ER_ext_plog_coefficients;
	// 
	std::vector<std::string> list_of_nominal_Beta_ext_plog_coefficients;
	std::vector<std::string> list_of_min_Beta_ext_plog_coefficients;
	std::vector<std::string> list_of_max_Beta_ext_plog_coefficients;


	// DIRECT REACTIONS
	std::vector<int> list_of_target_lnA;
	std::vector<int> list_of_target_Beta;
	std::vector<int> list_of_target_E_over_R;

	// INF REACTIONS
	std::vector<int> list_of_target_lnA_inf;
	std::vector<int> list_of_target_Beta_inf;
	std::vector<int> list_of_target_E_over_R_inf;

	// INF REACTIONS
	std::vector<int> list_of_target_thirdbody_reactions;
	std::vector<std::string> list_of_target_thirdbody_species;


	// AB // PLOG
	std::vector<unsigned int> indices_of_classic_plogs;
	std::vector<unsigned int> nominal_indices_of_classic_plogs;

	// Independent component analysis (ICA) for direct reactions
	// for input file
	std::vector<int> 			direct_reactions_indices_ica;
	std::vector<double> 		abs_max_direct_reactions_ica;
	std::vector<double> 		abs_min_direct_reactions_ica;
	std::vector<double> 		abs_init_direct_reactions_ica;

	// for code
	std::vector<std::string> 	mixing_matrices_files_ica	;
	std::vector<std::vector<std::vector<double>>> mixing_matrices_ica;

	// PCA direct reactions -- AB
	std::vector<int> number_of_eigenvector_for_reaction;

	std::vector<int> pca_direct_reactions;
	std::vector<double> pca_maxabs_direct_reactions;
	std::vector<double> pca_minabs_direct_reactions;
	std::vector<std::string> pca_files_direct_reactions;
	

	std::vector<std::vector<double>> pca_scaling_direct_reactions;
	std::vector<std::vector<std::vector<double>>> pca_eigenvectors_direct_reactions;
	std::vector<std::vector<double>> pca_centering_direct_reactions;

	//PCA Pinf reactions -- AB
	std::vector<int> pca_pinf_reactions;
	std::vector<double> pca_maxabs_pinf_reactions;
	std::vector<double> pca_minabs_pinf_reactions;
	std::vector<std::string> pca_files_pinf_reactions;

	std::vector<std::vector<double>> pca_scaling_pinf_reactions;
	std::vector<std::vector<std::vector<double>>> pca_eigenvectors_pinf_reactions;
	std::vector<std::vector<double>> pca_centering_pinf_reactions;

	//PCA PLOG reactions -- AB
	std::vector<int> pca_plog_reactions;
	std::vector<double> pca_maxabs_plog_reactions;
	std::vector<double> pca_minabs_plog_reactions;
	std::vector<std::string> pca_files_plog_reactions;

	std::vector<std::vector<double>> pca_scaling_plog_reactions;
	std::vector<std::vector<std::vector<double>>> pca_eigenvectors_plog_reactions;
	std::vector<std::vector<double>> pca_centering_plog_reactions;


	// PLOG -- AB
	std::vector<int> list_of_target_classic_plog_reactions;
	// AB // List of uncertainties for plogs
	std::vector<double> list_of_uncertainty_factors_classic_plog;
	// STRINGS of nominal, min and max value have to be defined here, as they are used in different scopes //
	std::vector<std::string> list_of_nominal_lnA_classic_plog_coefficients;
	std::vector<std::string> list_of_min_lnA_classic_plog_coefficients;
	std::vector<std::string> list_of_max_lnA_classic_plog_coefficients;

	std::vector<std::string> list_of_nominal_ER_classic_plog_coefficients;
	std::vector<std::string> list_of_min_ER_classic_plog_coefficients;
	std::vector<std::string> list_of_max_ER_classic_plog_coefficients;

	std::vector<std::string> list_of_nominal_Beta_classic_plog_coefficients;
	std::vector<std::string> list_of_min_Beta_classic_plog_coefficients;
	std::vector<std::string> list_of_max_Beta_classic_plog_coefficients;
	// PLOG -- AB

	// this is for scaled parametrization
	std::vector<std::string> list_of_nom_abs_A_scaled;
	std::vector<std::string> list_of_min_abs_A_scaled;
	std::vector<std::string> list_of_max_abs_A_scaled;


	std::vector<std::string> list_of_initial_lnA;
	std::vector<std::string> list_of_initial_Beta;
	std::vector<std::string> list_of_initial_E_over_R;
	std::vector<std::string> list_of_initial_lnA_inf;
	std::vector<std::string> list_of_initial_Beta_inf;
	std::vector<std::string> list_of_initial_E_over_R_inf;
	std::vector<std::string> list_of_initial_thirdbody_eff;

	std::vector<std::string> list_of_min_abs_lnA;
	std::vector<std::string> list_of_max_abs_lnA;
	std::vector<double> list_of_min_rel_lnA;
	std::vector<double> list_of_max_rel_lnA;

	std::vector<std::string> list_of_min_abs_Beta;
	std::vector<std::string> list_of_max_abs_Beta;
	std::vector<double> list_of_min_rel_Beta;
	std::vector<double> list_of_max_rel_Beta;

	std::vector<std::string> list_of_min_abs_E_over_R;
	std::vector<std::string> list_of_max_abs_E_over_R;
	std::vector<double> list_of_min_rel_E_over_R;
	std::vector<double> list_of_max_rel_E_over_R;

	std::vector<std::string> list_of_min_abs_lnA_inf;
	std::vector<std::string> list_of_max_abs_lnA_inf;
	std::vector<double> list_of_min_rel_lnA_inf;
	std::vector<double> list_of_max_rel_lnA_inf;

	std::vector<std::string> list_of_min_abs_Beta_inf;
	std::vector<std::string> list_of_max_abs_Beta_inf;
	std::vector<double> list_of_min_rel_Beta_inf;
	std::vector<double> list_of_max_rel_Beta_inf;

	std::vector<std::string> list_of_min_abs_E_over_R_inf;
	std::vector<std::string> list_of_max_abs_E_over_R_inf;
	std::vector<double> list_of_min_rel_E_over_R_inf;
	std::vector<double> list_of_max_rel_E_over_R_inf;

	std::vector<std::string> list_of_min_abs_thirdbody_eff;
	std::vector<std::string> list_of_max_abs_thirdbody_eff;
	std::vector<double> list_of_min_rel_thirdbody_eff;
	std::vector<double> list_of_max_rel_thirdbody_eff;
	
	std::vector<int> list_of_target_uncertainty_factors;
	std::vector<double> list_of_uncertainty_factors;

	std::vector<int> list_of_target_uncertainty_factors_inf;
	std::vector<double> list_of_uncertainty_factors_inf;


	boost::filesystem::path path_kinetics_output; 
	boost::filesystem::path path_nominal_kinetics_output; 

	std::string tabular_data_file;
	std::string method;
	std::string string_max_iterations;
	std::string string_max_function_evaluations;
	std::string string_convergence_tolerance;
	std::string string_solution_target;
	std::string string_seed;
	std::string string_population_size;
	std::string fitness_type;
	std::string mutation_type;
	std::string mutation_rate;
	std::string crossover_type;
	std::string crossover_rate;
	std::string replacement_type;
	std::string string_division;
	std::string max_boxsize_limit;
	std::string min_boxsize_limit;
	bool gradient_option;
	// String containing diverse input given to dakota
	std::vector<std::string> diverse_dakota_input;

	// For Dakota input
	std::string dakota_options_string;
	ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML;	
	TransportPropertiesMap_CHEMKIN* transportMapXML;	
	KineticsMap_CHEMKIN* kineticsMapXML;
	ThermodynamicsMap_CHEMKIN* nominalthermodynamicsMapXML;	
	KineticsMap_CHEMKIN* nominalkineticsMapXML;
	TransportPropertiesMap_CHEMKIN* nominaltransportMapXML;	
	//PreProcessorKinetics_CHEMKIN* preprocessor_kinetics;
	boost::filesystem::path path_input_thermodynamics;
	boost::filesystem::path path_input_kinetics;
	std::string path_folder;
	bool CKI_File_Read;
	// Functions

	typedef OpenSMOKE::PreProcessorKinetics< OpenSMOKE::PreProcessorKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_CHEMKIN> > PreProcessorKinetics_CHEMKIN;
 	PreProcessorKinetics_CHEMKIN* preprocessor_kinetics;

	void Initialize_the_pre_processor();
	void ReadInfo(const char* plugin_input_file, bool print_out);
	void ReadKinetics(bool print_out);
	void ReadNominalKinetics();
	void ReadExpData_xml();
	void ReadExpData();
	void ReadExpData_csv();
	void ParameterLimits();
	void DakotaInputString();
	//void BootStrapping_exp_data();
	//void PrintFinalMechanism(std::vector<double> best_parameters);
	void PrintFinalMechanism();
	// AB // Curve Matching
	/* Initializing some strings */
    std::string possibleNegativeOrdinatesString;
    std::string lineUpMaximaString = "true";
    std::string useSumOfIndexesForAlignmentString = "yes";
    std::string useIndexesToChooseExpNodesString;

	/* Number of plausible versions of the experimental data generated during the
	bootstrap procedure. numberOfBootstrapVariations = 1 : no bootstrap */

	// I'll put it here because 
	//int numberOfBootstrapVariations = 20;

	/* Specifies whether negative segments on the y-axis are admissible for the
	splines or whether they should be replaced with straight lines with ordinate 0
	*/
	bool possibleNegativeOrdinates = false;

	/* Specifies whether to calculate the exponential (base e) of the input
	ordinates before computing the splines */
	bool logScale = false;

	// initializes the vector of splines
	
	
	///////////////////////////////////////REACTIONS CLASSES////////////////////////////////////////
	
	bool Optimization4Classes = false;
	bool ScalingReactionClasses = false;
	boost::filesystem::path ReactionClassesPath;
	int numberOfReactionClasses;
	std::vector<std::vector<int>> matrixOfReactionIndex; // matrice dove mi salvo gli indici
    std::vector<std::vector<double>> matrixOfUnceratintyFactors; // matrice dove gli unc factor
    std::vector<std::vector<std::string>> matrixOfQOI; // matrice delle QOI forse la tolgo poi vediamo

	void ReadReactionClassesDefinition(boost::filesystem::path ReactionClassFile);
	};
}
#include "Read_Input.hpp"
#endif
