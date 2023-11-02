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

#ifdef OPTISMOKE_USE_MPI
#include <mpi.h>
#endif

// Standard Library
#include <string>
#include <iostream>
#include <numeric>
#include <random>
#include <chrono>
#include <algorithm>

// OpenSMOKE 
#include "OpenSMOKEpp"

// Thermodynamics
#include "kernel/thermo/Species.h"
#include "kernel/thermo/ThermoPolicy_CHEMKIN.h"
#include "kernel/thermo/ThermoReader.h"
#include "kernel/thermo/ThermoReaderPolicy_CHEMKIN.h"

// Transport
#include "kernel/transport/TransportPolicy_CHEMKIN.h"
#include "kernel/transport/TransportReader.h"
#include "kernel/transport/TransportReaderPolicy_CHEMKIN.h"

// Kinetics
#include "kernel/kinetics/ReactionPolicy_CHEMKIN.h"

// Preprocessing
#include "preprocessing/PreProcessorSpecies.h"
#include "preprocessing/PreProcessorKinetics.h"
#include "preprocessing/PreProcessorKineticsPolicy_CHEMKIN.h"
#include "preprocessing/PreProcessorSpeciesPolicy_CHEMKIN_WithTransport.h"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/TransportPropertiesMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"

//Typedefs
typedef OpenSMOKE::Species<OpenSMOKE::ThermoPolicy_CHEMKIN, OpenSMOKE::TransportPolicy_CHEMKIN> SpeciesCHEMKIN;
typedef OpenSMOKE::PreProcessorSpecies<OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<SpeciesCHEMKIN>> PreProcessorSpecies_CHEMKIN_WithoutTransport;
typedef OpenSMOKE::PreProcessorKinetics<OpenSMOKE::PreProcessorKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_CHEMKIN>> PreProcessorKinetics_CHEMKIN;
typedef OpenSMOKE::ThermoReader<OpenSMOKE::ThermoReaderPolicy_CHEMKIN<OpenSMOKE::ThermoPolicy_CHEMKIN>> ThermoReader_CHEMKIN;

// Boost Library
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
namespace fs = boost::filesystem;
namespace po = boost::program_options;

const double UNFEASIBLE_BIG_NUMBER = 1.e16;

// Dakota Library
#ifdef OPTISMOKE_USE_MPI
#define DAKOTA_HAVE_MPI
#include <ParallelLibrary.hpp>
#include <PluginParallelDirectApplicInterface.hpp>
#endif

#include "ProblemDescDB.hpp"
#include "LibraryEnvironment.hpp"
#include "DakotaModel.hpp"
#include "DakotaInterface.hpp"
#include "DakotaResponse.hpp"
#include "ParamResponsePair.hpp"
#include "DirectApplicInterface.hpp"

#ifdef OPTISMOKE_USE_NLOPT
#include <nlopt.hpp>
double NLOptFunction(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);
double OptFunction(const std::vector<double> &b, unsigned int fn_val);
#endif

// Curve Matching
#include "curve_matching/curve_matching.h"

// Header files
#include "utilities/OptiSMOKEUtilities"
#include "grammar/grammar.h"
#include "options/options.h"
#include "ideal_reactors/ideal_reactors.h"
#include "DataManager.h"
#include "InputManager.h"
#include "OptimizedKinetics.h"
#include "SerialDakotaInterface.h"
#include "SimulationsInterface.h"

// Run a Dakota LibraryEnvironment, mode 1: parsing an input file
void run_dakota_parse(const char* plugin_input_file, bool echo_dakota_string, int rank);

void opensmoke_interface_plugin(Dakota::LibraryEnvironment& env);

OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
OptiSMOKE::InputManager input(dictionaries);

#ifdef OPTISMOKE_USE_NLOPT
OptiSMOKE::SimulationsInterface* sim_iface_;
OptiSMOKE::OptimizedKinetics* opti_kinetics_;
unsigned int numberOfGradientEvaluations;
unsigned int numberOfFunctionEvaluations;
bool violated_uncertainty;
double prev_fn_val;
std::ofstream fOut;
# endif

#ifdef OPTISMOKE_USE_MPI
// void parallel_interface_plugin(Dakota::LibraryEnvironment& env);
#endif // OPTISMOKE_USE_MPI
