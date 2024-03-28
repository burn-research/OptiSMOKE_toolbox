/*-----------------------------------------------------------------------*\
|                                                                         |
|       ____            _  ______ __  __  ____  _  ________               |
|      / __ \       _  (_)/  ___ |  \/  |/ __ \| |/ /  ____|              |
|     | |  | |_ __ | |_ _|  (___ | \  / | |  | | ' /| |__    _     _      |
|     | |  | | '_ \|  _| |\___  \| |\/| | |  | |  < |  __| _| |_ _| |_    |
|     | |__| | |_) | |_| |____)  | |  | | |__| | . \| |___|_   _|_   _|   |
|      \____/| .__/\___|_|______/|_|  |_|\____/|_|\_\______||_|   |_|     |
|            | |                                                          |
|            |_|                                                          |
|                                                                         |
\*-----------------------------------------------------------------------*/

// Standard Library
#include <string>
#include <iostream>
#include <numeric>
#include <random>
#include <chrono>
#include <algorithm>

using std::string;
using std::vector;

// OpenSMOKE
#include <OpenSMOKEpp>

// Thermodynamics
#include <kernel/thermo/Species.h>
#include <kernel/thermo/ThermoPolicy_CHEMKIN.h>
#include <kernel/thermo/ThermoReader.h>
#include <kernel/thermo/ThermoReaderPolicy_CHEMKIN.h>

// Transport
#include <kernel/transport/TransportPolicy_CHEMKIN.h>
#include <kernel/transport/TransportReader.h>
#include <kernel/transport/TransportReaderPolicy_CHEMKIN.h>

// Kinetics
#include <kernel/kinetics/ReactionPolicy_CHEMKIN.h>

// Preprocessing
#include <preprocessing/PreProcessorSpecies.h>
#include <preprocessing/PreProcessorKinetics.h>
#include <preprocessing/PreProcessorKineticsPolicy_CHEMKIN.h>
#include <preprocessing/PreProcessorSpeciesPolicy_CHEMKIN_WithTransport.h>

// Maps
#include <maps/ThermodynamicsMap_CHEMKIN.h>
#include <maps/TransportPropertiesMap_CHEMKIN.h>
#include <maps/KineticsMap_CHEMKIN.h>

// Typedefs
typedef OpenSMOKE::Species<OpenSMOKE::ThermoPolicy_CHEMKIN,
                           OpenSMOKE::TransportPolicy_CHEMKIN>
    SpeciesCHEMKIN;

typedef OpenSMOKE::PreProcessorSpecies<
    OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<SpeciesCHEMKIN> >
    PreProcessorSpecies_CHEMKIN_WithoutTransport;

typedef OpenSMOKE::PreProcessorKinetics<
    OpenSMOKE::PreProcessorKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_CHEMKIN> >
    PreProcessorKinetics_CHEMKIN;

typedef OpenSMOKE::ThermoReader<
    OpenSMOKE::ThermoReaderPolicy_CHEMKIN<OpenSMOKE::ThermoPolicy_CHEMKIN> >
    ThermoReader_CHEMKIN;

// Boost Library
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
namespace fs = boost::filesystem;
namespace po = boost::program_options;

const double UNFEASIBLE_BIG_NUMBER = 1.e16;

// Dakota Library
#include <ParallelLibrary.hpp>
#include <ProblemDescDB.hpp>
#include <LibraryEnvironment.hpp>
#include <DakotaModel.hpp>
#include <DakotaInterface.hpp>
#include <DakotaResponse.hpp>
#include <ParamResponsePair.hpp>
#include <DirectApplicInterface.hpp>

#include <nlopt.hpp>

// Curve Matching
#include "curve_matching/curve_matching.h"

// Header files
#include "utilities/OptiSMOKEUtilities"
#include "grammar/grammar.h"
#include "options/options.h"
#include "ideal_reactors/ideal_reactors.h"
#include "1d_flames/1d_flames.h"
#include "DataManager.h"
#include "InputManager.h"
#include "OptimizedKinetics.h"
#include "SerialDakotaInterface.h"
#include "SimulationsInterface.h"

double NLOptFunction(const vector<double> &x, vector<double> &grad, void *my_func_data);
double OptFunction(const vector<double> &b, unsigned int fn_val);

#ifdef HAVE_AMPL
// Floating-point initialization from AMPL: switch to 53-bit rounding
// if appropriate, to eliminate some cross-platform differences.
extern "C" void fpinit_ASL();
#endif

#ifndef DAKOTA_HAVE_MPI
#define MPI_COMM_WORLD 0
#endif  // not DAKOTA_HAVE_MPI

// Run a Dakota LibraryEnvironment, mode 1: parsing an input file
void run_dakota_parse(const char *plugin_input_file, bool echo_dakota_string);

void opensmoke_interface_plugin(
    Dakota::LibraryEnvironment &env);  //,const char* plugin_input_file);

OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
OptiSMOKE::InputManager input(dictionaries);

// #if OPTISMOKE_USE_NLOPT
OptiSMOKE::SimulationsInterface *sim_iface_;
OptiSMOKE::OptimizedKinetics *opti_kinetics_;
unsigned int numberOfGradientEvaluations;
unsigned int numberOfFunctionEvaluations;
bool violated_uncertainty;
double prev_fn_val;
std::ofstream fOut;
// # endif
