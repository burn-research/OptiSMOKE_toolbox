/* ------------------------------------------------------------------------------- *\
|                                                                                   |
|             ____        __  _ _____ __  _______  __ __ ______                     |
|            / __ \____  / /_(_) ___//  |/  / __ \/ //_// ____/___  ____            |
|           / / / / __ \/ __/ /\__ \/ /|_/ / / / / ,<  / __/ / __ \/ __ \           |
|          / /_/ / /_/ / /_/ /___/ / /  / / /_/ / /| |/ /___/ /_/ / /_/ /           |
|          \____/ .___/\__/_//____/_/  /_/\____/_/ |_/_____/ .___/ .___/            |
|              /_/                                        /_/   /_/                 |
|                                                                                   |
| --------------------------------------------------------------------------------- |
|  Please refer to the copyright statement and license                              |
|  information at the end of this file.                                             |
| --------------------------------------------------------------------------------- |
|                                                                                   |
|           Author: Timoteo Dinelli <timoteo.dinelli@polimi.it>                     |
|              CRECK Modeling Group <http://creckmodeling.chem.polimi.it>           |
|              Department of Chemistry, Materials and Chemical Engineering          |
|              Politecnico di Milano                                                |
|              P.zza Leonardo da Vinci 32, 20133 Milano                             |
|                                                                                   |
\* ------------------------------------------------------------------------------- */
#pragma once

#include <string>
#include <iostream>
#include <numeric>
#include <random>
#include <chrono>
#include <algorithm>

// ==================================================
// Boost Headers
// ==================================================
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
namespace fs = boost::filesystem;
namespace po = boost::program_options;

// ==================================================
// Dakota Headers
// ==================================================
#include <ParallelLibrary.hpp>
#include <ProblemDescDB.hpp>
#include <LibraryEnvironment.hpp>
#include <DakotaModel.hpp>
#include <DakotaInterface.hpp>
#include <DakotaResponse.hpp>
#include <ParamResponsePair.hpp>
#include <DirectApplicInterface.hpp>
#ifdef HAVE_AMPL
// Floating-point initialization from AMPL: switch to 53-bit rounding if appropriate, to eliminate some cross-platform
// differences.
extern "C" void fpinit_ASL();
#endif
#ifndef DAKOTA_HAVE_MPI
#define MPI_COMM_WORLD 0
#endif  // not DAKOTA_HAVE_MPI

// #include <OpenSMOKEpp>
//
// // Thermodynamics
// #include <kernel/thermo/Species.h>
// #include <kernel/thermo/ThermoPolicy_CHEMKIN.h>
// #include <kernel/thermo/ThermoReader.h>
// #include <kernel/thermo/ThermoReaderPolicy_CHEMKIN.h>
//
// // Transport
// #include <kernel/transport/TransportPolicy_CHEMKIN.h>
// #include <kernel/transport/TransportReader.h>
// #include <kernel/transport/TransportReaderPolicy_CHEMKIN.h>
//
// // Kinetics
// #include <kernel/kinetics/ReactionPolicy_CHEMKIN.h>
//
// // Preprocessing
// #include <preprocessing/PreProcessorSpecies.h>
// #include <preprocessing/PreProcessorKinetics.h>
// #include <preprocessing/PreProcessorKineticsPolicy_CHEMKIN.h>
// #include <preprocessing/PreProcessorSpeciesPolicy_CHEMKIN_WithTransport.h>
//
// Maps
// #include <maps/ThermodynamicsMap_CHEMKIN.h>
// #include <maps/TransportPropertiesMap_CHEMKIN.h>
// #include <maps/KineticsMap_CHEMKIN.h>
//
// Typedefs
// typedef OpenSMOKE::Species<OpenSMOKE::ThermoPolicy_CHEMKIN, OpenSMOKE::TransportPolicy_CHEMKIN> SpeciesCHEMKIN;
//
// typedef OpenSMOKE::PreProcessorSpecies<OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<SpeciesCHEMKIN>>
//     PreProcessorSpecies_CHEMKIN_WithoutTransport;
//
// typedef OpenSMOKE::PreProcessorKinetics<
//     OpenSMOKE::PreProcessorKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_CHEMKIN>>
//     PreProcessorKinetics_CHEMKIN;
//
// typedef OpenSMOKE::ThermoReader<OpenSMOKE::ThermoReaderPolicy_CHEMKIN<OpenSMOKE::ThermoPolicy_CHEMKIN>>
//     ThermoReader_CHEMKIN;
//
//
// const double UNFEASIBLE_BIG_NUMBER = 1.e16;
//
//
// // NLopt++
// // #include <nlopt.hpp>
//
// // Curve Matching
// #include "curve_matching/curve_matching.h"
//
// // Header files
// #include "utilities/OptiSMOKEUtilities"
// #include "grammar/grammar.h"
// #include "options/options.h"
// #include "ideal_reactors/ideal_reactors.h"
// #include "1d_flames/1d_flames.h"
// #include "DataManager.h"
// #include "InputManager.h"
// #include "OptimizedKinetics.h"
// #include "SerialDakotaInterface.h"
// #include "SimulationsInterface.h"
//
// double NLOptFunction(const vector<double>& x, vector<double>& grad, void* my_func_data);
// double OptFunction(const vector<double>& b, unsigned int fn_val);
//
//
// // Run a Dakota LibraryEnvironment, mode 1: parsing an input file
// void run_dakota_parse(const char* plugin_input_file, bool echo_dakota_string);
//
// void opensmoke_interface_plugin(Dakota::LibraryEnvironment& env);  //,const char* plugin_input_file);
//
// OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
// OptiSMOKE::InputManager input(dictionaries);
//
// // #if OPTISMOKE_USE_NLOPT
// OptiSMOKE::SimulationsInterface* sim_iface_;
// OptiSMOKE::OptimizedKinetics* opti_kinetics_;
// unsigned int numberOfGradientEvaluations;
// unsigned int numberOfFunctionEvaluations;
// bool violated_uncertainty;
// double prev_fn_val;
// std::ofstream fOut;
// // # endif
/* ------------------------------------------------------------------------------- *\
|                                                                                   |
|   MIT License                                                                     |
|                                                                                   |
|   Copyright (c) 2025 Timoteo Dinelli                                              |
|                                                                                   |
|   Permission is hereby granted, free of charge, to any person obtaining a copy    |
|   of this software and associated documentation files (the "Software"), to deal   |
|   in the Software without restriction, including without limitation the rights    |
|   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       |
|   copies of the Software, and to permit persons to whom the Software is           |
|   furnished to do so, subject to the following conditions:                        |
|                                                                                   |
|   The above copyright notice and this permission notice shall be included in all  |
|   copies or substantial portions of the Software.                                 |
|                                                                                   |
|   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      |
|   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        |
|   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     |
|   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          |
|   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   |
|   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   |
|   SOFTWARE.                                                                       |
|                                                                                   |
\* ------------------------------------------------------------------------------- */
