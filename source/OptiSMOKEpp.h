/* ----------------------------------------------------------------------------------- *\
|                                                                                       |
|               ____        __  _ _____ __  _______  __ __ ______                       |
|              / __ \____  / /_(_) ___//  |/  / __ \/ //_// ____/___  ____              |
|             / / / / __ \/ __/ /\__ \/ /|_/ / / / / ,<  / __/ / __ \/ __ \             |
|            / /_/ / /_/ / /_/ /___/ / /  / / /_/ / /| |/ /___/ /_/ / /_/ /             |
|            \____/ .___/\__/_//____/_/  /_/\____/_/ |_/_____/ .___/ .___/              |
|                /_/                                        /_/   /_/                   |
|                                                                                       |
| ------------------------------------------------------------------------------------- |
|  See license and copyright at the end of this file.                                   |
| ------------------------------------------------------------------------------------- |
|                                                                                       |
|          Authors: Timoteo Dinelli  <timoteo.dinelli@polimi.it>                        |
|                   Andrea Bertolino <andrea.bertolino@ulb.be>                          |
|                   Magnus Fürst     <magnus.furst@ulb.ac.be>                           |
|                                                                                       |
|          [1] CRECK Modeling Lab <https://www.creckmodeling.polimi.it>                 |
|              Department of Chemistry, Materials and Chemical Engineering              |
|              Politecnico di Milano, P.zza Leonardo da Vinci 32, 20133 Milano          |
|                                                                                       |
|          [2] BRITE Research Group <https://brite-research.be>                         |
|              Brussels Institute for Thermal-fluid systems and clean Energy            |
|              Avenue F.D. Rooseveltlaan 50, Bruxelles 1050 Brussel                     |
|                                                                                       |
\* ----------------------------------------------------------------------------------- */
#pragma once

#include <string>
#include <iostream>
#include <numeric>
#include <random>
#include <chrono>
#include <algorithm>
#include <memory>

// ==================================================
// Boost Headers and definitions
// ==================================================
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/json.hpp>
#include <boost/optional.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;
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

// ==================================================
// NLOpt++ Headers
// ==================================================
// #include <nlopt.hpp>

// ==================================================
// OpenSMOKEpp Headers
// ==================================================
// Thermodynamics
#include <kernel/thermo/Species.h>
#include <kernel/thermo/ThermoPolicy_CHEMKIN.h>
#include <kernel/thermo/ThermoReader.h>
#include <kernel/thermo/ThermoReaderPolicy_CHEMKIN.h>
// ==================================================
// Transport
#include <kernel/transport/TransportPolicy_CHEMKIN.h>
#include <kernel/transport/TransportReader.h>
#include <kernel/transport/TransportReaderPolicy_CHEMKIN.h>
// ==================================================
// Kinetics
#include <kernel/kinetics/ReactionPolicy_CHEMKIN.h>
// ==================================================
// Preprocessing
#include <preprocessing/PreProcessorSpecies.h>
#include <preprocessing/PreProcessorKinetics.h>
#include <preprocessing/PreProcessorKineticsPolicy_CHEMKIN.h>
#include <preprocessing/PreProcessorSpeciesPolicy_CHEMKIN_WithTransport.h>
// ==================================================
// Maps
#include <maps/ThermodynamicsMap_CHEMKIN.h>
#include <maps/TransportPropertiesMap_CHEMKIN.h>
#include <maps/KineticsMap_CHEMKIN.h>
// ==================================================
// OpenSMOKE Dictionaries stuff
#include <dictionary/OpenSMOKE_DictionaryManager.h>
#include <dictionary/OpenSMOKE_DictionaryGrammar.h>
#include <dictionary/OpenSMOKE_DictionaryKeyWord.h>
// ==================================================
// Typedefs
typedef OpenSMOKE::Species<OpenSMOKE::ThermoPolicy_CHEMKIN, OpenSMOKE::TransportPolicy_CHEMKIN> SpeciesCHEMKIN;
typedef OpenSMOKE::PreProcessorSpecies<OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<SpeciesCHEMKIN>>
    PreProcessorSpecies_CHEMKIN_WithoutTransport;
typedef OpenSMOKE::PreProcessorKinetics<
    OpenSMOKE::PreProcessorKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_CHEMKIN>>
    PreProcessorKinetics_CHEMKIN;
typedef OpenSMOKE::ThermoReader<OpenSMOKE::ThermoReaderPolicy_CHEMKIN<OpenSMOKE::ThermoPolicy_CHEMKIN>>
    ThermoReader_CHEMKIN;

// ==================================================
// Internal Headers
// ==================================================
#include "DataStructures.h"
#include "utilities/OptiSMOKEUtilities"
#include "grammar/Grammar.h"
#include "options/Options.h"
// #include "ideal_reactors/ideal_reactors.h"
// #include "1d_flames/1d_flames.h"
#include "DataManager.h"
#include "InputManager.h"
// #include "OptimizedKinetics.h"
// #include "SerialDakotaInterface.h"
// #include "SimulationsInterface.h"
// #include "curve_matching/curve_matching.h"
/* ------------------------------------------------------------------------------- *\
|                                                                                   |
|   MIT License                                                                     |
|                                                                                   |
|   Copyright (c) 2025 Timoteo Dinelli, Andrea Bertolino, Magnus Fürst              |
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
