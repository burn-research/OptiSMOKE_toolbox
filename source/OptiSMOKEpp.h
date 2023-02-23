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

// Standard Library
#include <string>
#include <iostream>
#include <numeric>

// OpenSMOKE header
#include "OpenSMOKE_Definitions.h"
#include "dictionary/OpenSMOKE_DictionaryManager.h"

// Base classes
#include "kernel/thermo/Thermodynamics_CHEMKIN"
#include "kernel/transport/Transport_CHEMKIN"
#include "kernel/kinetics/ReactionPolicy_CHEMKIN.h"

// CHEMKIN maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/TransportPropertiesMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"

// Boost Library
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

// Eigen Library
#include <Eigen/Dense>

// Dakota Library
#include "ParallelLibrary.hpp"
#include "ProblemDescDB.hpp"
#include "LibraryEnvironment.hpp"
#include "DakotaModel.hpp"
#include "DakotaInterface.hpp"
#include "PluginSerialDirectApplicInterface.hpp"

// Header files
#include "utilities/OptiSMOKEUtilities"
#include "grammar/grammar.h"
#include "options/options.h"
#include "InputManager.h"

#ifdef HAVE_AMPL 
// Floating-point initialization from AMPL: switch to 53-bit rounding
// if appropriate, to eliminate some cross-platform differences.
extern "C" void fpinit_ASL(); 
#endif 

#ifndef DAKOTA_HAVE_MPI
#define MPI_COMM_WORLD 0
#endif // not DAKOTA_HAVE_MPI


/// Run a Dakota LibraryEnvironment, mode 1: parsing an input file
void run_dakota_parse(const char* plugin_input_file);

void opensmoke_interface_plugin(Dakota::LibraryEnvironment& env); //,const char* plugin_input_file);
