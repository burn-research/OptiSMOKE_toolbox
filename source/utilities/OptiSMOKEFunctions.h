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

#ifndef OPTISMOKE_OPTISMOKEFUNCTIONS_H
#define OPTISMOKE_OPTISMOKEFUNCTIONS_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <stdarg.h>
#include <assert.h>
#include <vector>
#include <map>
#include <stdlib.h>

#define __OPTISMOKE_VERSION__ "3.0.0-beta"
#define OPTISMOKE_FATAL_ERROR_EXIT -1
#define OPTISMOKE_SUCCESSFULL_EXIT  0

namespace OptiSMOKE{
    
    /**
	* Utility to print fatal error messages
	*/
	void ErrorMessage(const std::string functionName, const std::string errorMessage);

	/**
	* Utility to print fatal error messages
	*/
	int FatalErrorMessage(const std::string errorMessage);

    /**
	*@brief Writes the OptiSMOKE++ logo on the screen
	*@param application_name name of the OptiSMOKE++ solver
	*@param author_name author's name
	*/
	void OptiSMOKE_logo(const std::string application_name, const std::string author_name);
}

#include "OptiSMOKEFunctions.hpp"
#endif //OPTISMOKE_OPTISMOKEFUNCTIONS_H
