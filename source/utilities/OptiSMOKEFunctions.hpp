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

namespace OptiSMOKE{
    
    void ErrorMessage(const std::string functionName, const std::string errorMessage)
	{
		std::cout << "Function:     " << functionName << std::endl;
		std::cout << "Fatal error:  " << errorMessage << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(OPTISMOKE_FATAL_ERROR_EXIT);
	}

	int FatalErrorMessage(const std::string errorMessage)
	{
		std::cout << "Fatal error:  " << errorMessage << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(OPTISMOKE_FATAL_ERROR_EXIT);
		return OPTISMOKE_FATAL_ERROR_EXIT;
	}

	void OptiSMOKE_logo(const std::string application_name, const std::string author_name)
	{
		std::string current_time = __TIME__;
		std::string current_date = __DATE__;
		std::string author_complete = "Authors: " + author_name;
		std::string compilation_time = "Compilation date: " + current_date + " at " + current_time;
		std::string version = "Version: "; version += __OPTISMOKE_VERSION__;

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << std::endl;
    	std::cout << "           ____            _  ______ __  __  ____  _  ________              " << std::endl;
    	std::cout << "          / __ \\       _  (_)/  ___ |  \\/  |/ __ \\| |/ /  ____|             " << std::endl;
    	std::cout << "         | |  | |_ __ | |_ _ | (___ | \\  / | |  | | ' /| |__    _     _     " << std::endl;
    	std::cout << "         | |  | | '_ \\|  _| |\\___  \\| |\\/| | |  | |  < |  __| _| |_ _| |_   " << std::endl;
    	std::cout << "         | |__| | |_) | |_| |____)  | |  | | |__| | . \\| |___|_   _|_   _|  " << std::endl;
    	std::cout << "          \\____/| .__/\\___|_|______/|_|  |_|\\____/|_|\\_\\______||_|   |_|    " << std::endl;
    	std::cout << "                | |                                                         " << std::endl;
    	std::cout << "                |_|                                                         " << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "           Department of Chemistry, Materials and Chemical Engineering        " << std::endl;
		std::cout << "                              Politecnico di Milano                           " << std::endl;
		std::cout << "                         http://www.opensmoke.polimi.it/                      " << std::endl;
		std::cout << "                      http://creckmodeling.chem.polimi.it/                    " << std::endl;
		std::cout << std::endl;
		for (unsigned int i = 1; i <= 39 - (unsigned int)(application_name.size()) / 2; i++)	std::cout << " ";
		std::cout << application_name << std::endl;
		for (unsigned int i = 1; i <= 39 - (unsigned int)(version.size()) / 2; i++)	std::cout << " ";
		std::cout << version << std::endl;
		for (unsigned int i = 1; i <= 39 - (unsigned int)(author_complete.size()) / 2; i++)	std::cout << " ";
		std::cout << author_complete << std::endl;
		for (unsigned int i = 1; i <= 39 - (unsigned int)(compilation_time.size()) / 2; i++)	std::cout << " ";
		std::cout << compilation_time << std::endl;
		std::cout << std::endl;

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "                                  WARNING                                    " << std::endl;
		std::cout << "   This version of OpenSMOKE++ Suite can be used for educational purposes    " << std::endl;
		std::cout << "              only and cannot be distributed to third parties.               " << std::endl;
		std::cout << "       The software is and remains the sole property of Alberto Cuoci.       " << std::endl;
		std::cout << "      Whenever the OpenSMOKE++ Suite is used to produce any publication,     " << std::endl;
		std::cout << "       a detailed reference to the OpenSMOKE++ code should be reported       " << std::endl;
		std::cout << "                            (see User's Guide).                              " << std::endl;
		std::cout << "    Use for commercial purposes is not permitted. For any commercial issue   " << std::endl;
		std::cout << "         please contact Alberto Cuoci (email: alberto.cuoci@polimi.it)       " << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "                            LIMITED WARRANTY                                 " << std::endl;
		std::cout << "     This software is provided \"as is\" and without warranties as to        " << std::endl;
		std::cout << "  performance of merchantability or any other warranties whether expressed   " << std::endl;
		std::cout << "    or implied. Because of the various hardware and software environments    " << std::endl;
		std::cout << "   into which this library may be installed, no warranty of fitness for a    " << std::endl;
		std::cout << "   particular purpose is offered. The user must assume the entire risk of    " << std::endl;
		std::cout << "                          using  the library.                                " << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
	}
}