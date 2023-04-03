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
|					  Timoteo Dinelli <timoteo.dinelli@polimi.it>	      |
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

#ifndef OPTIONS_NLOPT_H
#define OPTIONS_NLOPT_H

namespace OptiSMOKE
{
	class options_nlopt
	{
	public:
		options_nlopt();
		
		~options_nlopt(){};

		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
								std::string dictionary_name);
    
		inline const int& algo_int() const {return algo_int_;};
		inline const int& max_function_evaluations() const {return max_function_evaluations_;};
    private:

		grammar_nlopt nlopt_options_grammar_;

		std::string algorithm_;
		std::string variant_;
		int max_iterations_;
		int max_function_evaluations_;
		double convergence_tolerance_;
		double solution_target_;

		bool global_;
		bool local_;

		int algo_int_;

		void SetupAlgorithm();
	};
}

#include "options_nlopt.hpp"
#endif // OPTIONS_NLOPT_H