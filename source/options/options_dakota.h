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
|					  Timoteo Dinelli <timoteo.dinelli@polimi.it>	                  |
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

#ifndef OPTIONS_DAKOTA_H
#define OPTIONS_DAKOTA_H

namespace OptiSMOKE
{
	class options_dakota
	{
	public:
		options_dakota();
		
		~options_dakota() {};

		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
								std::string dictionary_name);
    
		// Public access elements
    	inline const std::string& method() const {return method_;};
    	inline const std::string& population_size() const {return population_size_;};
    	inline const std::string& fitness_type() const {return fitness_type_;};
    	inline const std::string& mutation_type() const {return mutation_type_;};
    	inline const std::string& mutation_rate() const {return mutation_rate_;};
   		inline const std::string& crossover_type() const {return crossover_type_;};
    	inline const std::string& crossover_rate() const {return crossover_rate_;};
    	inline const std::string& replacement_type() const {return replacement_type_;};
    	inline const std::string& max_iterations() const {return max_iterations_;};
    	inline const std::string& max_function_evaluations() const {return max_function_evaluations_;};
    	inline const std::string& convergence_tolerance() const {return convergence_tolerance_;};
    	inline const std::string& solution_target() const {return solution_target_;};
    	inline const std::string& seed() const {return seed_;};
    	inline const std::vector<std::string>& diverse_dakota_input() const {return diverse_dakota_input_;};
    	inline const std::string& division() const {return division_;};
    	inline const std::string& max_boxsize_limit() const {return max_boxsize_limit_;};
    	inline const std::string& min_boxsize_limit() const {return min_boxsize_limit_;};
	    inline const bool& dakota_gradient() const {return dakota_gradient_;};
    	inline const bool& diverse_input() const {return diverse_input_;};
    	inline const std::string& tabular_data_file() const {return tabular_data_file_;};
    	inline const bool& echo_dakota_string() const {return echo_dakota_string_;};

    private:

    	grammar_dakota dakota_options_grammar_;

    	std::string method_;
		std::string population_size_;
		std::string fitness_type_;
		std::string mutation_type_;
		std::string mutation_rate_;
		std::string crossover_type_;
		std::string crossover_rate_;
		std::string replacement_type_;
		std::string max_iterations_;
    	std::string max_function_evaluations_;
    	std::string convergence_tolerance_;
    	std::string solution_target_;
    	std::string seed_;
    	std::vector<std::string> diverse_dakota_input_;
    	bool echo_dakota_string_;

 	   // Values for coliny_direct
		std::string division_;
		std::string max_boxsize_limit_;
		std::string min_boxsize_limit_;

    	bool dakota_gradient_;
    	bool diverse_input_;

	    std::string tabular_data_file_;
	};
}

#include "options_dakota.hpp"
#endif // OPTIONS_DAKOTA_H