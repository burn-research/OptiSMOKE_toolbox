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

#include "grammar/grammar_DakotaOptions.h"

namespace OptiSMOKE
{
  class options_dakota
  {
    public:
    options_dakota();
    ~options_dakota() {};

    void SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                            std::string dictionary_name);
    
    void WriteDakotaInputString();

    std::string dakota_input_string;
    
    private:

    grammar_DakotaOptions dakota_options_grammar_;

    std::string method;
		std::string string_population_size;
		std::string fitness_type;
		std::string mutation_type;
		std::string mutation_rate;
		std::string crossover_type;
		std::string crossover_rate;
		std::string replacement_type;
		
    // Values for coliny_direct
		std::string string_division;
		std::string max_boxsize_limit;
		std::string min_boxsize_limit;

    bool gradient_option;

    std::string tabular_data_file;
  };
}

#include "options_dakota.hpp"
#endif // OPTIONS_DAKOTA_H