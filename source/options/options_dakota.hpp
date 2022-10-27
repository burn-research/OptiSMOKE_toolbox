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

namespace OptiSMOKE
{
    options_dakota::options_dakota(){

        method = "coliny_ea";
		population_size = "50";
		fitness_type = "merit_function";
		mutation_type = "offset_normal";
		mutation_rate = "1.0";
		crossover_type = "two_point";
		crossover_rate = "0.0";
		replacement_type = "chc = 10";

		// Values for coliny_direct
		division = "major_dimension";
		max_boxsize_limit = "0.0";
		min_boxsize_limit = "1.0e-4";

        gradient_option = false;

        tabular_data_file = "tabulara_data.dat";
    }

    void options_dakota::SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager, 
											std::string dictionary_name)
	{

        dictionary_manager(dictionary_name).SetGrammar(dakota_options_grammar_);

        if (dictionary_manager(dictionary_name).CheckOption("@TabularDataFile") == true)
			dictionary_manager(dictionary_name).ReadString("@TabularDataFile", tabular_data_file);

		if (dictionary_manager(dictionary_name).CheckOption("@Method") == true)
			dictionary_manager(dictionary_name).ReadString("@Method", method);

		if (dictionary_manager(dictionary_name).CheckOption("@MaxIterations") == true)
			dictionary_manager(dictionary_name).ReadString("@MaxIterations", max_iterations);

		if (dictionary_manager(dictionary_name).CheckOption("@MaxFunctionEvaluations") == true)
			dictionary_manager(dictionary_name).ReadString("@MaxFunctionEvaluations", max_function_evaluations);

		if (dictionary_manager(dictionary_name).CheckOption("@ConvergenceTolerance") == true)
			dictionary_manager(dictionary_name).ReadString("@ConvergenceTolerance", convergence_tolerance);

		if (dictionary_manager(dictionary_name).CheckOption("@SolutionTarget") == true)
			dictionary_manager(dictionary_name).ReadString("@SolutionTarget", solution_target);

		if (dictionary_manager(dictionary_name).CheckOption("@Seed") == true)
			dictionary_manager(dictionary_name).ReadString("@Seed", seed);

		if (dictionary_manager(dictionary_name).CheckOption("@PopulationSize") == true)
			dictionary_manager(dictionary_name).ReadString("@PopulationSize", population_size);

		if (dictionary_manager(dictionary_name).CheckOption("@FitnessType") == true)
			dictionary_manager(dictionary_name).ReadString("@FitnessType", fitness_type);

		if (dictionary_manager(dictionary_name).CheckOption("@MutationType") == true)
			dictionary_manager(dictionary_name).ReadString("@MutationType", mutation_type);

		if (dictionary_manager(dictionary_name).CheckOption("@MutationRate") == true)
			dictionary_manager(dictionary_name).ReadString("@MutationRate", mutation_rate);

		if (dictionary_manager(dictionary_name).CheckOption("@CrossoverType") == true)
			dictionary_manager(dictionary_name).ReadString("@CrossoverType", crossover_type);

		if (dictionary_manager(dictionary_name).CheckOption("@CrossoverRate") == true)
			dictionary_manager(dictionary_name).ReadString("@CrossoverRate", crossover_rate);

		if (dictionary_manager(dictionary_name).CheckOption("@ReplacementType") == true)
			dictionary_manager(dictionary_name).ReadString("@ReplacementType", replacement_type);

		if (dictionary_manager(dictionary_name).CheckOption("@Division") == true)
			dictionary_manager(dictionary_name).ReadString("@Division", division);

		if (dictionary_manager(dictionary_name).CheckOption("@MaxBoxsizeLimit") == true)
			dictionary_manager(dictionary_name).ReadString("@MaxBoxsizeLimit", max_boxsize_limit);

		if (dictionary_manager(dictionary_name).CheckOption("@MinBoxsizeLimit") == true)
			dictionary_manager(dictionary_name).ReadString("@MinBoxsizeLimit", min_boxsize_limit);

		if (dictionary_manager(dictionary_name).CheckOption("@DiverseInput") == true)
			dictionary_manager(dictionary_name).ReadOption("@DiverseInput", diverse_dakota_input);
			
 		if (dictionary_manager(dictionary_name).CheckOption("@Gradient") == true)
            dictionary_manager(dictionary_name).ReadBool("@Gradient", gradient_option);
	}
}