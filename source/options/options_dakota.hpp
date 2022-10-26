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
		string_population_size = "50";
		fitness_type = "merit_function";
		mutation_type = "offset_normal";
		mutation_rate = "1.0";
		crossover_type = "two_point";
		crossover_rate = "0.0";
		replacement_type = "chc = 10";

		// Values for coliny_direct
		string_division = "major_dimension";
		max_boxsize_limit = "0.0";
		min_boxsize_limit = "1.0e-4";

        gradient_option = false;

        tabular_data_file = "tabulara_data.dat";
    }

    void options_dakota::SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager, 
    std::string dictionary_name){
        dictionary_manager(dictionary_name).SetGrammar(dakota_options_grammar_);

        if (dictionaries(dakota_options_dictionary).CheckOption("@TabularDataFile") == true)
			dictionaries(dakota_options_dictionary).ReadString("@TabularDataFile", tabular_data_file);

		if (dictionaries(dakota_options_dictionary).CheckOption("@Method") == true)
			dictionaries(dakota_options_dictionary).ReadString("@Method", method);

		if (dictionaries(dakota_options_dictionary).CheckOption("@MaxIterations") == true)
			dictionaries(dakota_options_dictionary).ReadString("@MaxIterations", string_max_iterations);

		if (dictionaries(dakota_options_dictionary).CheckOption("@MaxFunctionEvaluations") == true)
			dictionaries(dakota_options_dictionary).ReadString("@MaxFunctionEvaluations", string_max_function_evaluations);

		if (dictionaries(dakota_options_dictionary).CheckOption("@ConvergenceTolerance") == true)
			dictionaries(dakota_options_dictionary).ReadString("@ConvergenceTolerance", string_convergence_tolerance);

		if (dictionaries(dakota_options_dictionary).CheckOption("@SolutionTarget") == true)
			dictionaries(dakota_options_dictionary).ReadString("@SolutionTarget", string_solution_target);

		if (dictionaries(dakota_options_dictionary).CheckOption("@Seed") == true)
			dictionaries(dakota_options_dictionary).ReadString("@Seed", string_seed);

		if (dictionaries(dakota_options_dictionary).CheckOption("@PopulationSize") == true)
			dictionaries(dakota_options_dictionary).ReadString("@PopulationSize", string_population_size);

		if (dictionaries(dakota_options_dictionary).CheckOption("@FitnessType") == true)
			dictionaries(dakota_options_dictionary).ReadString("@FitnessType", fitness_type);

		if (dictionaries(dakota_options_dictionary).CheckOption("@MutationType") == true)
			dictionaries(dakota_options_dictionary).ReadString("@MutationType", mutation_type);

		if (dictionaries(dakota_options_dictionary).CheckOption("@MutationRate") == true)
			dictionaries(dakota_options_dictionary).ReadString("@MutationRate", mutation_rate);

		if (dictionaries(dakota_options_dictionary).CheckOption("@CrossoverType") == true)
			dictionaries(dakota_options_dictionary).ReadString("@CrossoverType", crossover_type);

		if (dictionaries(dakota_options_dictionary).CheckOption("@CrossoverRate") == true)
			dictionaries(dakota_options_dictionary).ReadString("@CrossoverRate", crossover_rate);

		if (dictionaries(dakota_options_dictionary).CheckOption("@ReplacementType") == true)
			dictionaries(dakota_options_dictionary).ReadString("@ReplacementType", replacement_type);

		if (dictionaries(dakota_options_dictionary).CheckOption("@Division") == true)
			dictionaries(dakota_options_dictionary).ReadString("@Division", string_division);

		if (dictionaries(dakota_options_dictionary).CheckOption("@MaxBoxsizeLimit") == true)
			dictionaries(dakota_options_dictionary).ReadString("@MaxBoxsizeLimit", max_boxsize_limit);

		if (dictionaries(dakota_options_dictionary).CheckOption("@MinBoxsizeLimit") == true)
			dictionaries(dakota_options_dictionary).ReadString("@MinBoxsizeLimit", min_boxsize_limit);

		if (dictionaries(dakota_options_dictionary).CheckOption("@DiverseInput") == true)
			dictionaries(dakota_options_dictionary).ReadOption("@DiverseInput", diverse_dakota_input);
			
 		if (dictionaries(dakota_options_dictionary).CheckOption("@Gradient") == true)
            dictionaries(dakota_options_dictionary).ReadBool("@Gradient", gradient_option);
	}

    void options_dakota::WriteDakotaInputString(){

        if(!tabular_data_file.empty()){
		    dakota_options_string = "     environment,"
		  	                        "\n      tabular_data";
		 	dakota_options_string.append("\n 		tabular_data_file '" + tabular_data_file + "'");
		}
		dakota_options_string.append("\n	method,"); 
		dakota_options_string.append("\n 		" + method);
		if(!string_max_iterations.empty()){
		  	dakota_options_string.append("\n		  max_iterations = " + string_max_iterations);
		}
		if(!string_max_function_evaluations.empty()){
		    dakota_options_string.append("\n 		  max_function_evaluations = " + string_max_function_evaluations);
		}
		if(!string_convergence_tolerance.empty()){
		    dakota_options_string.append("\n		  convergence_tolerance = " + string_convergence_tolerance);
		}
		if(!string_solution_target.empty()){
	        dakota_options_string.append("\n 		  solution_target = " + string_solution_target);
		}
		if(!string_seed.empty())
		{
		    dakota_options_string.append("\n 		  seed = " + string_seed);
		}
		if(!diverse_dakota_input.empty()){
		    dakota_options_string.append("\n");
			for (int i = 0; i < diverse_dakota_input.size(); i++){
				dakota_options_string.append( " " + diverse_dakota_input[i]);
			}
		} 
        else if(method == "coliny_ea"){
		    dakota_options_string.append( "\n 		  population_size = " + string_population_size);
		    dakota_options_string.append( "\n		  fitness_type " + fitness_type);
		    dakota_options_string.append( "\n		  mutation_type " + mutation_type);
		    dakota_options_string.append( "\n		  mutation_rate " + mutation_rate);
		    dakota_options_string.append( "\n		  crossover_type " + crossover_type);
		    dakota_options_string.append( "\n		  crossover_rate " + crossover_rate);
		    dakota_options_string.append( "\n		  replacement_type " + replacement_type);
		}
        else if(method == "coliny_direct"){
		    dakota_options_string.append( "\n                 division " + string_division);
		    dakota_options_string.append( "\n                 max_boxsize_limit " + max_boxsize_limit);
		    dakota_options_string.append( "\n                 min_boxsize_limit " + min_boxsize_limit);
		}
		dakota_options_string.append( "\n	variables,");
		if( distribution == "uniform"){
			dakota_options_string.append( "\n		continuous_design = " + std::to_string(number_of_parameters));
			dakota_options_string.append( "\n		  descriptors " + param_name_string);
			dakota_options_string.append( "\n 		  initial_point " + initial_values_string);
			dakota_options_string.append( "\n 	 	  lower_bounds " + lower_bounds_string);
			dakota_options_string.append( "\n 		  upper_bounds " + upper_bounds_string);
		} 
        else if (distribution == "normal"){
			dakota_options_string.append( "\n               active uncertain " );
			dakota_options_string.append( "\n		normal_uncertain = " + std::to_string(number_of_parameters));
			dakota_options_string.append( "\n		  descriptors " + param_name_string);
			dakota_options_string.append( "\n 		  means " + initial_values_string);
			dakota_options_string.append( "\n 		  std_deviations " + std_deviations_string);
		}

		dakota_options_string.append( "\n	interface,");
		dakota_options_string.append( "\n		direct");
		dakota_options_string.append( "\n		  analysis_driver = 'plugin_opensmoke'");
		dakota_options_string.append( "\n	responses,");
		dakota_options_string.append( "\n		num_objective_functions = 1");
			
		// Options to use other optimization method (gradient-based)

		if (gradient_option == true){
			dakota_options_string.append( "\n               numerical_gradients");
			dakota_options_string.append( "\n               	method_source dakota");
			dakota_options_string.append( "\n               	interval_type forward");
			dakota_options_string.append( "\n               	fd_step_size = 1.e-5");
		}
		else{
			dakota_options_string.append( "\n		no_gradients");
		}

		dakota_options_string.append( "\n		no_hessians");

    }
}