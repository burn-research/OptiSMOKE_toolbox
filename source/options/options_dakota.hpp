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

        method_ = "coliny_ea";
		population_size_ = "50";
		fitness_type_ = "merit_function";
		mutation_type_ = "offset_normal";
		mutation_rate_ = "1.0";
		crossover_type_ = "two_point";
		crossover_rate_ = "0.0";
		replacement_type_ = "chc = 10";

		// Values for coliny_direct
		division_ = "major_dimension";
		max_boxsize_limit_ = "0.0";
		min_boxsize_limit_ = "1.0e-4";

        dakota_gradient_ = false;
		diverse_input_ = false;

        tabular_data_file_ = "tabular_data.dat";

		echo_dakota_string_ = false;
    }

    void options_dakota::SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager, 
											std::string dictionary_name)
	{

        dictionary_manager(dictionary_name).SetGrammar(dakota_options_grammar_);

        if (dictionary_manager(dictionary_name).CheckOption("@TabularDataFile"))
			dictionary_manager(dictionary_name).ReadString("@TabularDataFile", tabular_data_file_);

		if (dictionary_manager(dictionary_name).CheckOption("@Method"))
			dictionary_manager(dictionary_name).ReadString("@Method", method_);

		if (dictionary_manager(dictionary_name).CheckOption("@MaxIterations"))
			dictionary_manager(dictionary_name).ReadString("@MaxIterations", max_iterations_);

		if (dictionary_manager(dictionary_name).CheckOption("@MaxFunctionEvaluations"))
			dictionary_manager(dictionary_name).ReadString("@MaxFunctionEvaluations", max_function_evaluations_);

		if (dictionary_manager(dictionary_name).CheckOption("@ConvergenceTolerance"))
			dictionary_manager(dictionary_name).ReadString("@ConvergenceTolerance", convergence_tolerance_);

		if (dictionary_manager(dictionary_name).CheckOption("@SolutionTarget"))
			dictionary_manager(dictionary_name).ReadString("@SolutionTarget", solution_target_);

		if (dictionary_manager(dictionary_name).CheckOption("@Seed"))
			dictionary_manager(dictionary_name).ReadString("@Seed", seed_);

		if (dictionary_manager(dictionary_name).CheckOption("@PopulationSize"))
			dictionary_manager(dictionary_name).ReadString("@PopulationSize", population_size_);

		if (dictionary_manager(dictionary_name).CheckOption("@FitnessType"))
			dictionary_manager(dictionary_name).ReadString("@FitnessType", fitness_type_);

		if (dictionary_manager(dictionary_name).CheckOption("@MutationType"))
			dictionary_manager(dictionary_name).ReadString("@MutationType", mutation_type_);

		if (dictionary_manager(dictionary_name).CheckOption("@MutationRate"))
			dictionary_manager(dictionary_name).ReadString("@MutationRate", mutation_rate_);

		if (dictionary_manager(dictionary_name).CheckOption("@CrossoverType"))
			dictionary_manager(dictionary_name).ReadString("@CrossoverType", crossover_type_);

		if (dictionary_manager(dictionary_name).CheckOption("@CrossoverRate"))
			dictionary_manager(dictionary_name).ReadString("@CrossoverRate", crossover_rate_);

		if (dictionary_manager(dictionary_name).CheckOption("@ReplacementType"))
			dictionary_manager(dictionary_name).ReadString("@ReplacementType", replacement_type_);

		if (dictionary_manager(dictionary_name).CheckOption("@Division"))
			dictionary_manager(dictionary_name).ReadString("@Division", division_);

		if (dictionary_manager(dictionary_name).CheckOption("@MaxBoxsizeLimit"))
			dictionary_manager(dictionary_name).ReadString("@MaxBoxsizeLimit", max_boxsize_limit_);

		if (dictionary_manager(dictionary_name).CheckOption("@MinBoxsizeLimit"))
			dictionary_manager(dictionary_name).ReadString("@MinBoxsizeLimit", min_boxsize_limit_);

		if (dictionary_manager(dictionary_name).CheckOption("@DiverseInput")){
			diverse_input_ = true;
			dictionary_manager(dictionary_name).ReadOption("@DiverseInput", diverse_dakota_input_);
		}

 		if (dictionary_manager(dictionary_name).CheckOption("@Gradient"))
            dictionary_manager(dictionary_name).ReadBool("@Gradient", dakota_gradient_);

		if(dictionary_manager(dictionary_name).CheckOption("@EchoDakotaInput"))
			dictionary_manager(dictionary_name).ReadBool("@EchoDakotaInput", echo_dakota_string_);
	}
}