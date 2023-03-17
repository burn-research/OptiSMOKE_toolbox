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
                                          
#ifndef GRAMMAR_DAKOTAOPTIONS_H
#define GRAMMAR_DAKOTAOPTIONS_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OptiSMOKE
{
	class grammar_dakota : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TabularDataFile",
                OpenSMOKE::SINGLE_STRING,
                "Specific name for tabular data file",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Method",
                OpenSMOKE::SINGLE_STRING,
                "Method for optimization",
                true));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaxIterations",
                OpenSMOKE::SINGLE_STRING,
                "Maximum number of iterations",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaxFunctionEvaluations",
                OpenSMOKE::SINGLE_STRING,
                "Maximum number of function evaluations",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ConvergenceTolerance",
                OpenSMOKE::SINGLE_STRING,
                "Convergence tolerance",
                false));
	
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SolutionTarget",
                OpenSMOKE::SINGLE_STRING,
                "Solution target for the objective function",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Seed",
                OpenSMOKE::SINGLE_STRING,
                "Seed for the random sampling",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PopulationSize",
                OpenSMOKE::SINGLE_STRING,
                "Population size for the random samples (EA)",
                false));
		
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FitnessType",
                OpenSMOKE::SINGLE_STRING,
                "Deciding the fitness type for the objective function (linear_rank | merit_function)",
                false));
		
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MutationType",
                OpenSMOKE::SINGLE_STRING,
                "The mutation_type controls what approach is employed in randomly modifying continuous design variables within the EA population.",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MutationRate",
                OpenSMOKE::SINGLE_STRING,
                "The mutation_rate controls the probability of mutation being performed on an individual, both for new individuals generated by crossover (if crossover occurs) and for individuals from the existing population. It is the fraction of trial points that are mutated in a given iteration and therefore must be specified to be between 0 and 1.",
                false,"none","none","none"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CrossoverType",
                OpenSMOKE::SINGLE_STRING,
                "The crossover_type controls what approach is employed for combining parent genetic information to create offspring. The SCOLIB EA method supports three forms of crossover, two_point, blend, and uniform, which generate a new individual through combinations of two parent individuals.",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CrossoverRate",
                OpenSMOKE::SINGLE_STRING,
                "The crossover_type controls what approach is employed for combining parent genetic information to create offspring, and the crossover_rate specifies the probability of a crossover operation being performed to generate a new offspring.",
                false,"none","none","none"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ReplacementType",
                OpenSMOKE::SINGLE_STRING,
                "The replacement_type controls how current populations and newly generated individuals are combined to create a new population. Each of the replacement_type selections accepts an associated integer value. (random | chc | elitist)",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Division",
                OpenSMOKE::SINGLE_STRING,
                "The division specification determines how DIRECT subdivides each subregion of the search space. (major_dimension | all_dimensions)",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaxBoxsizeLimit",
                OpenSMOKE::SINGLE_STRING,
                "max_boxsize_limit specification terminates DIRECT if the size of the largest subregion falls below this threshold.",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MinBoxsizeLimit",
                OpenSMOKE::SINGLE_STRING,
                "min_boxsize_limit specification terminates DIRECT if the size of the smallest subregion falls below this threshold.",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@GlobalBalanceParameter",
                OpenSMOKE::SINGLE_STRING,
                "The global_balance_parameter controls how much global search is performed by only allowing a subregion to be subdivided if the size of the subregion divided by the size of the largest subregion is at least global_balance_parameter. Intuitively, this forces large subregions to be subdivided before the smallest subregions are refined.",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DiverseInput",
                OpenSMOKE::VECTOR_STRING,
                "List of diverse inputs that can be given to other Dakota optimization methods.",
                false));

	        AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Gradient",
                OpenSMOKE::SINGLE_BOOL,
                "Activate gradient for the Optimization method (default: false)",
                false));

		}
	};
}

#endif // GRAMMAR_DAKOTAOPTIONS_H