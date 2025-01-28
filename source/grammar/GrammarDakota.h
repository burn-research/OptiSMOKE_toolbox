/* ------------------------------------------------------------------------------- *\
|                                                                                   |
|             ____        __  _ _____ __  _______  __ __ ______                     |
|            / __ \____  / /_(_) ___//  |/  / __ \/ //_// ____/___  ____            |
|           / / / / __ \/ __/ /\__ \/ /|_/ / / / / ,<  / __/ / __ \/ __ \           |
|          / /_/ / /_/ / /_/ /___/ / /  / / /_/ / /| |/ /___/ /_/ / /_/ /           |
|          \____/ .___/\__/_//____/_/  /_/\____/_/ |_/_____/ .___/ .___/            |
|              /_/                                        /_/   /_/                 |
|                                                                                   |
| --------------------------------------------------------------------------------- |
|  Please refer to the copyright statement and license                              |
|  information at the end of this file.                                             |
| --------------------------------------------------------------------------------- |
|                                                                                   |
|           Authors: Timoteo Dinelli  <timoteo.dinelli@polimi.it>                   |
|                    Andrea Bertolino <andrea.bertolino@ulb.be>                     |
|                    Magnus Fürst     <magnus.furst@ulb.ac.be>                      |
|                                                                                   |
|             [1] CRECK Modeling Lab <https://www.creckmodeling.polimi.it>          |
|                 Department of Chemistry, Materials and Chemical Engineering       |
|                 Politecnico di Milano, P.zza Leonardo da Vinci 32, 20133 Milano   |
|                                                                                   |
|             [2] BRITE Research Group <https://brite-research.be>                  |
|                 Brussels Institute for Thermal-fluid systems and clean Energy     |
|                 Avenue F.D. Rooseveltlaan 50, Bruxelles 1050 Brussel              |
|                                                                                   |
\* ------------------------------------------------------------------------------- */
#pragma once

namespace OptiSMOKE {
class GrammarDakota : public OpenSMOKE::OpenSMOKE_DictionaryGrammar {
 protected:
  virtual void DefineRules() {
    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TabularDataFile",
                                                      OpenSMOKE::SINGLE_STRING,
                                                      "Specific name for tabular data file",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EchoDakotaInput",
                                                      OpenSMOKE::SINGLE_BOOL,
                                                      "Print on the log file the dakota input string",
                                                      false));

    AddKeyWord(
        OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Method", OpenSMOKE::SINGLE_STRING, "Method for optimization", true));

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

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@FitnessType",
        OpenSMOKE::SINGLE_STRING,
        "Deciding the fitness type for the objective function (linear_rank | merit_function)",
        false));

    AddKeyWord(
        OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MutationType",
                                               OpenSMOKE::SINGLE_STRING,
                                               "The mutation_type controls what approach is employed in randomly "
                                               "modifying continuous design variables within the EA population.",
                                               false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@MutationRate",
        OpenSMOKE::SINGLE_STRING,
        "The mutation_rate controls the probability of mutation being performed on an individual, both for new "
        "individuals generated by crossover (if crossover occurs) and for individuals from the existing population. It "
        "is the fraction of trial points that are mutated in a given iteration and therefore must be specified to be "
        "between 0 and 1.",
        false,
        "none",
        "none",
        "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@CrossoverType",
        OpenSMOKE::SINGLE_STRING,
        "The crossover_type controls what approach is employed for combining parent genetic information to create "
        "offspring. The SCOLIB EA method supports three forms of crossover, two_point, blend, and uniform, which "
        "generate a new individual through combinations of two parent individuals.",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@CrossoverRate",
        OpenSMOKE::SINGLE_STRING,
        "The crossover_type controls what approach is employed for combining parent genetic information to create "
        "offspring, and the crossover_rate specifies the probability of a crossover operation being performed to "
        "generate a new offspring.",
        false,
        "none",
        "none",
        "none"));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@ReplacementType",
        OpenSMOKE::SINGLE_STRING,
        "The replacement_type controls how current populations and newly generated individuals are combined to create "
        "a new population. Each of the replacement_type selections accepts an associated integer value. (random | chc "
        "| elitist)",
        false));

    AddKeyWord(
        OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Division",
                                               OpenSMOKE::SINGLE_STRING,
                                               "The division specification determines how DIRECT subdivides each "
                                               "subregion of the search space. (major_dimension | all_dimensions)",
                                               false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaxBoxsizeLimit",
                                                      OpenSMOKE::SINGLE_STRING,
                                                      "max_boxsize_limit specification terminates DIRECT if the size "
                                                      "of the largest subregion falls below this threshold.",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MinBoxsizeLimit",
                                                      OpenSMOKE::SINGLE_STRING,
                                                      "min_boxsize_limit specification terminates DIRECT if the size "
                                                      "of the smallest subregion falls below this threshold.",
                                                      false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@GlobalBalanceParameter",
        OpenSMOKE::SINGLE_STRING,
        "The global_balance_parameter controls how much global search is performed by only allowing a subregion to be "
        "subdivided if the size of the subregion divided by the size of the largest subregion is at least "
        "global_balance_parameter. Intuitively, this forces large subregions to be subdivided before the smallest "
        "subregions are refined.",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        "@DiverseInput",
        OpenSMOKE::VECTOR_STRING,
        "List of diverse inputs that can be given to other Dakota optimization methods.",
        false));

    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Gradient",
                                                      OpenSMOKE::SINGLE_BOOL,
                                                      "Activate gradient for the Optimization method (default: false)",
                                                      false));
  }
};
}  // namespace OptiSMOKE

/* ------------------------------------------------------------------------------- *\
|                                                                                   |
|   MIT License                                                                     |
|                                                                                   |
|   Copyright (c) 2025 Timoteo Dinelli, Andrea Bertolino, Magnus Fürst              |
|                                                                                   |
|   Permission is hereby granted, free of charge, to any person obtaining a copy    |
|   of this software and associated documentation files (the "Software"), to deal   |
|   in the Software without restriction, including without limitation the rights    |
|   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       |
|   copies of the Software, and to permit persons to whom the Software is           |
|   furnished to do so, subject to the following conditions:                        |
|                                                                                   |
|   The above copyright notice and this permission notice shall be included in all  |
|   copies or substantial portions of the Software.                                 |
|                                                                                   |
|   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      |
|   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        |
|   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     |
|   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          |
|   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   |
|   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   |
|   SOFTWARE.                                                                       |
|                                                                                   |
\* ------------------------------------------------------------------------------- */
