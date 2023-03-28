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
    options_optimization_setup::options_optimization_setup()
    {
        penalty_function_ = true;
        iReactionClasses_ = false;
    }
    
    options_optimization_setup::~options_optimization_setup(){}

    void options_optimization_setup::SetupFromDictionary
                (OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                std::string dictionary_name)
    {
        dictionary_manager(dictionary_name).SetGrammar(optimization_setup_grammar_);
        
        if (dictionary_manager(dictionary_name).CheckOption("@ParametersBoundaries"))
			dictionary_manager(dictionary_name).ReadString("@ParametersBoundaries", parameter_boundaries_);

        if (dictionary_manager(dictionary_name).CheckOption("@SigmaExpDistribution"))
			dictionary_manager(dictionary_name).ReadInt("@SigmaExpDistribution", sigma_exp_ditribution_);

        if (dictionary_manager(dictionary_name).CheckOption("@AcceptedSigmaInKDistribution"))
			dictionary_manager(dictionary_name).ReadInt("@AcceptedSigmaInKDistribution", sigma_k_distribution_);

        if (dictionary_manager(dictionary_name).CheckOption("@Parameters_Distribution"))
			dictionary_manager(dictionary_name).ReadString("@Parameters_Distribution", parameter_distribution_);

        if (dictionary_manager(dictionary_name).CheckOption("@PenaltyFunction"))
			dictionary_manager(dictionary_name).ReadBool("@PenaltyFunction", penalty_function_);

        if (dictionary_manager(dictionary_name).CheckOption("@ObjectiveFunctionType"))
			dictionary_manager(dictionary_name).ReadString("@ObjectiveFunctionType", objective_function_type_);

        if (dictionary_manager(dictionary_name).CheckOption("@ReactionsClasses"))
			dictionary_manager(dictionary_name).ReadBool("@ReactionsClasses", iReactionClasses_);
        
        // if(penalty_function_ == false)
        // {
        //     std::cout << " * Penalty Function is turned off!" << std::endl;
        // }
    }
} // namespace OptiSMOKE
