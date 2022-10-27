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

#ifndef OPTIONS_OPTIMIZATION_SETUP_H
#define OPTIONS_OPTIMIZATION_SETUP_H

namespace OptiSMOKE
{
    class options_optimization_setup
    {

    public:
        options_optimization_setup();
        
        ~options_optimization_setup();
        
        void SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                                std::string dictionary_name);
    
        inline const std::string& parameter_boundaries() {return parameter_boundaries_;};
        inline const int& sigma_exp_ditribution() {return sigma_exp_ditribution_;};
        inline const std::string& sigma_k_distribution() {return sigma_k_distribution_;};
        inline const std::string& parameter_distribution() {return parameter_distribution_;};
        inline const bool& penalty_function() {return penalty_function_;};
        inline const std::string& objective_function_type() {return objective_function_type_;};
        inline const bool& iReactionClasses() {return iReactionClasses_;};
        inline const bool& iScalingClasses() {return iScalingClasses_;};

    private:

        grammar_optimization_setup optimization_setup_grammar_;
        std::string parameter_boundaries_;
        int sigma_exp_ditribution_;
        std::string sigma_k_distribution_;
        std::string parameter_distribution_;
        bool penalty_function_;
        std::string objective_function_type_;
        bool iReactionClasses_;
        bool iScalingClasses_;
    };
    
} // namespace OptiSMOKE

#include "options_optimization_setup.hpp"
#endif // OPTIONS_OPTIMIZATION_SETUP_H