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

#ifndef OPTIONS_CURVEMATCHING_H
#define OPTIONS_CURVEMATCHING_H

namespace OptiSMOKE
{
    class options_curvematching
    {  
    public:
        
        options_curvematching();
        ~options_curvematching();
        void SetupFromDictionary(OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager,
                                std::string dictionary_name);

        inline const int& number_of_bootstrap() {return number_of_bootstrap_;};
        inline const bool& line_up_maxima() {return line_up_maxima_;};
        inline const bool& use_index_for_alignement() {return use_index_for_alignement_;};
        inline const bool& use_bootstrap() {return use_bootstrap_;};
        inline const double& fraction_of_exp_for_model_extrapolation() {return fraction_of_exp_for_model_extrapolation_;};

    private:

        grammar_curve_matching grammar_curve_matching_;
        int number_of_bootstrap_;
        bool line_up_maxima_;
        bool use_index_for_alignement_;
        bool use_bootstrap_;
        double fraction_of_exp_for_model_extrapolation_;
        
    };
} // namespace OptiSMOKE

#include "options_curvematching.hpp"
#endif // OPTIONS_CURVEMATCHING_H