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
    options_curvematching::options_curvematching(){
        use_bootstrap_ = false;
        number_of_bootstrap_ = 0;
    }
    
    options_curvematching::~options_curvematching() {}

    void options_curvematching::SetupFromDictionary
                                (OpenSMOKE::OpenSMOKE_DictionaryManager& dictionary_manager, 
                                std::string dictionary_name) 
    {
        dictionary_manager(dictionary_name).SetGrammar(grammar_curve_matching_);

        if (dictionary_manager(dictionary_name).CheckOption("@NumberOfBootstrapVariations"))
			dictionary_manager(dictionary_name).ReadInt("@NumberOfBootstrapVariations", number_of_bootstrap_);

        if (dictionary_manager(dictionary_name).CheckOption("@UseBootStrap"))
			dictionary_manager(dictionary_name).ReadBool("@UseBootStrap", use_bootstrap_);

        CheckCurveMatchingOptions();
    }
    
    void options_curvematching::CheckCurveMatchingOptions(){
        if(use_bootstrap_){
            if(number_of_bootstrap_ < 2)
                OptiSMOKE::FatalErrorMessage("Use bootstrap true implies that the number of bootstrap variations is greter than one!");
        }
    }


} // namespace OptiSMOKE
