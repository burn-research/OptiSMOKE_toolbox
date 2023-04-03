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
|                     Timoteo Dinelli <timoteo.dinelli@polimi.it>         |
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

#ifndef GRAMMAR_OPTISMOKEPP_H
#define GRAMMAR_OPTISMOKEPP_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OptiSMOKE
{
	class grammar_optismoke : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

	virtual void DefineRules()
	{

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@DakotaOptions", 
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary with options for Dakota", 
                        false,
                        "none",
                        "@OptimizationLibrary",
                        "@NLOPTOptions"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NLOPTOptions", 
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary with options for NLOPT", 
                        false,
                        "none",
                        "@OptimizationLibrary",
                        "@DakotaOptions"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CurveMatchingOptions", 
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary with options for setting Curve Matching as the objective function", 
                        false));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OptimizationSetup",
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary to set the options for the optimization",
                        true));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OptimizationTarget",
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary to set the targets for the optimization",
                        true));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsFolder",
                        OpenSMOKE::SINGLE_PATH,
                        "Name of the folder containing the kinetic scheme (XML Version)",
                        true,
                        "@KineticsPreProcessor",
                        "none",
                        "none"));
                
                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor",
                        OpenSMOKE::SINGLE_DICTIONARY,
                        "Name of the dictionary containing the list of kinetic files to be interpreted",
                        true,
                        "@KineticsFolder",
                        "none",
                        "none"));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OutputFolder",
                        OpenSMOKE::SINGLE_PATH,
                        "Name of the folder where to write all the output file of the program",
                        true));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfExperimentalDataFiles",
                        OpenSMOKE::VECTOR_STRING,
                        "List of Experimental data files",
                        true));

                AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OptimizationLibrary",
                        OpenSMOKE::SINGLE_STRING,
                        "Name of the library containing the optimization routines",
                        true));                
	}
    };
}

#endif // GRAMMAR_OPTISMOKEPP_H
