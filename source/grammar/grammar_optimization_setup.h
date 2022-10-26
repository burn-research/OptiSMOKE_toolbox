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
|                     Timoteo Dinelli <timoteo.dinelli@polimi.it          |
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

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OptiSMOKE
{
	class grammar_optimization_setup : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SigmaExpDistribution",
                OpenSMOKE::SINGLE_INT,
                "How many sigma is the relative error",
                false,
                "none",
                "none",
                "none"));
	
	        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@AcceptedSigmaInKDistribution",
                OpenSMOKE::SINGLE_INT,
                "How many sigma we assume to accept in our k distribution",
                false,
                "none",
                "none",
                "none"));

    	    AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Parameters_Distribution",
                OpenSMOKE::SINGLE_STRING,
                "Specifies the parameters distribution: Uniform or Normal",
                true,
                "none",
                "none",
                "none"));

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PenaltyFunction",
                OpenSMOKE::SINGLE_BOOL,
		        "Use penalty function for checking the rate parameters (default: true)",
		        false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ObjectiveFunctionType",
                OpenSMOKE::SINGLE_STRING,
                "Specify how the objective function is calculated.",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ParametersBoundaries",
                OpenSMOKE::SINGLE_STRING,
                "Specify the method to compute the parameters boundaries",
                false));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ScalingReactionClasses",
                OpenSMOKE::SINGLE_BOOL,
                "If the reaction class optimization has to be done with scaling or with substitution",
                false,
                "none",
                "@ReactionsClasses",
                "none"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ReactionsClasses",
                OpenSMOKE::SINGLE_BOOL,
                "Enable the optimisation by reaction classes",
                false, 
                "none",
		        "@ReactionsClassesDefinitions",
		        "none"));
        }
    };
}