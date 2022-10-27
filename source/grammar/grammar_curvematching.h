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

#ifndef GRAMMAR_CURVEMATCHINGOPTIONS_H
#define GRAMMAR_CURVEMATCHINGOPTIONS_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OptiSMOKE
{
	class grammar_curve_matching : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
    
            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfBootstrapVariations",
		        OpenSMOKE::SINGLE_INT,
	            "Number of Bootstraps variations",
                false,
                "none",
                "none",
                "none"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@LineUpMaxima",
                OpenSMOKE::SINGLE_BOOL,
	            "Number of Bootstraps variations",
                false,
                "none",
                "none",
                "none"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@UseSumOfIndexesForAlignment",
                OpenSMOKE::SINGLE_BOOL,
	            "Number of Bootstraps variations",
                false,
                "none",
                "none",
                "none"));

            AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FractionOfExpRangeForModelExtrapolation",
                OpenSMOKE::VECTOR_DOUBLE,
	            "Number of Bootstraps variations",
                false,
                "none",
                "none",
                "none"));

            // BOOTSTRAP//		
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@UseBootStrap",
                OpenSMOKE::SINGLE_BOOL,
                "Use Bootstrap technique in Curve Matching Index calculations (default: false)",
                false) );

            // Da vedere se metterli
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PrintIndexes",
                OpenSMOKE::SINGLE_BOOL,
                "Decide wheter to print out the indexes (default: false)",
                false));

	        AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PrintSplines",
                OpenSMOKE::SINGLE_BOOL,
                "Decide wheter to print out the splines (default: false)",
                false));

	        AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PrintBootstrap",
                OpenSMOKE::SINGLE_BOOL,
                "Decide wheter to print out the bootstrap  (default: false)",
                false));
		}
	};
}

#endif // GRAMMAR_CURVEMATCHINGOPTIONS_H