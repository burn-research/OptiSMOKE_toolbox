/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2016  Alberto Cuoci                                      |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_Grammar_LewisNumbers_H
#define	OpenSMOKE_Grammar_LewisNumbers_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{

	class Grammar_LewisNumbers : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Default",
							OpenSMOKE::SINGLE_DOUBLE,
							"Default Lewis numbers for species not available in the provided list",
							true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Species",
							OpenSMOKE::VECTOR_STRING_DOUBLE,
							"List of Lewis numbers for each species (@Species CH4 0.90 N2 1.00;)",
							false ));
		}
	};

	void GetLewisNumbersFromDictionary(	OpenSMOKE::OpenSMOKE_Dictionary& dictionary, 
										OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXML,
										std::vector<double>& lewis_numbers)
	{
		// Set the grammar
		Grammar_LewisNumbers grammar_lewis_numbers;
		dictionary.SetGrammar(grammar_lewis_numbers);

		// Resize the vector
		lewis_numbers.resize(thermodynamicsMapXML.NumberOfSpecies());

		// Read the default value
		double default_value;
		dictionary.ReadDouble("@Default", default_value);
		for (unsigned int i = 0; i<thermodynamicsMapXML.NumberOfSpecies(); i++)
			lewis_numbers[i] = default_value;

		// Read the list of Lewis numbers
		if (dictionary.CheckOption("@Species"))
		{
			std::vector<std::string> names;
			std::vector<double> values;

			dictionary.ReadOption("@Species", names, values);
			for (unsigned int i = 0; i < values.size(); i++)
			{
				if (values[i] <= 0.)
					OpenSMOKE::FatalErrorMessage("The provided Lewis numbers must be strictly positive");

				lewis_numbers[thermodynamicsMapXML.IndexOfSpecies(names[i])-1] = values[i];
			}
		}

		// Lewis numbers on the screen
		{
			std::cout << std::endl;
			std::cout << "------------------------------------------" << std::endl;
			std::cout << "                Lewis numbers             " << std::endl;
			std::cout << "------------------------------------------" << std::endl;
			
			std::cout << std::setw(20) << std::left << " * Default" << std::fixed << std::setprecision(4) << default_value << std::endl;
			for (unsigned int i = 0; i < thermodynamicsMapXML.NumberOfSpecies(); i++)
			{
				if (lewis_numbers[i] != default_value)
				{
					std::string tag = " * " + thermodynamicsMapXML.NamesOfSpecies()[i];
					std::cout << std::setw(20) << std::left << tag << std::fixed << std::setprecision(4) << lewis_numbers[i] << std::endl;
				}
			}

			std::cout << "------------------------------------------" << std::endl;
			std::cout << std::endl;
		}
	}
}

#endif	/* OpenSMOKE_Grammar_LewisNumbers_H */

