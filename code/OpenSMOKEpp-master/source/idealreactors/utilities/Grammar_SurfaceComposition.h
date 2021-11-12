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
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
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

#ifndef OpenSMOKE_Grammar_SurfaceComposition_H
#define	OpenSMOKE_Grammar_SurfaceComposition_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{

	class Grammar_SurfaceComposition : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SurfaceFractions", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Site fractions of the surface (i.e. OH(s) 0.60 Rh(s) 0.40)", 
																true,
																"@SurfaceCompositions",
																"none",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SurfaceCompositions", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Mass fractions of the mixture (i.e. CH4 6 H2 4)", 
																true,
																"@SurfaceFractions",
																"none",
																"none") );

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SiteNonConservation",
															   OpenSMOKE::SINGLE_BOOL,
															   "If true site non conservation is allowed. If false, site density is forced to stay constant. Default: true",
															   false));
		}
	};

	template<typename ThermodynamicsMap>
	void GetSurfaceCompositionFromDictionary(	OpenSMOKE::OpenSMOKE_Dictionary& dictionary, ThermodynamicsMap& thermodynamicsMap, OpenSMOKE::OpenSMOKEVectorDouble& Z, bool& site_non_conservation)
	{
		Grammar_SurfaceComposition grammar_surface_status;
		dictionary.SetGrammar(grammar_surface_status);
		
		std::string name_of_site_phase = dictionary.name().substr(8, dictionary.name().size()-8);

		int index_site_phase = -1;
		for(unsigned int k=0;k<thermodynamicsMap.matrix_names_site_phases()[0].size();k++)
			if (thermodynamicsMap.matrix_names_site_phases()[0][k] == name_of_site_phase)
				index_site_phase = k;

		// Site non conservation
		{
			site_non_conservation = true;
			if (dictionary.CheckOption("@SiteNonConservation") == true)
				dictionary.ReadBool("@SiteNonConservation", site_non_conservation);
		}

		// Composition
		{
			std::vector<std::string> names;
			std::vector<double> values;
		
			if (dictionary.CheckOption("@SurfaceFractions") == true)
			{
				dictionary.ReadOption("@SurfaceFractions", names, values);
				
				const double sum =std::accumulate(values.begin(),values.end(),0.);
				if (sum<(1.-1e-6) || sum>(1.+1e-6))
					OpenSMOKE::FatalErrorMessage("The surface fractions must sum to 1.");

				//ChangeDimensions(thermodynamicsMap.number_of_site_species(), &Z, true);
				for(unsigned int i=0;i<names.size();i++)
				{
					const unsigned int index = thermodynamicsMap.IndexOfSpecies(names[i]) -  thermodynamicsMap.number_of_gas_species();
					if (thermodynamicsMap.vector_site_phases_belonging()[index-1] != index_site_phase)
					{
						std::string message = "The following species do not belong to the " + name_of_site_phase + " site phase:" + names[i];
						OpenSMOKE::FatalErrorMessage(message);
					}
					Z[index] = values[i]/sum;
				}
			}
			else if (dictionary.CheckOption("@SurfaceCompositions") == true)
			{
				dictionary.ReadOption("@SurfaceCompositions", names, values);
				const double sum =std::accumulate(values.begin(),values.end(),0.);

				//ChangeDimensions(thermodynamicsMap.number_of_site_species(), &Z, true);
				for(unsigned int i=0;i<names.size();i++)
				{
					const unsigned int index = thermodynamicsMap.IndexOfSpecies(names[i]) -  thermodynamicsMap.number_of_gas_species();
					if (thermodynamicsMap.vector_site_phases_belonging()[index] != index_site_phase)
					{
						std::string message = "The following species do not belong to the " + name_of_site_phase + " site phase:" + names[i];
						OpenSMOKE::FatalErrorMessage(message);
					}
					Z[index] = values[i]/sum;
				}
			}
		}
	}
}

#endif	/* OpenSMOKE_Grammar_SurfaceStatus_Options_H */

