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

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
	class Grammar_PlugFlowReactor_Options : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OutputFolder", 
																OpenSMOKE::SINGLE_PATH, 
																"Name of the folder where to write the output data", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@StepsVideo", 
																OpenSMOKE::SINGLE_INT, 
																"Parameter governing the frequency of output on video", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@StepsFile", 
																OpenSMOKE::SINGLE_INT, 
																"Parameter governing the frequency of output on file", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OutputSpecies", 
																OpenSMOKE::VECTOR_STRING, 
																"List of species which will be written on ASCII file", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Verbose", 
																OpenSMOKE::SINGLE_BOOL, 
																"If set false means that video info and output files will not be written", 
																false) );

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@VerboseVideo",
															   OpenSMOKE::SINGLE_BOOL,
															   "If set false means that output video file will not be written",
															   false));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@VerboseASCIIFile", 
																OpenSMOKE::SINGLE_BOOL, 
																"If set false means that output ASCII file will not be written", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@VerboseXMLFile", 
																OpenSMOKE::SINGLE_BOOL, 
																"If set false means that output XML file will not be written", 
																false) );
		}
	};

	PlugFlowReactor_Options::PlugFlowReactor_Options()
	{
		verbose_output_ = true;
		verbose_video_ = true;
		verbose_ascii_file_ = true;
		verbose_xml_file_ = true;
		output_path_ = "Output";
		n_step_video_ = 50;
		n_step_file_ = 5;
		sensitivity_analysis_ = false;
	}

	void PlugFlowReactor_Options::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_PlugFlowReactor_Options grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@OutputFolder") == true)
			dictionary.ReadPath("@OutputFolder", output_path_);

		if (dictionary.CheckOption("@StepsVideo") == true)
			dictionary.ReadInt("@StepsVideo", n_step_video_);

		if (dictionary.CheckOption("@StepsFile") == true)
			dictionary.ReadInt("@StepsFile", n_step_file_);

		if (dictionary.CheckOption("@OutputSpecies") == true)
			dictionary.ReadOption("@OutputSpecies", output_species_);

		if (dictionary.CheckOption("@Verbose") == true)
			dictionary.ReadBool("@Verbose", verbose_output_);

		if (dictionary.CheckOption("@VerboseVideo") == true)
			dictionary.ReadBool("@VerboseVideo", verbose_video_);

		if (dictionary.CheckOption("@VerboseASCIIFile") == true)
			dictionary.ReadBool("@VerboseASCIIFile", verbose_ascii_file_);

		if (dictionary.CheckOption("@VerboseXMLFile") == true)
			dictionary.ReadBool("@VerboseXMLFile", verbose_xml_file_);

		// Cross-Check
		if (verbose_xml_file_ == true)
			verbose_output_ = true;

		// Cross-Check
		if (verbose_ascii_file_ == true)
			verbose_output_ = true;

		// Cross-Check
		if (sensitivity_analysis_ == true)
		{
			verbose_output_ = true;
			verbose_xml_file_ = true;
		}
	}

}
