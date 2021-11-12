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

#ifndef OpenSMOKE_PlugFlowReactor_Options_H
#define	OpenSMOKE_PlugFlowReactor_Options_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
	class PlugFlowReactor_Options
	{
	public:

		/**
		*@brief Default constructor
		*/
		PlugFlowReactor_Options();

		/**
		*@brief Initializes the object from an external dictionary
		*@param dictionary external dictionary
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Sets the verbose output
		*/
		void SetVerboseOutput(const bool verbose_output)				{ verbose_output_ = verbose_output; }

		/**
		*@brief Sets the verbose output
		*/
		void SetVerboseVideo(const bool verbose_video)				{ verbose_video_ = verbose_video; }

		/**
		*@brief Sets the verbose ASCII output
		*/
		void SetVerboseASCIIFile(const bool verbose_ascii_file)			{ verbose_ascii_file_ = verbose_ascii_file; }

		/**
		*@brief Sets the verbose XML output
		*/
		void SetVerboseXMLFile(const bool verbose_xml_file)				{ verbose_xml_file_ = verbose_xml_file; }
		
		/**
		*@brief Sets the output folder
		*/
		void SetOutputPath(const boost::filesystem::path output_path)	{ output_path_ = output_path; }

		/**
		*@brief Sets the output folder
		*/
		void SetOutputPath(const std::string output_path)				{ output_path_ = output_path; }

		/**
		*@brief Sets the parameter governing the output flow on the video
		*/
		void SetNumberOfSteps_Video(unsigned int n_step_video)		{ n_step_video_  = n_step_video; }

		/**
		*@brief Sets the parameter governing the output flow on files
		*/
		void SetNumberOfSteps_File(unsigned int n_step_file)		{ n_step_file_   = n_step_file; }

		/**
		*@brief Sets the sensitivity analysis
		*/
		void SetSensitivityAnalysis(const bool sensitivity_analysis)	{ sensitivity_analysis_   = sensitivity_analysis; }
		
		/**
		*@brief Sets the species for which the profiles will be written on output
		*/		
		void SetOutputSpecies(std::vector<std::string>& output_species)				{ output_species_   = output_species; }

		/**
		*@brief Returns the verbose output
		*/
		bool verbose_output() const				{ return verbose_output_; }

		/**
		*@brief Returns the verbose video
		*/
		bool verbose_video() const				{ return verbose_video_; }

		/**
		*@brief Returns the verbose ASCII output
		*/
		bool verbose_ascii_file() const			{ return verbose_ascii_file_; }

		/**
		*@brief Returns the verbose XML output
		*/
		bool verbose_xml_file() const			{ return verbose_xml_file_; }

		/**
		*@brief Returns the output folder
		*/
		boost::filesystem::path output_path() const		{ return output_path_; }

		/**
		*@brief Returns the parameter governing the output flow on the video
		*/
		unsigned int n_step_video() const				{ return n_step_video_; }

		/**
		*@brief Returns the parameter governing the output flow on the files
		*/
		unsigned int n_step_file() const				{ return n_step_file_; }

		/**
		*@brief Returns the sensitivity analysis
		*/
		bool sensitivity_analysis() const		{ return sensitivity_analysis_; }

		/**
		*@brief Returns the species for which the profiles will be written on output
		*/
		const std::vector<std::string>& output_species() const			{ return output_species_; }

	private:

		bool verbose_output_;								//!< verbose output (both files and video)
		bool verbose_video_;								//!< verbose output (only video)
		bool verbose_ascii_file_;							//!< verbose ASCII file
		bool verbose_xml_file_;								//!< verbose XML file
		boost::filesystem::path output_path_;				//!< output folder where the output data will be written
		int n_step_video_;									//!< parameter governing the output flow on the video
		int n_step_file_;									//!< parameter governing the output flow on the file
		bool sensitivity_analysis_;							//!< sensitivity analysis 
		std::vector<std::string> output_species_;			//!< list of species for which the profiles will be written on output
	};
}

#include "PlugFlowReactor_Options.hpp"

#endif	/* OpenSMOKE_PlugFlowReactor_Options_H */

