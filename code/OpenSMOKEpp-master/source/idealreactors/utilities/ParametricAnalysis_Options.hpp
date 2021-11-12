/*-----------------------------------------------------------------------*\
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
#include <iostream>
#include <algorithm>
#include <vector>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"


namespace OpenSMOKE
{
	class Grammar_ParametricAnalysis_Options : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type", 
																OpenSMOKE::SINGLE_STRING, 
																"Type of parameter to be analyzed: residence-time | temperature | pressure | temperature-pressure | moles | masses", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfValues", 
																OpenSMOKE::VECTOR_STRING, 
																"List of values of parameter 1 (together with units of measure, if any)",
																true,
																"@MinimumValue @ListOfProfiles",
																"none",
																"none"));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfValues2", 
																OpenSMOKE::VECTOR_STRING, 
																"List of values of parameter 2 (together with units of measure, if any)",
																false,
																"none",
																"@ListOfValues",
																"none"));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MinimumValue", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Minimum value of parameter 1", 
																true,
																"@ListOfValues @ListOfProfiles",
																"@MaximumValue",
																"none"));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MinimumValue2", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Minimum value of parameter 2", 
																false,
																"none",
																"@MaximumValue2 @MinimumValue",
																"none"));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaximumValue", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Maximum value of parameter 1", 
																false,
																"none",
																"@NumberOfPoints",
																"none"));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaximumValue2", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Maximum value of parameter 2", 
																false,
																"none",
																"@MaximumValue @NumberOfPoints",
																"none"));

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfPoints", 
																OpenSMOKE::SINGLE_INT,
																"Number of points", 
																false,
																"none",
																"@MinimumValue",
																"none"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ListOfProfiles",
																OpenSMOKE::VECTOR_STRING,
																"List of files containing profiles",
																true,
																"@ListOfValues @MinimumValue",
																"none",
																"@Combinations"));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfThreads",
															   OpenSMOKE::SINGLE_STRING,
															   "Number of threads (in case OpenMP is enabled)",
															   false) );

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Combinations",
																OpenSMOKE::SINGLE_BOOL,
																"Combination of list of values",
																false,
																"none",
																"none",
																"@ListOfProfiles"));
		}
	};

	Eigen::VectorXd ReadListOfValues(const std::vector<std::string> values, const ParametricAnalysis_Options::ParametricAnalysisType parameter_type);
	
	double ReadScalarValue(const double input_value, const std::string units, const ParametricAnalysis_Options::ParametricAnalysisType parameter_type);

	ParametricAnalysis_Options::ParametricAnalysis_Options()
	{
		enabled_        = false;						//!< enabled on/off
		parameter_type_ = PARAMETER_TYPE_TEMPERATURE;	//!< parameter type
		number_of_points_ = 10;							//!< number of points (default 10)
		number_of_parameters_ = 0;						//!< number of parameters
		maximum_values_.resize(0);						//!< minimum value of parameter
		minimum_values_.resize(0);						//!< maximum value of parameter
		list_of_values_.resize(0,0);					//!< list of values
		number_of_threads_ = 1;							//!< number of threads (is OpenMP is enabled)
		combinations_ = false;							//!< combination of list of values
	}
	
	void ParametricAnalysis_Options::SetNumberOfPoints(const int number_of_points)
	{
		number_of_points_ = number_of_points;
	}

	void ParametricAnalysis_Options::SetNumberOfParameters(const int number_of_parameters)
	{
		number_of_parameters_ = number_of_parameters;
	}
	
	void ParametricAnalysis_Options::SetMinimumValues(const Eigen::VectorXd& minimum_values)
	{
		minimum_values_ = minimum_values;
	}	
	
	void ParametricAnalysis_Options::SetMaximumValues(const Eigen::VectorXd& maximum_values)
	{
		maximum_values_ = maximum_values;
	}

	void ParametricAnalysis_Options::SetNumberOfThreads(const int number_of_threads)
	{
		number_of_threads_ = number_of_threads;
	}

	void ParametricAnalysis_Options::SetNumberOfThreads(const std::string number_of_threads)
	{
		if (boost::iequals(number_of_threads, "max"))
		{
			#if defined(_OPENMP)
				SetNumberOfThreads(omp_get_max_threads());
			#endif
		}
		else
		{
			SetNumberOfThreads(boost::lexical_cast<int>(number_of_threads));
		}		
	}

	void ParametricAnalysis_Options::SetListOfValues(const Eigen::MatrixXd& list_of_values)
	{
		list_of_values_ = list_of_values;
		number_of_points_ = boost::lexical_cast<int>(list_of_values_.rows());
		number_of_parameters_ = boost::lexical_cast<int>(list_of_values_.cols());
		
		minimum_values_.resize(number_of_parameters_);
		maximum_values_.resize(number_of_parameters_);
		for (unsigned int i = 0; i < number_of_parameters_; i++)
		{
			minimum_values_(i) = list_of_values_.col(i).minCoeff();
			maximum_values_(i) = list_of_values_.col(i).maxCoeff();
		}
	}		

	void ParametricAnalysis_Options::AddListOfValues(const Eigen::VectorXd& list_of_values)
	{
		if (list_of_values_.rows() == 0)
		{
			Eigen::MatrixXd temp(list_of_values.size(), 1);
			for (unsigned int i = 0; i < list_of_values.size(); i++)
				temp(i, 0) = list_of_values(i);
			SetListOfValues(temp);
		}
		else
		{
			if (list_of_values.size() != number_of_points_)
				OpenSMOKE::FatalErrorMessage("ParametricAnalysis_Options: Unconsistent number of values in the list");

			Eigen::MatrixXd temp = list_of_values_;
			list_of_values_.resize(number_of_points_, number_of_parameters_ + 1);
			for (unsigned int i = 0; i < number_of_points_; i++)
				for (unsigned int j = 0; j < number_of_parameters_; j++)
					list_of_values_(i, j) = temp(i, j);
			for (unsigned int i = 0; i < number_of_points_; i++)
				list_of_values_(i, number_of_parameters_) = list_of_values(i);

			number_of_parameters_++;

			minimum_values_.resize(number_of_parameters_);
			maximum_values_.resize(number_of_parameters_);
			for (unsigned int i = 0; i < number_of_parameters_; i++)
			{
				minimum_values_(i) = list_of_values_.col(i).minCoeff();
				maximum_values_(i) = list_of_values_.col(i).maxCoeff();
			}
		}
	}

	void ParametricAnalysis_Options::AddMinimumValue(const double value)
	{
		Eigen::VectorXd temp = minimum_values_;
		minimum_values_.resize(minimum_values_.size() + 1);
		for (unsigned int i = 0; i < temp.size(); i++)
			minimum_values_(i) = temp(i);

		minimum_values_(minimum_values_.size()-1) = value;
	}

	void ParametricAnalysis_Options::AddMaximumValue(const double value)
	{
		Eigen::VectorXd temp = maximum_values_;
		maximum_values_.resize(maximum_values_.size() + 1);
		for (unsigned int i = 0; i < temp.size(); i++)
			maximum_values_(i) = temp(i);

		maximum_values_(maximum_values_.size() - 1) = value;
	}

	void ParametricAnalysis_Options::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_ParametricAnalysis_Options grammar;

		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@Type") == true)
		{
			std::string value;
			dictionary.ReadString("@Type", value);
			if (value == "temperature")
				parameter_type_ = PARAMETER_TYPE_TEMPERATURE;
			else if (value == "pressure")
				parameter_type_ = PARAMETER_TYPE_PRESSURE;
			else if (value == "residence-time")
				parameter_type_ = PARAMETER_TYPE_TIME;
			else if (value == "masses")
				parameter_type_ = PARAMETER_TYPE_MASSES;
			else if (value == "moles")
				parameter_type_ = PARAMETER_TYPE_MOLES;
			else if (value == "temperature-pressure")
				parameter_type_ = PARAMETER_TYPE_TEMPERATURE_PRESSURE;
			else
				OpenSMOKE::FatalErrorMessage("ParametricAnalysis_Options: Unknown parametric type: " + value);
		}

		if (dictionary.CheckOption("@Combinations") == true)
			dictionary.ReadBool("@Combinations", combinations_);

		if (parameter_type_ == PARAMETER_TYPE_MASSES || parameter_type_ == PARAMETER_TYPE_MOLES)
		{
			if (dictionary.CheckOption("@ListOfValues") == true)
			{
				std::vector<std::string> values;
				dictionary.ReadOption("@ListOfValues", values);

				number_of_points_ = 1;
				for (unsigned int i = 0; i < values.size(); i++)
				if (values[i] == ":")
					number_of_points_++;

				list_of_species_.resize(number_of_points_);
				list_of_compositions_.resize(number_of_points_);

				unsigned int index = 0;
				unsigned int i = 0;
				for (;;)
				{
					if (i >= values.size())
						break;

					if (values[i] == ":")
					{
						index++;
						i++;
						continue;
					}
					else
					{
						list_of_species_[index].push_back(values[i]);
						list_of_compositions_[index].push_back(boost::lexical_cast<double>(values[i + 1]));
						i += 2;
					}
				}

				// Print on the screen the list of values to be analyzed
				std::cout << std::endl;
				std::cout << "-----------------------------------------------------------" << std::endl;
				std::cout << "Parametric Analysis:" << std::endl;
				std::cout << "-----------------------------------------------------------" << std::endl;
				for (int i = 0; i < number_of_points_; i++)
				{
					std::cout << i + 1 << "\t";
					for (unsigned int j = 0; j < list_of_species_[i].size(); j++)
						std::cout << list_of_species_[i][j] << " " << list_of_compositions_[i][j] << " ";
					std::cout << std::endl;
				}
				std::cout << "-----------------------------------------------------------" << std::endl;
				std::cout << std::endl;
			}
			else
				OpenSMOKE::FatalErrorMessage("ParametricAnalysis_Options: Parameteric compositions can be specified only using the @ListOfValues options");
		}
		else
		{
			if (dictionary.CheckOption("@ListOfValues") == true)
			{
				std::vector<std::string> values;
				dictionary.ReadOption("@ListOfValues", values);
				
				Eigen::VectorXd list_of_values;
				if (parameter_type_ == PARAMETER_TYPE_TEMPERATURE ||
					parameter_type_ == PARAMETER_TYPE_TEMPERATURE_PRESSURE)
					list_of_values = ReadListOfValues(values, PARAMETER_TYPE_TEMPERATURE);
				else if (parameter_type_ == PARAMETER_TYPE_PRESSURE)
					list_of_values = ReadListOfValues(values, PARAMETER_TYPE_PRESSURE);
				else if (parameter_type_ == PARAMETER_TYPE_TIME)
					list_of_values = ReadListOfValues(values, PARAMETER_TYPE_TIME);
				
				AddListOfValues(list_of_values);
				
				// Second set of parameters
				if (dictionary.CheckOption("@ListOfValues2") == true)
				{
					std::vector<std::string> values;
					dictionary.ReadOption("@ListOfValues2", values);

					Eigen::VectorXd list_of_values;
					if (parameter_type_ == PARAMETER_TYPE_PRESSURE ||
						parameter_type_ == PARAMETER_TYPE_TEMPERATURE_PRESSURE)
						list_of_values = ReadListOfValues(values, PARAMETER_TYPE_PRESSURE);
					else if (parameter_type_ == PARAMETER_TYPE_TIME)
						list_of_values = ReadListOfValues(values, PARAMETER_TYPE_TIME);

					AddListOfValues(list_of_values);
				}
			}
			else if (dictionary.CheckOption("@ListOfProfiles") == true)
			{
				std::vector<std::string> values;
				dictionary.ReadOption("@ListOfProfiles", values);

				if (parameter_type_ == PARAMETER_TYPE_TEMPERATURE_PRESSURE)
				{
					list_of_profiles_.resize(values.size());
					for (unsigned int i = 0; i < values.size(); i++)
					{
						std::cout << "Reading profile in file: " << values[i] << std::endl;
						list_of_profiles_[i] = new FixedProfileEnriched(values[i]);
					}

					Eigen::MatrixXd temp(values.size(), 2);
					for (unsigned int i = 0; i < values.size(); i++)
					{
						temp(i, 0) = list_of_profiles_[i]->temperature();
						temp(i, 1) = list_of_profiles_[i]->pressure();
					}

					SetListOfValues(temp);
				}
				else
					OpenSMOKE::FatalErrorMessage("ParametricAnalysis_Options: The @ListOfProfiles option can be used (currently) only with @Type equal to 'temperature-pressure'");
			}
			else
			{
				if (dictionary.CheckOption("@NumberOfPoints") == true)
				{
					dictionary.ReadInt("@NumberOfPoints", number_of_points_);
					if (number_of_points_ <= 1)
						OpenSMOKE::FatalErrorMessage("The number of points must be at least equal to 2");
				}

				if (dictionary.CheckOption("@MinimumValue") == true)
				{
					std::string units;
					double value;
					dictionary.ReadMeasure("@MinimumValue", value, units);

					if (parameter_type_ == PARAMETER_TYPE_TEMPERATURE ||
						parameter_type_ == PARAMETER_TYPE_TEMPERATURE_PRESSURE)
						value = ReadScalarValue(value, units, PARAMETER_TYPE_TEMPERATURE);
					else if (parameter_type_ == PARAMETER_TYPE_PRESSURE)
						value = ReadScalarValue(value, units, PARAMETER_TYPE_PRESSURE);
					else if (parameter_type_ == PARAMETER_TYPE_TIME)
						value = ReadScalarValue(value, units, PARAMETER_TYPE_TIME);
					
					AddMinimumValue(value);

					if (dictionary.CheckOption("@MinimumValue2") == true)
					{
						std::string units;
						double value;
						dictionary.ReadMeasure("@MinimumValue2", value, units);

						if (parameter_type_ == PARAMETER_TYPE_PRESSURE ||
							parameter_type_ == PARAMETER_TYPE_TEMPERATURE_PRESSURE)
							value = ReadScalarValue(value, units, PARAMETER_TYPE_PRESSURE);
						else if (parameter_type_ == PARAMETER_TYPE_TIME)
							value = ReadScalarValue(value, units, PARAMETER_TYPE_TIME);

						AddMinimumValue(value);
					}
				}

				if (dictionary.CheckOption("@MaximumValue") == true)
				{
					std::string units;
					double value;
					dictionary.ReadMeasure("@MaximumValue", value, units);

					if (parameter_type_ == PARAMETER_TYPE_TEMPERATURE ||
						parameter_type_ == PARAMETER_TYPE_TEMPERATURE_PRESSURE)
						value = ReadScalarValue(value, units, PARAMETER_TYPE_TEMPERATURE);
					else if (parameter_type_ == PARAMETER_TYPE_PRESSURE)
						value = ReadScalarValue(value, units, PARAMETER_TYPE_PRESSURE);
					else if (parameter_type_ == PARAMETER_TYPE_TIME)
						value = ReadScalarValue(value, units, PARAMETER_TYPE_TIME);

					AddMaximumValue(value);

					if (dictionary.CheckOption("@MaximumValue2") == true)
					{
						std::string units;
						double value;
						dictionary.ReadMeasure("@MaximumValue2", value, units);

						if (parameter_type_ == PARAMETER_TYPE_PRESSURE ||
							parameter_type_ == PARAMETER_TYPE_TEMPERATURE_PRESSURE)
							value = ReadScalarValue(value, units, PARAMETER_TYPE_PRESSURE);
						else if (parameter_type_ == PARAMETER_TYPE_TIME)
							value = ReadScalarValue(value, units, PARAMETER_TYPE_TIME);

						AddMaximumValue(value);
					}
				}

				// Final Checks
				for (unsigned int j = 0; j < minimum_values_.size(); j++)
					if (minimum_values_(j) >= maximum_values_(j))
						OpenSMOKE::FatalErrorMessage("ParametricAnalysis_Options: Parameter minimum value must be strictly smaller than maximum value");

				// Populating the list of parameters
				number_of_parameters_ = minimum_values_.size();
				list_of_values_.resize(number_of_points_, number_of_parameters_);
				for (unsigned int j = 0; j < number_of_parameters_; j++)
				{
					const double step = (maximum_values_(j) - minimum_values_(j)) / double(number_of_points_ - 1);
					list_of_values_(0,j) = minimum_values_(j);
					for (unsigned int i = 1; i < list_of_values_.rows(); i++)
						list_of_values_(i,j) = list_of_values_(i-1,j) + step;
				}
			}
		}

		if (dictionary.CheckOption("@NumberOfThreads") == true)
		{
			std::string value;
			dictionary.ReadString("@NumberOfThreads", value);
			SetNumberOfThreads(value);
		}

		// Clean
		const bool clean = false;
		if (number_of_parameters_ == 1 && clean == true)
		{
			std::vector<double> temp(number_of_points_);
			for (unsigned int i = 0; i < number_of_points_; i++)
				temp[i] = list_of_values_(i, 0);

			// reorder and remove duplicates
			std::sort(temp.begin(), temp.end());
			std::vector<double>::iterator last = std::unique(temp.begin(), temp.end());
			temp.erase(last, temp.end());

			number_of_points_ = boost::lexical_cast<int>(temp.size());
			list_of_values_.resize(number_of_parameters_, 1);

			for (unsigned int i = 0; i < number_of_points_; i++)
				list_of_values_(i, 0) = temp[i];
		}

		// Combinations
		if (combinations_ == true && number_of_parameters_ == 2)
		{
			Eigen::MatrixXd temp(number_of_points_ * number_of_points_, number_of_parameters_);

			unsigned int count = 0;
			for (unsigned int i = 0; i < number_of_points_; i++)
			{
				for (unsigned int j = 0; j < number_of_points_; j++)
				{
					temp(count, 0) = list_of_values_(i, 0);
					temp(count, 1) = list_of_values_(j, 1);
					count++;
				}
			}

			SetListOfValues(temp);
		}

		// Print on the screen the list of values to be analyzed
		std::cout << std::endl;
		std::cout << "-----------------------------------------------------------" << std::endl;
		std::cout << "Parametric Analysis:" << std::endl;
		std::cout << "-----------------------------------------------------------" << std::endl;
		for (unsigned int i = 0; i < list_of_values_.rows(); i++)
		{
			std::cout << i << " ";
			for (unsigned int j = 0; j < number_of_parameters_; j++)
				std::cout << list_of_values_(i, j) << " ";
			std::cout << std::endl;
		}
		std::cout << "-----------------------------------------------------------" << std::endl;
		std::cout << std::endl;

		enabled_ = true;
	}

	Eigen::VectorXd ReadListOfValues(const std::vector<std::string> values, const ParametricAnalysis_Options::ParametricAnalysisType parameter_type)
	{
		Eigen::VectorXd list_of_values;
		std::string units = values.back();

		if (parameter_type == ParametricAnalysis_Options::PARAMETER_TYPE_TEMPERATURE)
		{
			double conversion_factor = 0.;
			if (values.back() == "K") conversion_factor = 0.;
			else if (values.back() == "C") conversion_factor = 273.15;
			else OpenSMOKE::FatalErrorMessage("ParametricAnalysis_Options: Wrong units for temperature. Allowed units: K | C");

			list_of_values.resize(values.size() - 1);
			for (unsigned int i = 0; i < values.size() - 1; i++)
				list_of_values(i) = (boost::lexical_cast<double>(values[i]) + conversion_factor);
		}
		else if (parameter_type == ParametricAnalysis_Options::PARAMETER_TYPE_PRESSURE)
		{
			double conversion_factor = 1.;
			if (values.back() == "Pa") conversion_factor = 1.;
			else if (values.back() == "atm") conversion_factor = 101325.;
			else if (values.back() == "bar") conversion_factor = 100000.;
			else OpenSMOKE::FatalErrorMessage("ParametricAnalysis_Options: Wrong units for pressure. Allowed units: Pa | atm | bar");

			list_of_values.resize(values.size() - 1);
			for (unsigned int i = 0; i < values.size() - 1; i++)
				list_of_values(i) = (boost::lexical_cast<double>(values[i])*conversion_factor);
		}
		else if (parameter_type == ParametricAnalysis_Options::PARAMETER_TYPE_TIME)
		{
			double conversion_factor = 1.;
			if (values.back() == "s") conversion_factor = 1.;
			else if (values.back() == "ms") conversion_factor = 1e-3;
			else if (values.back() == "min") conversion_factor = 60.;
			else OpenSMOKE::FatalErrorMessage("ParametricAnalysis_Options: Wrong units for time. Allowed units: s | ms | min");

			list_of_values.resize(values.size() - 1);
			for (unsigned int i = 0; i < values.size() - 1; i++)
				list_of_values(i) = (boost::lexical_cast<double>(values[i])*conversion_factor);
		}

		return list_of_values;
	}

	double ReadScalarValue(const double input_value, const std::string units, const ParametricAnalysis_Options::ParametricAnalysisType parameter_type)
	{
		double value = input_value;

		if (parameter_type == ParametricAnalysis_Options::PARAMETER_TYPE_TEMPERATURE)
		{
			if (units == "K")			value += 0.;
			else if (units == "C")		value += 273.15;
			else OpenSMOKE::FatalErrorMessage("ParametricAnalysis_Options: Unknown temperature units. Allowed units: K | C");
		}
		else if (parameter_type == ParametricAnalysis_Options::PARAMETER_TYPE_PRESSURE)
		{
			if (units == "Pa")			value *= 1.;
			else if (units == "atm")	value *= 101325.;
			else if (units == "bar")	value *= 100000.;
			else OpenSMOKE::FatalErrorMessage("ParametricAnalysis_Options: Unknown pressure units. Allowed units: Pa | atm | bar");
		}
		else if (parameter_type == ParametricAnalysis_Options::PARAMETER_TYPE_TIME)
		{
			if (units == "s")			value *= 1.;
			else if (units == "ms")		value *= 1e-3;
			else if (units == "min")	value *= 60.;
			else OpenSMOKE::FatalErrorMessage("ParametricAnalysis_Options: Unknown time units. Allowed units: s | ms | min");
		}

		return value;
	}

}
