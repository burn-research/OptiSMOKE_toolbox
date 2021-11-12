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

#ifndef OpenSMOKE_Grammar_XYProfile_H
#define	OpenSMOKE_Grammar_XYProfile_H

#include <string>
#include "boost/filesystem.hpp"
#include "math/OpenSMOKEVector.h"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
	class Grammar_XYProfile : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@XUnits", 
																OpenSMOKE::SINGLE_STRING, 
																"Units of x variable", 
																true) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@YUnits", 
																OpenSMOKE::SINGLE_STRING, 
																"Units of y variable", 
																true) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@XVariable", 
																OpenSMOKE::SINGLE_STRING, 
																"X variable", 
																true) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@YVariable", 
																OpenSMOKE::SINGLE_STRING, 
																"X variable", 
																true) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Profile", 
																OpenSMOKE::VECTOR_DOUBLE, 
																"X,Y profile", 
																true) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FilterWidth", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Width of filter (same units of @XUnits)", 
																false,
																"none",
																"none",
																"none") );
		}
	};

	enum ConversionType { CONVERSIONTYPE_X , CONVERSIONTYPE_Y} ;

	double Conversion(const ConversionType index, const std::string variable, const std::string units)
	{
		double conversion = 1.;

		if (index == ConversionType::CONVERSIONTYPE_X)
		{
			if (variable == "length")
			{
				if (units == "m")			conversion = 1.;
				else if (units == "cm")		conversion = 0.01;
				else if (units == "mm")		conversion = 0.001;
				else OpenSMOKE::FatalErrorMessage("Unknown length unit. Allowed units: m | cm | mm");
			}
			else if (variable == "time")
			{
				if (units == "s")			conversion = 1.;
				else if (units == "min")	conversion = 60.;
				else if (units == "h")		conversion = 3600.;
				else if (units == "ms")		conversion = 0.001;
				else OpenSMOKE::FatalErrorMessage("Unknown time unit. Allowed units: s | min | h | ms");
			}
			else OpenSMOKE::FatalErrorMessage("Unknown x variable for xy profile");

			return conversion;
		}
		else if (index == ConversionType::CONVERSIONTYPE_Y)
		{
			if (variable == "temperature")
			{
				if (units == "K")			conversion = 1.;
				else OpenSMOKE::FatalErrorMessage("Unknown temperature units. Allowed units: K");
			}
			else if (variable == "pressure")
			{
				if (units == "Pa")			conversion = 1.;
				else if (units == "atm")	conversion = 101325.;
				else if (units == "bar")	conversion = 100000.;
				else OpenSMOKE::FatalErrorMessage("Unknown pressure units. Allowed units: Pa | atm | bar");
			}
			else if (variable == "volume")
			{
				if (units == "m3")			conversion = 1.;
				else if (units == "l")	conversion = 1.e-3;
				else if (units == "dm3")	conversion = 1.e-3;
				else if (units == "cm3")	conversion = 1e-6;
				else if (units == "mm3")	conversion = 1e-9;
				else OpenSMOKE::FatalErrorMessage("Unknown volume units. Allowed units: m3 | cm3 | mm3 | l | dm3");
			}
			else if (variable == "specific-mass-flow-rate")
			{
				if (units == "kg/m2/s")			conversion = 1.;
				else if (units == "g/cm2/s")	conversion = 10.;
				else OpenSMOKE::FatalErrorMessage("Unknown specific-mass-flow-rate units. Allowed units: kg/m2/s | g/cm2/s");
			}
			else if (variable == "area-per-unit-of-volume")
			{
				if (units == "1/m")				conversion = 1.;
				else if (units == "1/cm")		conversion = 1.e2;
				else if (units == "1/mm")		conversion = 1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown area-per-unit-of-volume units. Allowed units: 1/m | 1/cm | 1/mm");
			}
			else if (variable == "dimensionless")
			{
				if (units == "dimensionless")	conversion = 1.;
				else if (units == "none")		conversion = 1.;
				else if (units == "-")			conversion = 1.;
				else OpenSMOKE::FatalErrorMessage("Unknown dimensionless units. Allowed units: dimensionless | none | -");
			}

			else OpenSMOKE::FatalErrorMessage("Unknown y variable for xy profile");
		}

		return conversion;
	}

	void GetXYProfileFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, 
									OpenSMOKE::OpenSMOKEVectorDouble& x, OpenSMOKE::OpenSMOKEVectorDouble& y, 
									std::string& x_variable, std::string& y_variable)
	{
		Grammar_XYProfile grammar_xyprofile_status;
		dictionary.SetGrammar(grammar_xyprofile_status);
		
		dictionary.ReadString("@XVariable", x_variable);
		dictionary.ReadString("@YVariable", y_variable);

		std::string x_units;
		dictionary.ReadString("@XUnits", x_units);
		std::string y_units;
		dictionary.ReadString("@YUnits", y_units);

		std::vector<double> xy_profile;
		dictionary.ReadOption("@Profile", xy_profile);

		const double conversion_x = Conversion(ConversionType::CONVERSIONTYPE_X, x_variable, x_units);
		const double conversion_y = Conversion(ConversionType::CONVERSIONTYPE_Y, y_variable, y_units);

		const int n = boost::lexical_cast<int>(xy_profile.size()) / 2;
		ChangeDimensions(n, &x, true);
		ChangeDimensions(n, &y, true);
		
		unsigned int count = 0;
		for(int i=1;i<=n;i++)
		{
			x[i] = xy_profile[count++] * conversion_x;
			y[i] = xy_profile[count++] * conversion_y;
		}

		// Filtering
		if (dictionary.CheckOption("@FilterWidth") == true)
		{
			double width = 0.;
			std::string width_units = "dimensionless";
			dictionary.ReadMeasure("@FilterWidth", width, width_units);
			width *= Conversion(ConversionType::CONVERSIONTYPE_X, x_variable, width_units);

			OpenSMOKE::OpenSMOKEVectorDouble x_new;
			OpenSMOKE::OpenSMOKEVectorDouble y_new;
			ApplyFilter(width, 5, x_new, y_new, x, y);
			x = x_new;
			y = y_new;
		}
	}

	void GetXYProfileFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary,
		Eigen::VectorXd& x, Eigen::VectorXd& y,
		std::string& x_variable, std::string& y_variable)
	{
		Grammar_XYProfile grammar_xyprofile_status;
		dictionary.SetGrammar(grammar_xyprofile_status);

		dictionary.ReadString("@XVariable", x_variable);
		dictionary.ReadString("@YVariable", y_variable);

		std::string x_units;
		dictionary.ReadString("@XUnits", x_units);
		std::string y_units;
		dictionary.ReadString("@YUnits", y_units);

		std::vector<double> xy_profile;
		dictionary.ReadOption("@Profile", xy_profile);

		const double conversion_x = Conversion(ConversionType::CONVERSIONTYPE_X, x_variable, x_units);
		const double conversion_y = Conversion(ConversionType::CONVERSIONTYPE_Y, y_variable, y_units);

		const int n = boost::lexical_cast<int>(xy_profile.size()) / 2;
		x.resize(n);
		y.resize(n);

		unsigned int count = 0;
		for (int i = 1; i <= n; i++)
		{
			x(i-1) = xy_profile[count++] * conversion_x;
			y(i-1) = xy_profile[count++] * conversion_y;
		}

		// Filtering
		if (dictionary.CheckOption("@FilterWidth") == true)
		{
			double width = 0.;
			std::string width_units = "dimensionless";
			dictionary.ReadMeasure("@FilterWidth", width, width_units);
			width *= Conversion(ConversionType::CONVERSIONTYPE_X, x_variable, width_units);

			Eigen::VectorXd x_new;
			Eigen::VectorXd y_new;
			ApplyFilter(width, 5, x.size(), x_new.data(), y_new.data(), x.data(), y.data());
			x = x_new;
			y = y_new;
		}
	}
}

#endif	/* OpenSMOKE_Grammar_XYProfile_H */

