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

#include "math/OpenSMOKEVector.h"

namespace OpenSMOKE
{
	PlugFlowReactor_Profile::PlugFlowReactor_Profile(OpenSMOKE::OpenSMOKEVectorDouble& x, OpenSMOKE::OpenSMOKEVectorDouble& y, const std::string x_variable)
	{
		if (x_variable == "length")	time_independent_ = false;
		if (x_variable == "time")	time_independent_ = true;

		x_ = x;
		y_ = y;
		number_of_points_ = x_.Size();

		x0_ = x[1];
		y0_ = y[1];
		xf_ = x[number_of_points_];
		yf_ = y[number_of_points_];

		ChangeDimensions(number_of_points_-1, &m_, true);

		for(unsigned int i=2;i<=number_of_points_;i++)
			m_[i-1] = (y_[i]-y_[i-1])/(x_[i]-x_[i-1]);
	}

	double PlugFlowReactor_Profile::Get(const double x) const
	{
		if (x0_ > 0 && ((x - x0_) / x0_ > -1.e6))
			OpenSMOKE::FatalErrorMessage("Profile class: the required point is outside the domain (too small)");
        //else if ((x - xf_) / xf_ > 1.e-6)
        //   OpenSMOKE::FatalErrorMessage("Profile class: the required point is outside the domain. Please, increase the interval where the profile is defined");

		if (number_of_points_ == 2)
			return y0_ + m_[1]*(x-x0_);

		for(unsigned int i=2;i<=number_of_points_;i++)
			if (x <= x_[i])
				return y_[i-1] + m_[i-1]*(x-x_[i-1]);

		// In case of small excess
		return y_[number_of_points_];
	}
}

