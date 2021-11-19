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

FixedProfile::FixedProfile(const unsigned int n, const double* x, const double* y)
{
	n_ = n;
	x_.resize(n_);
	y_.resize(n_);
	for (unsigned int i = 0; i < n_; i++)
	{
		x_(i) = x[i];
		y_(i) = y[i];
	}
}

void FixedProfile::Interpolate(const Eigen::VectorXd& xx, Eigen::VectorXd& yy)
{
	yy.resize(xx.size());

	if (std::fabs(x_(0) - xx(0)) > 1.e-12)
		OpenSMOKE::FatalErrorMessage("Interpolating fixed profile: the requested coordinate is smaller than the minimum available coordinate");

	if (std::fabs(xx(xx.size() - 1) - x_(n_ - 1)) > 1.e-12)
		OpenSMOKE::FatalErrorMessage("Interpolating fixed profile: the requested coordinate is larger than the maximum available coordinate");

	yy(0) = y_(0);
	yy(yy.size() - 1) = y_(y_.size() - 1);

	for (unsigned int i = 1; i < xx.size() - 1; i++)
	{
		for (unsigned int j = 1; j < n_; j++)
		{
			if (xx(i) <= x_(j))
			{
				yy(i) = y_(j - 1) + (y_(j) - y_(j - 1)) / (x_(j) - x_(j - 1)) * (xx(i) - x_(j - 1));
				break;
			}
		}
	}
}