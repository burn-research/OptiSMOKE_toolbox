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

namespace OpenSMOKE
{
	BatchReactor_VolumeProfile::BatchReactor_VolumeProfile()
	{
		type_ = BatchReactor_VolumeProfile_PressureCoefficient;
		pressureCoefficient_ = 0;
	}
	
	void BatchReactor_VolumeProfile::SetProfile(const OpenSMOKE::OpenSMOKEVectorDouble& x, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		if (x[1] !=0 )
			OpenSMOKE::FatalErrorMessage("The user-defined profile must start from time equal to 0");
		for(int i=2;i<=x.Size();i++)
			if (x[i] <= x[i-1])
				OpenSMOKE::FatalErrorMessage("The sequence of times in the profile must monotonically increase");
		
		n_ = x.Size();
		x_ = x;
		y_ = y;
		x0_ = x_[1];
		xf_ = x_[n_];
		y0_ = y_[1];
		yf_ = y_[n_];
		
		ChangeDimensions(n_-1, &m_, true);

		for(unsigned int i=2;i<=n_;i++)
			m_[i-1] = (y_[i]-y_[i-1])/(x_[i]-x_[i-1]);

		type_ = BatchReactor_VolumeProfile_VolumeHistory;
	}

	void BatchReactor_VolumeProfile::SetPressureProfile(const OpenSMOKE::OpenSMOKEVectorDouble& x, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		SetProfile(x, y);
		type_ = BatchReactor_VolumeProfile_FromPressureProfile;
	}

	void BatchReactor_VolumeProfile::SetPressureCoefficient(const double pressureCoefficient)
	{
		type_ = BatchReactor_VolumeProfile_PressureCoefficient;
		pressureCoefficient_ = pressureCoefficient;
	}
			
	double BatchReactor_VolumeProfile::operator () (const double x) const
	{
		if (x<x_[1])
			OpenSMOKE::FatalErrorMessage("A request for a negative time was performed");
		
		if (x>=x_[n_])
			return yf_;
			
		if (n_ == 2)
			return y0_ + m_[1]*(x-x0_);

		for(unsigned int i=2;i<=n_;i++)
			if (x <= x_[i])
				return y_[i-1] + m_[i-1]*(x-x_[i-1]);

		return 0;
	}
}

