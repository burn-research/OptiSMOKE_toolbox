/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alessandro Stagni <alessandro.stagni@polimi.it>               |
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
|   Copyright(C) 2016, 2015 Alessandro Stagni & Alberto Cuoci             |
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

FugacityGIdealGas_mix::FugacityGIdealGas_mix(const unsigned int NS) :
FugacityModel_mix(NS)
{  
}

FugacityGIdealGas_mix::~FugacityGIdealGas_mix() 
{ 
}

std::vector<double> FugacityGIdealGas_mix::fugacity_coefficient(const double T, const double P, const std::vector<double>& x)
{
	std::vector<double> phi(NS_);

	for (unsigned int i = 0; i < NS_; i++)
		phi[i] = 1.;

	return phi;
}

std::vector<double> FugacityGIdealGas_mix::fugacity(const double T, const double P, const std::vector<double>& x)
{
	std::vector<double> phi = fugacity_coefficient(T, P, x);

	std::vector<double> f(NS_);
	for (unsigned int i = 0; i < NS_; i++)
		f[i] = P * x[i] * phi[i];

	return f;
}

std::vector<double> FugacityGIdealGas_mix::fugacity_coefficient_pure(const double T, const double P)
{
	OpenSMOKE::FatalErrorMessage("vector<double> FugacityGIdealGas_mix::fugacity_coefficient_pure(double T, double P)");
	std::vector<double> dummy;
	return dummy;
}
