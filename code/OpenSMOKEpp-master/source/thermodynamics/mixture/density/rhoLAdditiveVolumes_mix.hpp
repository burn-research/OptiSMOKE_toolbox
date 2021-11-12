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

rhoLAdditiveVolumes_mix::rhoLAdditiveVolumes_mix(unsigned int NS,
						 std::vector<DensityModel*>& rhoL,
						 const std::vector<double>& MW) :
densityLModel_mix(NS),
rhoL_(rhoL),
MW_(MW)
{
	y_.resize(NS);
}

rhoLAdditiveVolumes_mix::~rhoLAdditiveVolumes_mix()
{

}

double rhoLAdditiveVolumes_mix::rho(const double T, const double P, const std::vector<double>& x)
{

	//Mixture molecular weight
	double MWmix = 0.;
	for (unsigned int i = 0; i < NS_; i++)
		MWmix += x[i] * MW_[i];

	for (unsigned int i = 0; i < NS_; i++)
		y_[i] = x[i] * MW_[i] / MWmix;

	double urhoL_mix = 0.;

	for (unsigned int i = 0; i < NS_; i++)
	{
		double rhoLspecies = rhoL_[i]->rho(T,P);
		urhoL_mix += y_[i]/rhoLspecies;
	}

	return 1/urhoL_mix;
}
