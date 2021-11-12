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

DinfSL::DinfSL(const unsigned int NS, const std::vector<double>& MW,
	      const std::vector<double>& Tnbp, std::vector<DensityModel*>& rhoL) :
MW_(MW),
Tnbp_(Tnbp),
rhoL_(rhoL),
DinfModel(NS)
{
	Dinf_.resize(NS_, NS_);

	for (unsigned int i = 0; i < NS_; i++)
		for (unsigned int j = 0; j < NS_; j++)
			Dinf_(i, j) = 0.;

}

void DinfSL::UpdateInfiniteDiffusivities(const std::vector<double>& etaL, const double T, const double P)
{
	for (unsigned int i = 0; i < NS_; i++)
	{
		//Calculate i-th molar volume [cm3/mol]
		double ith_mv = MW_[i]/rhoL_[i]->rho(Tnbp_[i], 101325.) * 1.e6;
		for (unsigned int j = 0; j < NS_; j++)
		{
			//Calculate j-th molar volume [cm3/mol]
			double jth_mv = MW_[j]/rhoL_[j]->rho(Tnbp_[j], 101325.) * 1.e6;
			//Siddiqi - Lucas (1986)
			Dinf_(i, j) =	correction_ * 9.89e-8 * std::pow(etaL[j]*1.e3, -0.907) *
							std::pow(ith_mv, -0.45) *
							std::pow(jth_mv, 0.265) *
							T * 1.e-4;
		}
	}
}
