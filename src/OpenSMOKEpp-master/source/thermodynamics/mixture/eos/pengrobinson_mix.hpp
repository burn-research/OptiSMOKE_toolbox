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

pengrobinson_mix::pengrobinson_mix(	const std::vector<double>& Tc,
									const std::vector<double>& Pc, 
									const std::vector<double>& omega,
									const std::vector<double>& MW) :
EosModel_mix(Tc, Pc, omega, MW)
{
	//Setup EOS parameters
	R = PhysicalConstants::R_J_mol;

	a_species_.resize(NS_);
	b_species_.resize(NS_);
	A_species_.resize(NS_);
	B_species_.resize(NS_);
}

pengrobinson_mix::~pengrobinson_mix() 
{
}

void pengrobinson_mix::Solve(const double T, const double P, const std::vector<double>& x)
{
	//Setup EOS parameters
	for (unsigned int i = 0; i < NS_; i++)
	{
		const double Tr = T / Tc_[i];
		const double Pr = P / Pc_[i];

		const double m = 0.37464 + 1.54226 * omega_[i] - 0.26992 * omega_[i] * omega_[i];

		const double alfa = std::pow(1. + m * (1. - std::sqrt(Tr)), 2.);
		a_species_[i] = 0.45724 * (std::pow((R * Tc_[i]), 2.) / Pc_[i]) * alfa;
		b_species_[i] = 0.0778 * R * Tc_[i] / Pc_[i];
		A_species_[i] = a_species_[i] * P / std::pow((R * T), 2.);
		B_species_[i] = b_species_[i] * P / (R * T);
	}

	a_ = a_mix(x);
	b_ = b_mix(x);

	A_ = a_ * P / std::pow((R * T), 2.);
	B_ = b_ * P / (R * T);

	double a1 = -(1. - B_);
	double a2 = A_ - 3. * B_ * B_ - 2. * B_;
	double a3 = -(A_ * B_ - B_ * B_ - B_ * B_ * B_);

	ZR_ = OpenSMOKE::CubicRootsReal(1., a1, a2, a3);
}
