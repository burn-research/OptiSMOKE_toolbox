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

FugacityGPengRobinson_mix::FugacityGPengRobinson_mix(	const unsigned int NS,
														const std::vector<double>& Tc,
														const std::vector<double>& Pc,
														const std::vector<double>& omega,
														const std::vector<double>& MW) :
FugacityModel_mix(NS),
Tc_(Tc),
Pc_(Pc),
omega_(omega),
MW_(MW)
{
	pr = new pengrobinson_mix(Tc, Pc, omega, MW);
}

FugacityGPengRobinson_mix::~FugacityGPengRobinson_mix()
{
	delete pr;
}

std::vector<double> FugacityGPengRobinson_mix::fugacity_coefficient(const double T, const double P, const std::vector<double>& x)
{
	std::vector<double> phi(NS_);
	std::vector<double> lnphi(NS_);

	double Zgas;
	pr->Solve(T, P, x);
	Zgas = pr->Zmax();

	for (unsigned int i = 0; i < NS_; i++)
	{
		pengrobinson species_pr(Tc_[i], Pc_[i], omega_[i], MW_[i]);
		species_pr.Solve(T, P);

		//First term
		lnphi[i] = (species_pr.B() / pr->B()) * (Zgas - 1.);

		//Second term
		lnphi[i] += ((pr->A() / (2. * std::sqrt(2.) * pr->B())) * (species_pr.B() / pr->B() - 2.*std::sqrt(species_pr.A() / pr->A())) *
					std::log((Zgas + pr->B() * (1. + std::sqrt(2.))) / (Zgas + pr->B() * (1.-std::sqrt(2.)))));

		//Third term
		lnphi[i] -= std::log(Zgas - pr->B());
		phi[i] = std::exp(lnphi[i]);
	}

	lnphi.clear();

	return phi;
}

std::vector<double> FugacityGPengRobinson_mix::fugacity(const double T, const double P, const std::vector<double>& x)
{
	std::vector<double> phi = fugacity_coefficient(T, P, x);

	std::vector<double> f(NS_);
	for (unsigned int i = 0; i < NS_; i++)
		f[i] = P * x[i] * phi[i];

	return f;
}

std::vector<double> FugacityGPengRobinson_mix::fugacity_coefficient_pure(const double T, const double P) 
{ 
	OpenSMOKE::FatalErrorMessage("std::vector<double> FugacityGPengRobinson_mix::fugacity_coefficient_pure(double T, double P)");
	std::vector<double> dummy;
	return dummy;
}
