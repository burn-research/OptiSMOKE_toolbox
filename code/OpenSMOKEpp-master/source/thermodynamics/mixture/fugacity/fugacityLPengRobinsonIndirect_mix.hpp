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

FugacityLPengRobinsonIndirect_mix::FugacityLPengRobinsonIndirect_mix(	const unsigned int NS,
																		const std::vector<double>& Tc,
																		const std::vector<double>& Pc,
																		const std::vector<double>& omega,
																		const std::vector<double>& MW,
																		std::vector<DensityModel*>& density,
																		std::vector<VaporPressureModel*>& vaporPressure,
																		ActivityModel& activity) :
FugacityModel_mix(NS),
Tc_(Tc),
Pc_(Pc),
omega_(omega),
MW_(MW),
density_(density),
vaporPressure_(vaporPressure),
activity_(activity)
{
	pr = new pengrobinson_mix(Tc_, Pc_, omega_, MW_);
	species_fugacity_.resize(NS_);

	for (unsigned int i = 0; i < NS_; i++)
		species_fugacity_[i] = new FugacityGPengRobinson(Tc_[i], Pc_[i], omega_[i], MW_[i]);

	iPoyntingCorrection_ = true;
}

FugacityLPengRobinsonIndirect_mix::~FugacityLPengRobinsonIndirect_mix()
{
	delete pr;

	for (unsigned int i = 0; i < NS_; i++)
		delete species_fugacity_[i];
}

std::vector<double> FugacityLPengRobinsonIndirect_mix::fugacity(const double T, const double P, const std::vector<double>& x)
{
	std::vector<double> f(NS_);
	std::vector<double> phi(NS_);

	for (unsigned int i = 0; i < NS_; i++)
		phi[i] = species_fugacity_[i]->fugacity_coefficient(T, vaporPressure_[i]->pVap(T));

	std::vector<double> vm(NS_);
	std::vector<double> poynting(NS_);

	activity_.UpdateGamma(T, P, x);
	
	std::vector<double> gamma = activity_.gamma();

	//Calculating Specific volume
	if (iPoyntingCorrection_ == true)
	{
		for (unsigned int i = 0; i < NS_; i++)
		{
			vm[i] = MW_[i] / density_[i]->rho(T, P);
			poynting[i] = std::exp(vm[i] * (P - vaporPressure_[i]->pVap(T)) /(8.3144621 * T));
		}
	}
	else
	{
		for (unsigned int i = 0; i < NS_; i++)
			poynting[i] = 1.;
	}

	for (unsigned int i = 0; i < NS_; i++)
	{
		f[i] = vaporPressure_[i]->pVap(T) * phi[i] * poynting[i] *
		x[i] * gamma[i];
	}

	return f;
}

std::vector<double> FugacityLPengRobinsonIndirect_mix::fugacity_coefficient_pure(const double T, const double P)
{
	std::vector<double> phi(NS_);

	for (unsigned int i = 0; i < NS_; i++)
	{
		double Zliq;
		pengrobinson species_pr(Tc_[i], Pc_[i], omega_[i], MW_[i]);
		species_pr.Solve(T, P);
		Zliq = species_pr.Zmin();

		phi[i] =	std::exp(Zliq - 1. - (species_pr.A() / (2. * std::sqrt(2.) * species_pr.B()))
					* std::log((Zliq + species_pr.B() * (1. + std::sqrt(2.)))
					/ ((Zliq + species_pr.B() * (1. - std::sqrt(2.)))))
					- std::log(Zliq - species_pr.B()));
	}

	return phi;
}

std::vector<double> FugacityLPengRobinsonIndirect_mix::fugacity_coefficient(const double T, const double P, const std::vector<double>& x)
{
	std::vector<double> phi(NS_), lnphi(NS_);

	double Zliq;
	pr->Solve(T, P, x);
	Zliq = pr->Zmin();

	for (unsigned int i = 0; i < NS_; i++)
	{
		pengrobinson species_pr(Tc_[i], Pc_[i], omega_[i], MW_[i]);
		species_pr.Solve(T, P);

		//First term
		lnphi[i] = (species_pr.B() / pr->B()) * (Zliq - 1.);

		//Second term
		lnphi[i] += ((pr->A() / (2. * std::sqrt(2.) * pr->B())) * (species_pr.B() / pr->B() -
					2. * std::sqrt(species_pr.A() / pr->A())) *
					std::log((Zliq + pr->B() * (1. + sqrt(2.))) / (Zliq + pr->B() * (1. - sqrt(2.)))));

		//Third term
		lnphi[i] -= std::log(Zliq - pr->B());
		phi[i] = std::exp(lnphi[i]);
	}

	lnphi.clear();

	return phi;
}
