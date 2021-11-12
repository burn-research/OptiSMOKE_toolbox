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

stefanmaxwell::stefanmaxwell(	DinfModel& dinf_model,
								ActivityModel& activity_model) :
dinf_model_(dinf_model),
activity_model_(activity_model)
{
	NS_ = dinf_model.NS();
	D_sm_.resize(NS_, NS_);
	B_.resize(NS_ - 1, NS_ - 1);
	tm_.resize(NS_ - 1, NS_ - 1);

	A_.resize(NS_ - 1, NS_ - 1);
	b_.resize(NS_ - 1);

	Dinf_ = &dinf_model.Dinf();

	last_index_ = NS_ - 1;

	for (unsigned int i = 0; i < NS_; i++)
		for (unsigned int j = 0; j < NS_; j++)
			D_sm_(i, j) = 0.;
}

void stefanmaxwell::UpdateBinaryDiffusivities(const std::vector<double>& x_original)
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	for (unsigned int i = 0; i < NS_; i++)
		for (unsigned int j = 0; j < NS_; j++)
		{
			D_sm_(i, j) =	std::pow(Dinf_->coeff(i, j),
							(1. + x[j] - x[i]) / 2.) *
							std::pow(Dinf_->coeff(j, i),
							(1. + x[i] - x[j]) / 2.);
		}
}

void stefanmaxwell::UpdateB(const double T, const double P, const std::vector<double>& x_original)
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	for (unsigned int i = 0; i < NS_ - 1; i++)
		for (unsigned int j = 0; j < NS_ - 1; j++)
			B_(i, j) = 0.;

	//Diagonal elements
	int p, q;
	for (unsigned int i = 0; i < NS_; i++)
	{
		if (i < last_index_)
			p = i;
		else if (i == last_index_)
			continue;
		else if (i > last_index_)
			p = i - 1;

		B_(p, p) = x[i] / D_sm_(i, last_index_);
		
		for (unsigned int j = 0; j < NS_; j++)
			if (j != i)
				B_(p, p) += x[j] / D_sm_(i, j);
	}

	//Non-diagonal elements
	for (unsigned int i = 0; i < NS_; i++)
		for (unsigned int k = 0; k < NS_; k++)
		{
			if (i < last_index_)
				p = i;
			else if (i == last_index_)
				continue;
			else if (i > last_index_)
				p = i - 1;

			if (k < last_index_)
				q = k;
			else if (k == last_index_)
				continue;
			else if (k > last_index_)
				q = k - 1;

			if (p != q) // k != i or p != q??
				B_(p, q) = -x[i] * (1./D_sm_(i, k) - 1./D_sm_(i, last_index_));
		}

}

void stefanmaxwell::UpdateThermodynamicMatrix(const double T, const double P, const std::vector<double>& x_original)
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	double eps;

	for (unsigned int i = 0; i < NS_ - 1; i++)
		for (unsigned int j = 0; j < NS_ - 1; j++)
			tm_(i, j) = 0.;


	for (unsigned int i = 0; i < NS_ - 1; i++)
		tm_(i, i) += 1.;

	//Add the activity part
	int p, q;
	for (unsigned int i = 0; i < NS_; i++)
		for (unsigned int k = 0; k < NS_; k++)
		{
			eps = 0.001 * x[k];
			std::vector<double> x_plus, x_minus;
			std::vector<double> lngamma_plus(NS_), lngamma_minus(NS_);

			x_plus = x;
			x_minus = x;

			if (i < last_index_)
				p = i;
			else if (i == last_index_)
				continue;
			else if (i > last_index_)
				p = i - 1;

			if (k < last_index_)
				q = k;
			else if (k == last_index_)
				continue;
			else if (k > last_index_)
				q = k - 1;

			double eps_plus = std::min(eps, 1. - x[k]);
			eps_plus = std::min(eps_plus, x[last_index_]);
			double eps_minus = std::min(eps, x[k]);
			eps_minus = std::min(eps_minus, 1. - x[last_index_]);

			x_plus[k] += eps_plus;
			x_minus[k] -= eps_minus;

			x_plus[last_index_] -= eps_plus;
			x_minus[last_index_] += eps_minus;

			if (x_plus[k] - x_minus[k] == 0.)
				continue;

			activity_model_.UpdateGamma(T, P, x_plus);
			lngamma_plus = activity_model_.lngamma();

			activity_model_.UpdateGamma(T, P, x_minus);
			lngamma_minus = activity_model_.lngamma();

			tm_(p, q) += x[i] * (lngamma_plus[i] - lngamma_minus[i]) /
						 (x_plus[k] - x_minus[k]);
		}
}

void stefanmaxwell::Update(const double T, const double P, const std::vector<double>& x)
{
	UpdateBinaryDiffusivities(x);
	UpdateB(T, P, x);
	UpdateThermodynamicMatrix(T, P, x);
}

Eigen::VectorXd stefanmaxwell::molar_diffusion_fluxes(const double ctot, const std::vector<double>& dx)
{
	Eigen::VectorXd x(NS_ - 1);
	Eigen::VectorXd dX_dx(NS_ - 1);

	for (unsigned int i = 0; i < NS_ - 1; i++)
		dX_dx(i) = dx[i];


	A_ = -B_;
	b_ = tm_ * dX_dx;
	//  b_ = dX_dx;
	b_ *= ctot;

	x = A_.lu().solve(b_);

	return x;
}

Eigen::VectorXd stefanmaxwell::molar_diffusion_fluxes(const double T, const double P, const std::vector<double>& x, const double ctot, const std::vector<double>& dx)
{
	Update(T, P, x);
	return molar_diffusion_fluxes(ctot, dx);
}
