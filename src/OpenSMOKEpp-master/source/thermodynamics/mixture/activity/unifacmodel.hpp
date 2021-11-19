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

UnifacModel::UnifacModel(	const unsigned int NS,
							const std::vector< std::vector<double> >& R,
							const std::vector< std::vector<double> >& Q,
							const std::vector< std::vector<int> >& subgroup,
							const std::vector< std::vector<int> >& maingroup,
							const std::vector< std::vector<int> >& number,
							const Eigen::MatrixXd& interactiontable) :
ActivityModel(NS),
R_(R),
Q_(Q),
subgroup_(subgroup),
maingroup_(maingroup),
number_(number),
interactiontable_(interactiontable)
{
	lngammaR_.resize(NS_);
	lngammaC_.resize(NS_);
	lngamma_.resize(NS_);
	gamma_.resize(NS_);

	MemoryAllocation();
}

void UnifacModel::MemoryAllocation()
{
	nGroups_.resize(NS_);
	for (unsigned int i = 0; i < NS_; i++)
		nGroups_[i] = R_[i].size();

	nTotalGroups_ = 0;
	for (unsigned int i = 0; i < NS_; i++)
	{
		for (int j = 0; j < nGroups_[i]; j++)
		{
			int dummy = subgroup_[i][j];

			if (!OpenSMOKE::IsValuePresent(dummy, TotalSubGroups_))
			{
				nTotalGroups_++;
				TotalSubGroups_.push_back(subgroup_[i][j]);
				TotalMainGroups_.push_back(maingroup_[i][j]);
				QG_.push_back(Q_[i][j]);
			}
		}
	}

	thetaGm_.resize(nTotalGroups_);
	xgroupG_.resize(nTotalGroups_);
	lnGammaGk_.resize(nTotalGroups_);
	thetaPsiGmk_.resize(nTotalGroups_);
	thetaPsiGkm_.resize(nTotalGroups_);

	thetaPm_.resize(NS_);
	xgroupP_.resize(NS_);
	lnGammaPk_.resize(NS_);
	thetaPsiPmk_.resize(NS_);
	thetaPsiPkm_.resize(NS_);

	for (unsigned int i = 0; i < NS_; i++)
	{
		thetaPm_[i].resize(nGroups_[i], 0.);
		xgroupP_[i].resize(nGroups_[i], 0.);
		lnGammaPk_[i].resize(nGroups_[i], 0.);
		thetaPsiPmk_[i].resize(nGroups_[i], 0.);
		thetaPsiPkm_[i].resize(nGroups_[i], 0.);
	}

	QgroupsumP_.resize(NS_);

	ri_.resize(NS_);
	qi_.resize(NS_);
	li_.resize(NS_);
	theta_.resize(NS_);
	phi_.resize(NS_);
	z = 10;
}

void UnifacModel::UpdateGamma(const double T, const double P, const std::vector<double>& x_original)
{
	std::vector<double> x = x_original;
	OpenSMOKE::CheckAndCorrectSumOfFractions(x);

	UpdateGammaC(T, P, x);
	UpdateGammaR(T, P, x);

	for (unsigned int i = 0; i < NS_; i++)
	{
		lngamma_[i] = lngammaC_[i] + lngammaR_[i];
		gamma_[i] = exp(lngamma_[i]);
	}
}

void UnifacModel::UpdateGammaC(const double T, const double P, const std::vector<double>& x)
{
	Rsum_ = 0.;
	Qsum_ = 0.;
	Lsum_ = 0.;

	memset(&ri_[0], 0, ri_.size() * sizeof ri_[0]);
	memset(&qi_[0], 0, qi_.size() * sizeof qi_[0]);


	for (unsigned int i = 0; i < NS_; i++)
	{
		for (int j = 0; j < nGroups_[i]; j++)
		{
			ri_[i] += R_[i][j] * number_[i][j];
			qi_[i] += Q_[i][j] * number_[i][j];
		}
	}

	for (unsigned int i = 0; i < NS_; i++)
		li_[i] = (z / 2.) * (ri_[i] - qi_[i]) - (ri_[i] - 1.);

	for (unsigned int i = 0; i < NS_; i++)
	{
		Rsum_ += ri_[i] * x[i];
		Qsum_ += qi_[i] * x[i];
		Lsum_ += li_[i] * x[i];
	}

	for (unsigned int i = 0; i < NS_; i++)
	{
		phi_[i] = ri_[i] * x[i] / Rsum_;
		theta_[i] = qi_[i] * x[i] / Qsum_;
	}

	for (unsigned int i = 0; i < NS_; i++)
	{
		lngammaC_[i] =	1. - ri_[i] / Rsum_ + std::log(ri_[i] / Rsum_) -
						5. * qi_[i] * (1. - ri_[i] * Qsum_ / (qi_[i] * Rsum_) +
						std::log(ri_[i] * Qsum_ / (qi_[i] * Rsum_)));
	}

}

void UnifacModel::UpdateGammaR(const double T, const double P, const std::vector<double>& x)
{
	std::memset(&QgroupsumP_[0], 0, QgroupsumP_.size() * sizeof QgroupsumP_[0]);
	std::memset(&thetaPsiGkm_[0], 0, thetaPsiGkm_.size() * sizeof thetaPsiGkm_[0]);
	std::memset(&thetaPsiGmk_[0], 0, thetaPsiGkm_.size() * sizeof thetaPsiGkm_[0]);
	std::memset(&lngammaR_[0], 0, lngammaR_.size() * sizeof lngammaR_[0]);
	std::memset(&xgroupG_[0], 0, xgroupG_.size() * sizeof xgroupG_[0]);
  
	for (unsigned int i = 0; i < NS_; i++)
	{
		std::memset(&thetaPm_[i][0], 0, thetaPm_[i].size() * sizeof thetaPm_[i][0]);
		std::memset(&xgroupP_[i][0], 0, xgroupP_[i].size() * sizeof xgroupP_[i][0]);
		std::memset(&lnGammaPk_[i][0], 0, lnGammaPk_[i].size() * sizeof lnGammaPk_[i][0]);
		std::memset(&thetaPsiPmk_[i][0], 0, thetaPsiPmk_[i].size() * sizeof thetaPsiPmk_[i][0]);
		std::memset(&thetaPsiPkm_[i][0], 0, thetaPsiPkm_[i].size() * sizeof thetaPsiPkm_[i][0]);
	}


	//Pure components
	for (unsigned int i = 0; i < NS_; i++)
	{
		for (int j = 0; j < nGroups_[i]; j++)
		{
			xgroupP_[i][j] = double(number_[i][j]) /
			double(std::accumulate(number_[i].begin(), number_[i].end(), 0.));
		}
	}

	for (unsigned int i = 0; i < NS_; i++)
		for (int j = 0; j < nGroups_[i]; j++)
			QgroupsumP_[i] += xgroupP_[i][j] * Q_[i][j];

	for (unsigned int i = 0; i < NS_; i++)
		for (int j = 0; j < nGroups_[i]; j++)
			thetaPm_[i][j] = Q_[i][j] * xgroupP_[i][j] / QgroupsumP_[i];

	for (unsigned int i = 0; i < NS_; i++)
		for (int j = 0; j < nGroups_[i]; j++)
			for (int k = 0; k < nGroups_[i]; k++)
			{
				thetaPsiPmk_[i][j] +=	(thetaPm_[i][k] *
										std::exp(-interactiontable_(maingroup_[i][k] - 1, maingroup_[i][j] - 1) / T));
			}

	for (unsigned int i = 0; i < NS_; i++)
		for (int j = 0; j < nGroups_[i]; j++)
		{
			for (int k = 0; k < nGroups_[i]; k++)
			{
				double Psinm = 0;
				for (int m = 0; m < nGroups_[i]; m++)
					Psinm +=	(thetaPm_[i][m] *
								std::exp(-interactiontable_(maingroup_[i][m] - 1, maingroup_[i][k] - 1) / T));

				thetaPsiPkm_[i][j] +=	(thetaPm_[i][k] *
										std::exp(-interactiontable_(maingroup_[i][j] - 1, maingroup_[i][k] - 1) / T)) /
										Psinm;
			}
		}

	for (unsigned int i = 0; i < NS_; i++)
		for (int j = 0; j < nGroups_[i]; j++)
		{
			lnGammaPk_[i][j] = Q_[i][j] * (1. - std::log(thetaPsiPmk_[i][j]) - thetaPsiPkm_[i][j]);
		}

	//Mixtures
	//Weighted sum of groups
	double group_sum = 0.;
	for (unsigned int i = 0; i < NS_; i++)
		group_sum += x[i] * std::accumulate(number_[i].begin(), number_[i].end(), 0.);

	for (unsigned int i = 0; i < NS_; i++)
		for (int j = 0; j < nGroups_[i]; j++)
			xgroupG_[OpenSMOKE::Index(int(subgroup_[i][j]), TotalSubGroups_)] += x[i] * number_[i][j];

	for (int i = 0; i < nTotalGroups_; i++)
		xgroupG_[i] /= group_sum;

	QgroupsumG_ = 0.;
	for (int i = 0; i < nTotalGroups_; i++)
		QgroupsumG_ += xgroupG_[i] * QG_[i];

	for (int i = 0; i < nTotalGroups_; i++)
		thetaGm_[i] = xgroupG_[i] * QG_[i] / QgroupsumG_;

	for (int j = 0; j < nTotalGroups_; j++)
		for (int k = 0; k < nTotalGroups_; k++)
		{
			thetaPsiGmk_[j] += (thetaGm_[k] * std::exp(-interactiontable_(TotalMainGroups_[k] - 1, TotalMainGroups_[j] - 1) / T));
		}

	for (int j = 0; j < nTotalGroups_; j++)
		for (int k = 0; k < nTotalGroups_; k++)
		{
			double Psinm = 0;
			for (int m = 0; m < nTotalGroups_; m++)
				Psinm +=	(thetaGm_[m] *
							std::exp(-interactiontable_(TotalMainGroups_[m] - 1, TotalMainGroups_[k] - 1) / T));

			thetaPsiGkm_[j] += (thetaGm_[k] *
								std::exp(-interactiontable_(TotalMainGroups_[j] - 1, TotalMainGroups_[k] - 1) / T)) /
								Psinm;
		}

	for (int j = 0; j < nTotalGroups_; j++)
		lnGammaGk_[j] = QG_[j] * (1. - log(thetaPsiGmk_[j]) - thetaPsiGkm_[j]);


	//Residual activity coefficients
	for (unsigned int i = 0; i < NS_; i++)
		for (int j = 0; j < nTotalGroups_; j++)
			if (OpenSMOKE::IsValuePresent(TotalSubGroups_[j], subgroup_[i]))
			{
				int group_index = OpenSMOKE::Index(TotalSubGroups_[j], subgroup_[i]);
				lngammaR_[i] += number_[i][group_index] * (lnGammaGk_[j] - lnGammaPk_[i][group_index]);
			}

}
