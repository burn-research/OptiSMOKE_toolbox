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

/**
 \class UnifacModel
 \brief Class managing UNIFAC activity model
 
 Class managing UNIFAC activity model
 */

#ifndef OPENSMOKEPP_UNIFACMODEL_H
#define OPENSMOKEPP_UNIFACMODEL_H

#include "activitymodel.h"

class UnifacModel : public ActivityModel
{
	public:

		UnifacModel(	const unsigned int NS,
						const std::vector< std::vector<double> >& R,
						const std::vector< std::vector<double> >& Q,
						const std::vector< std::vector<int> >& subgroup,
						const std::vector< std::vector<int> >& maingroup,
						const std::vector< std::vector<int> >& number,
						const Eigen::MatrixXd& interactiontable);
  
		void UpdateGamma(const double T, const double P, const std::vector<double>& x);

	private:

		void UpdateGammaC(const double T, const double P, const std::vector<double>& x);
		void UpdateGammaR(const double T, const double P, const std::vector<double>& x);

		void MemoryAllocation();

	private:

		const std::vector< std::vector<double> >& R_;
		const std::vector< std::vector<double> >& Q_;
		const std::vector< std::vector<int> >& subgroup_;
		const std::vector< std::vector<int> >& maingroup_;
		const std::vector< std::vector<int> >& number_;
		const Eigen::MatrixXd& interactiontable_;

		//Service variables
		std::vector<int> nGroups_;
  
		//Combinatorial part
		std::vector<double> ri_, qi_, li_;
		double Rsum_, Qsum_, Lsum_;
		std::vector<double> theta_, phi_;
		double z;
  
		std::vector<int> TotalSubGroups_;
		std::vector<int> TotalMainGroups_;
		std::vector<double> QG_;
		int nTotalGroups_;
  
		std::vector<double> thetaGm_;
		std::vector<double> xgroupG_;
		std::vector<double> lnGammaGk_;
		std::vector<double> thetaPsiGmk_;
		std::vector<double> thetaPsiGkm_;
		double QgroupsumG_;
  
		std::vector<double> QgroupsumP_;
		std::vector< std::vector<double> > thetaPm_;
		std::vector< std::vector<double> > xgroupP_;
		std::vector< std::vector<double> > lnGammaPk_;
		std::vector< std::vector<double> > thetaPsiPmk_;
		std::vector< std::vector<double> > thetaPsiPkm_;
};

#include "unifacmodel.hpp"

#endif	// OPENSMOKEPP_UNIFACMODEL_H
