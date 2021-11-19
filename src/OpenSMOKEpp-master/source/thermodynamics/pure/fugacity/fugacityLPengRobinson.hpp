/*-----------------------------------------------------------------------*\
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
|	License                                                           |
|                                                                         |
|   Copyright(C) 2016, 2015 Alberto Cuoci                                 |
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

FugacityLPengRobinson::FugacityLPengRobinson
(const double Tc,
        const double Pc,
        const double omega,
        const double MW) :
FugacityModel(),
Tc_(Tc),
Pc_(Pc),
omega_(omega),
MW_(MW)
{
  pr_ = new pengrobinson(Tc_, Pc_, omega_, MW_);
}

FugacityLPengRobinson::~FugacityLPengRobinson()
{
  delete pr_;
}

double FugacityLPengRobinson::fugacity_coefficient(double T, double P)
{
  double phi;
  double Zliq;
  pr_->Solve(T, P);
  Zliq = pr_->Zmin();

  phi = exp(Zliq - 1 - (pr_->A() / (2 * sqrt(2) * pr_->B()))
          * log((Zliq + pr_->B() * (1 + sqrt(2)))
          / ((Zliq + pr_->B() * (1 - sqrt(2)))))
          - log(Zliq - pr_->B()));
  
  return phi;

}
