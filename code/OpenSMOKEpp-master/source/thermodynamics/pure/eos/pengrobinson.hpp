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

pengrobinson::pengrobinson(const double Tc,
        const double Pc, const double omega,
        const double MW) :
EosModel(Tc, Pc, omega, MW)
{
  Tc_ = Tc;
  Pc_ = Pc;
  omega_ = omega;
  MW_ = MW;

  //Setup EOS parameters
  R = 8.314;
}

pengrobinson::~pengrobinson() { }

void
pengrobinson::Solve(const double T, const double P)
{
  Tr = T / Tc_;
  Pr = P / Pc_;

  m = 0.37464 + 1.54226 * omega_ - 0.26992 * omega_*omega_;

  alfa = std::pow(1 + m * (1 - std::sqrt(Tr)), 2);
  a_ = 0.45724 * (std::pow((R * Tc_), 2.) / Pc_) * alfa;
  b_ = 0.0778 * R * Tc_ / Pc_;

  A_ = a_ * P / std::pow((R * T), 2.);
  B_ = b_ * P / (R * T);

  a1 = -(1 - B_);
  a2 = A_ - 3 * B_ * B_ - 2 * B_;
  a3 = -(A_ * B_ - B_ * B_ - B_ * B_ * B_);

  ZR_ = OpenSMOKE::CubicRootsReal(1., a1, a2, a3);
}
