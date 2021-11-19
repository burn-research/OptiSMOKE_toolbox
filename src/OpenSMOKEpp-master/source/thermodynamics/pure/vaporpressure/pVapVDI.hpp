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

PVapVDI::PVapVDI(double* coefficients, double Tc, double pc)
:
VaporPressureModel(coefficients, Tc),
_pc(pc)
{
  
}

double PVapVDI::pVap(double T)
{
  double TTc = min(T, _Tc) / _Tc;
  return _pc * exp(1.0 / TTc * (_coefficients[0]*(1.0 - TTc) + _coefficients[1] * pow(1.0 - TTc, 1.5) + _coefficients[2] * pow(1.0 - TTc, 3.0) + _coefficients[3] * pow(1.0 - TTc, 6.0)));
}

double PVapVDI::dpVapdT(double T)
{
  double TTc = min(T, _Tc) / _Tc;
  double dpdT = _pc * exp(1.0 / TTc * (_coefficients[0]*(1.0 - TTc) + _coefficients[1] * pow(1.0 - TTc, 1.5) + _coefficients[2] * pow(1.0 - TTc, 3.0) + _coefficients[3] * pow(1.0 - TTc, 6.0)));
  dpdT *= 1.0 / TTc * (
          -_coefficients[0] / T
          - _coefficients[1]*(
          pow(1.0 - TTc, 1.5) / T
          + pow(1.0 - TTc, 0.5)*1.5 / _Tc
          )
          - _coefficients[2]*(
          pow(1.0 - TTc, 3.0) / T
          + pow(1.0 - TTc, 2.0)*3.0 / _Tc
          )
          - _coefficients[3]*(
          pow(1.0 - TTc, 6.0) / T
          + pow(1.0 - TTc, 5.0)*6.0 / _Tc
          )
          );
  return dpdT;
}
