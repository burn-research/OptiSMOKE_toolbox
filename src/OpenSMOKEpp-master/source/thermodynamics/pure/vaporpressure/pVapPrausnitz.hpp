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

PVapPrausnitz::PVapPrausnitz(double* coefficients)
:
VaporPressureModel(coefficients, 0.0)
{
    
}

double PVapPrausnitz::pVap(double T)
{
    return exp(_coefficients[0]+_coefficients[1]/(T+_coefficients[2])+_coefficients[3]*T+_coefficients[4]*log(T)+_coefficients[5]*pow(T,_coefficients[6]))*10e+04;
}

double PVapPrausnitz::dpVapdT(double T)
{
    double dpdT = exp(_coefficients[0]+_coefficients[1]/(T+_coefficients[2])+_coefficients[3]*T+_coefficients[4]*log(T)+_coefficients[5]*pow(T,_coefficients[6]))*10e+04;
    dpdT *= (-_coefficients[1]/pow(T + _coefficients[2], 2.0) + _coefficients[3] + _coefficients[4]/T + _coefficients[6]*_coefficients[5]*pow(T, _coefficients[6] - 1.0));
    return dpdT;
}
