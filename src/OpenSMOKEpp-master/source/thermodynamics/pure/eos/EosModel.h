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

/**
 \class EoSModel
 \brief Abstract base class to define eos properties for cubic equations

 */
#include <iostream>

#ifndef _EOSMODEL_H_
#define _EOSMODEL_H_


using namespace std;

class EosModel {
public:
    EosModel(const double Tc,
                const double Pc, const double omega,
                const double MW);
  
  inline double a() const                                 {return a_;};
  inline double b() const                                 {return b_;};
  inline double A() const                                 {return A_;};
  inline double B() const                                 {return B_;};
  inline const  vector<double>& ZR() const                       {return ZR_;};
  
  double Zmin() const;
  double Zmax() const;
  
private:
  friend class pengrobinson;
    
  double Tc_;
  double Pc_;
  double omega_;
  double MW_;

  double R; // gas constant [J/(mol K)]

  vector<double> ZR_;
  double Tr; //Reduced temperature
  double Pr; //Reduced pressure

  double m;
  double alfa;
  double a_;
  double b_;
  double A_;
  double B_;
  
  double a1,a2,a3;
  
  virtual void Solve(const double T, const double P)                {};
};

#include "EosModel.hpp"

#endif
